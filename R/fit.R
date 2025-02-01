#' @rdname sequenza
#' @export
sequenza.fit <- function(sequenza.extract, female = TRUE, N.ratio.filter = 10,
    N.BAF.filter = 1, segment.filter = 3e+06, mufreq.threshold = 0.1,
    XY = c(X = "X", Y = "Y"), cellularity = seq(0.1, 1, 0.01),
    ploidy = seq(1, 7, 0.1), ratio.priority = FALSE, method = "baf",
    priors.table = data.frame(CN = 2, value = 2), chromosome.list = 1:24,
    mc.cores = getOption("mc.cores", 2L), verbose = TRUE) {
    # Validate input
    validate_fit_input(sequenza.extract, method, ploidy, cellularity)

    # Add progress tracking
    if (verbose)
        message("Starting sequenza fit analysis...")

    # Prepare input data with progress tracking
    data <- tryCatch({
        prepare_fit_data(sequenza.extract, chromosome.list)
    }, error = function(e) {
        stop("Error preparing data: ", e$message)
    })

    avg.depth.ratio <- sequenza.extract$avg.depth.ratio

    # Process method with better error handling
    result <- tryCatch({
        if (method == "baf") {
            process_baf_method(data$segs, data$segs_len, filters = list(N.ratio = N.ratio.filter,
                N.BAF = N.BAF.filter, segment = segment.filter),
                params = list(female = female, XY = XY, cellularity = cellularity,
                  ploidy = ploidy, ratio.priority = ratio.priority,
                  priors.table = priors.table, mc.cores = mc.cores,
                  verbose = verbose), avg.depth.ratio = avg.depth.ratio)
        } else if (method == "mufreq") {
            process_mufreq_method(data$mutations, params = list(female = female,
                XY = XY, cellularity = cellularity, ploidy = ploidy,
                threshold = mufreq.threshold, priors.table = priors.table,
                mc.cores = mc.cores), avg.depth.ratio = avg.depth.ratio)
        }
    }, error = function(e) {
        stop("Error in ", method, " method: ", e$message)
    }, finally = {
        if (verbose)
            message("Analysis completed")
    })

    # Add metadata to results
    result$params <- list(method = method, female = female, filters = list(N.ratio = N.ratio.filter,
        N.BAF = N.BAF.filter, segment = segment.filter), chromosomes = chromosome.list)

    return(result)
}

# Add input validation function
validate_fit_input <- function(sequenza.extract, method, ploidy,
    cellularity) {
    if (!is.list(sequenza.extract)) {
        stop("sequenza.extract must be a list object")
    }

    if (!method %in% c("baf", "mufreq")) {
        stop("method must be either 'baf' or 'mufreq'")
    }

    if (any(ploidy < 0)) {
        stop("ploidy values must be positive")
    }

    if (any(cellularity < 0 | cellularity > 1)) {
        stop("cellularity values must be between 0 and 1")
    }
}

# Add progress tracking for long operations
track_fit_progress <- function(total, verbose = TRUE) {
    if (!verbose) {
        return(NULL)
    }
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    function(i) setTxtProgressBar(pb, i)
}

# Helper function to prepare input data
prepare_fit_data <- function(sequenza.extract, chromosome.list) {
    # Extract and combine data for specified chromosomes
    if (is.null(chromosome.list)) {
        mutations <- do.call(rbind, sequenza.extract$mutations)
        segments <- do.call(rbind, sequenza.extract$segments)
    } else {
        mutations <- do.call(rbind, sequenza.extract$mutations[chromosome.list])
        segments <- do.call(rbind, sequenza.extract$segments[chromosome.list])
    }

    # Clean and process data
    mutations <- na.exclude(mutations)
    segs_len <- segments$end.pos - segments$start.pos

    list(mutations = mutations, segs = segments, segs_len = segs_len)
}

# Process data using BAF method
process_baf_method <- function(segs, segs_len, filters, params,
    avg.depth.ratio) {
    # Process in chunks if data is large
    if (nrow(segs) > 1e+05) {
        chunk_size <- 1e+05
        n_chunks <- ceiling(nrow(segs)/chunk_size)

        results <- vector("list", n_chunks)
        progress <- track_fit_progress(n_chunks)

        for (i in seq_len(n_chunks)) {
            idx <- ((i - 1) * chunk_size + 1):min(i * chunk_size,
                nrow(segs))
            results[[i]] <- process_baf_chunk(segs[idx, ], segs_len[idx],
                filters, params, avg.depth.ratio)
            progress(i)
        }

        # Combine results
        return(combine_baf_results(results))
    }

    # Calculate average standard deviations
    sd_stats <- calculate_segment_sds(segs)
    segs <- adjust_zero_sds(segs, sd_stats)

    # Apply filters
    filtered <- filter_segments(segs, segs_len, filters, params$female,
        params$XY)

    # Calculate segment lengths in megabases
    seg_len_mb <- segs_len[filtered$mask]/1e+06

    # Use baf.model.fit from bayes.R
    baf.model.fit(Bf = filtered$data$Bf, depth.ratio = filtered$data$depth.ratio,
        sd.ratio = filtered$data$sd.ratio, weight.ratio = seg_len_mb,
        sd.Bf = filtered$data$sd.BAF, weight.Bf = seg_len_mb,
        avg.depth.ratio = avg.depth.ratio, cellularity = params$cellularity,
        ploidy = params$ploidy, priors.table = params$priors.table,
        mc.cores = params$mc.cores, ratio.priority = params$ratio.priority)
}

# Process data using mutation frequency method
process_mufreq_method <- function(mutations, params, avg.depth.ratio) {
    # Apply filters
    mask <- mutations$F >= params$threshold
    if (params$female) {
        xy_mask <- mutations$chromosome == params$XY["Y"]
    } else {
        xy_mask <- mutations$chromosome %in% params$XY
    }

    filtered <- mutations[mask & !xy_mask, ]
    weights <- round(filtered$good.reads, 0)

    mufreq.model.fit(mufreq = filtered$F, depth.ratio = filtered$adjusted.ratio,
        weight.ratio = 2 * weights, weight.mufreq = weights,
        avg.depth.ratio = avg.depth.ratio, cellularity = params$cellularity,
        ploidy = params$ploidy, priors.table = params$priors.table,
        mc.cores = params$mc.cores)
}

# Calculate average standard deviations for segments
calculate_segment_sds <- function(segs) {
    avg.sd.ratio <- sum(segs$sd.ratio * segs$N.ratio, na.rm = TRUE)/sum(segs$N.ratio,
        na.rm = TRUE)
    avg.sd.Bf <- sum(segs$sd.BAF * segs$N.BAF, na.rm = TRUE)/sum(segs$N.BAF,
        na.rm = TRUE)
    list(sd.ratio = avg.sd.ratio, sd.Bf = avg.sd.Bf)
}

# Adjust segments with zero standard deviations
adjust_zero_sds <- function(segs, sd_stats) {
    segs$sd.BAF[segs$sd.BAF == 0] <- max(segs$sd.BAF, na.rm = TRUE)
    segs$sd.ratio[segs$sd.ratio == 0] <- max(segs$sd.ratio, na.rm = TRUE)
    segs
}

# Improve filter_segments with better error handling
filter_segments <- function(segs, segs_len, filters, female,
    XY) {
    # Add input validation
    if (is.null(segs) || nrow(segs) == 0) {
        stop("No segments to filter")
    }

    segs.filt <- segs$N.ratio > filters$N.ratio & segs$N.BAF >
        filters$N.BAF
    segs.filt <- segs_len >= filters$segment & segs.filt
    if (female) {
        segs.is.xy <- segs$chromosome == XY["Y"]
    } else {
        segs.is.xy <- segs$chromosome %in% XY
    }
    filt.test <- segs.filt & !segs.is.xy
    list(data = segs[filt.test, ], mask = filt.test)
}
