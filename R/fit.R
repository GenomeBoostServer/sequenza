#' @rdname sequenza
#' @export
sequenza.fit <- function(sequenza.extract, female = TRUE, N.ratio.filter = 10,
    N.BAF.filter = 1, segment.filter = 3e+06, mufreq.threshold = 0.1,
    XY = c(X = "X", Y = "Y"), cellularity = seq(0.1, 1, 0.01),
    ploidy = seq(1, 7, 0.1), ratio.priority = FALSE, method = "baf",
    priors.table = data.frame(CN = 2, value = 2), chromosome.list = 1:24,
    mc.cores = getOption("mc.cores", 2L), verbose = TRUE) {

    # Validate sequenza.extract is a list with required
    # components
    if (!is.list(sequenza.extract)) {
        stop("sequenza.extract must be a list object from sequenza.extract()")
    }

    required_components <- c("mutations", "segments", "avg.depth.ratio")
    missing <- setdiff(required_components, names(sequenza.extract))
    if (length(missing) > 0) {
        stop("Missing required components in sequenza.extract: ",
            paste(missing, collapse = ", "))
    }

    # Add progress tracking
    if (verbose)
        message("Starting sequenza fit analysis...")

    # Prepare input data with error handling
    data <- safely_execute({
        prepare_fit_data(sequenza.extract, chromosome.list)
    }, error_message = "Error preparing data")

    if (is.null(data)) {
        stop("Failed to prepare input data")
    }

    avg.depth.ratio <- sequenza.extract$avg.depth.ratio

    # Process method with better error handling and
    # progress tracking
    result <- safely_execute({
        if (method == "baf") {
            message("Processing BAF method...")
            process_baf_method(data$segs, data$segs_len, filters = list(N.ratio = N.ratio.filter,
                N.BAF = N.BAF.filter, segment = segment.filter),
                params = list(female = female, XY = XY, cellularity = cellularity,
                  ploidy = ploidy, ratio.priority = ratio.priority,
                  priors.table = priors.table, mc.cores = mc.cores,
                  verbose = verbose), avg.depth.ratio = avg.depth.ratio)
        } else if (method == "mufreq") {
            message("Processing mutation frequency method...")
            process_mufreq_method(data$mutations, params = list(female = female,
                XY = XY, cellularity = cellularity, ploidy = ploidy,
                threshold = mufreq.threshold, priors.table = priors.table,
                mc.cores = mc.cores), avg.depth.ratio = avg.depth.ratio)
        }
    }, error_message = sprintf("Error in %s method", method))

    if (is.null(result)) {
        stop("Method processing failed")
    }

    # Add metadata to results with validation
    result$params <- list(method = method, female = female, filters = list(N.ratio = N.ratio.filter,
        N.BAF = N.BAF.filter, segment = segment.filter), chromosomes = chromosome.list)

    if (verbose)
        message("Analysis completed successfully")
    return(result)
}

# Add improved data preparation with validation
prepare_fit_data <- function(sequenza.extract, chromosome.list) {
    # Extract and combine data for specified chromosomes
    mutations <- if (is.null(chromosome.list)) {
        do.call(rbind, sequenza.extract$mutations)
    } else {
        do.call(rbind, sequenza.extract$mutations[chromosome.list])
    }

    segments <- if (is.null(chromosome.list)) {
        do.call(rbind, sequenza.extract$segments)
    } else {
        do.call(rbind, sequenza.extract$segments[chromosome.list])
    }

    # Validate segments data
    segments <- validate_data_frame(segments, c("start.pos",
        "end.pos", "chromosome"), "Segments")
    if (is.null(segments)) {
        stop("Invalid segments data")
    }

    # Clean and process data
    mutations <- na.exclude(mutations)
    segs_len <- segments$end.pos - segments$start.pos

    list(mutations = mutations, segs = segments, segs_len = segs_len)
}

# Process BAF method with better validation
process_baf_method <- function(segs, segs_len, filters, params,
    avg.depth.ratio) {
    # Validate required columns for BAF method
    segs <- validate_data_frame(segs, c("sd.ratio", "sd.BAF",
        "N.ratio", "N.BAF", "Bf", "depth.ratio"), "BAF segments")
    if (is.null(segs)) {
        stop("Invalid segment data for BAF method")
    }

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

# Process mutation frequency method with validation
process_mufreq_method <- function(mutations, params, avg.depth.ratio) {
    # Validate mutations data
    mutations <- validate_data_frame(mutations, c("F", "good.reads",
        "adjusted.ratio", "chromosome"), "Mutations")
    if (is.null(mutations)) {
        stop("Invalid mutation data")
    }

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
