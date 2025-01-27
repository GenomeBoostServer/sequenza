find_breaks <- function(seqz.baf, slide_win, peak_win, arms, chr_name, verbose) {
    chromosome <- gsub(x = seqz.baf$chromosome, pattern = "chr", replacement = "")
    chromosome <- paste0("chr", chromosome)
    if (verbose) {
        message("Segmenting depth ratios")
    }
    ratio_diffs <- slide_matrix(seqz.baf$adjusted.ratio, w = slide_win, position = seqz.baf$position,
        verbose = verbose)
    if (verbose) {
        message("Segmenting allele frequencies")
    }
    bf_diffs <- slide_matrix(seqz.baf$Bf, w = slide_win, position = seqz.baf$position,
        verbose = verbose)

    peaks_both <- get_gaps_peaks(x = (ratio_diffs$y + bf_diffs$y)/2, position = ratio_diffs$x,
        w = peak_win, arms = arms)
    breaks <- lapply(peaks_both, FUN = function(peaks) {
        coords <- peaks
        if (is.null(coords)) {
            NULL
        } else {
            pos_start <- coords[-length(coords)]
            pos_end <- coords[-1]
            data.frame(chrom = chr_name, start.pos = pos_start, end.pos = pos_end)
        }
    })
    breaks <- do.call(rbind, breaks)
    not.uniq <- which(breaks$end.pos == c(breaks$start.pos[-1], 0))
    breaks$end.pos[not.uniq] <- breaks$end.pos[not.uniq] - 1
    breaks
}

slide_tracks <- function(seqz.baf, slide_win, signal_out = c("both", "ratio", "baf"),
    verbose = TRUE) {
    signal_out <- match.arg(arg = signal_out, choices = signal_out)
    chromosome <- gsub(x = seqz.baf$chromosome, pattern = "chr", replacement = "")
    chromosome <- paste0("chr", chromosome)
    if (signal_out %in% c("both", "ratio")) {
        if (verbose) {
            message("Segmenting depth ratios")
        }
        ratio_diffs <- slide_matrix(seqz.baf$adjusted.ratio, w = slide_win, position = seqz.baf$position,
            verbose = verbose)
    }

    if (signal_out %in% c("both", "baf")) {
        if (verbose) {
            message("Segmenting allele frequencies")
        }
        bf_diffs <- slide_matrix(seqz.baf$Bf, w = slide_win, position = seqz.baf$position,
            verbose = verbose)
    }

    if (signal_out == "both") {
        data.frame(y = (ratio_diffs$y + bf_diffs$y)/2, x = ratio_diffs$x)
    } else if (signal_out == "baf") {
        bf_diffs
    } else {
        ratio_diffs
    }
}

peaks_tracks <- function(diff_track, peak_win, arms, chr_name, verbose) {
    peaks <- get_gaps_peaks(x = diff_track$y, position = diff_track$x, w = peak_win,
        arms = arms)

    breaks <- lapply(peaks, FUN = function(peaks) {
        coords <- peaks
        if (is.null(coords)) {
            NULL
        } else {
            pos_start <- coords[-length(coords)]
            pos_end <- coords[-1]
            data.frame(chrom = chr_name, start.pos = pos_start, end.pos = pos_end)
        }
    })
    breaks <- do.call(rbind, breaks)
    not_uniq <- which(breaks$end.pos == c(breaks$start.pos[-1], 0))
    breaks$end.pos[not_uniq] <- breaks$end.pos[not_uniq] - 1
    breaks
}


extract_breaks <- function(data, data_het, breaks, slide_win, peak_win, assembly,
    chromosome, verbose = TRUE) {
    if (is.null(breaks)) {
        golden_path <- paste("http://hgdownload.cse.ucsc.edu", "goldenPath", assembly,
            "database", "cytoBandIdeo.txt.gz", sep = "/")
        arms <- get_assembly(url = golden_path, prefix = "chr")
        chr_arm <- gsub(x = chromosome, pattern = "chr", replacement = "")
        chr_arm <- paste0("chr", chr_arm)

        arms_i <- arms[arms$chromosome == chr_arm, ]
        data_het <- data_het[data_het$chromosome == chromosome, ]
        find_breaks(data_het, slide_win, peak_win, arms_i, chromosome, verbose)
    } else {
        breaks
    }
}

extract_breaks_tracks <- function(track, breaks, slide_win, peak_win, assembly, chromosome) {
    if (is.null(breaks)) {
        golden_path <- paste("http://hgdownload.cse.ucsc.edu", "goldenPath", assembly,
            "database", "cytoBand.txt.gz", sep = "/")
        arms <- get_assembly(url = golden_path, prefix = "chr")
        chr_arm <- gsub(x = chromosome, pattern = "chr", replacement = "")
        chr_arm <- paste0("chr", chr_arm)

        arms_i <- arms[arms$chromosome == chr_arm, ]
        peaks_tracks(track, peak_win, arms_i, chromosome)
    } else {
        breaks
    }
}

segment.breaks <- function(seqz.tab, breaks, min.reads.baf = 1, weighted.mean = TRUE) {
    # Input validation
    if (missing(seqz.tab) || missing(breaks)) {
        stop("Both seqz.tab and breaks arguments are required")
    }

    # Check required columns
    required_cols <- c("chromosome", "position", "zygosity.normal", "good.reads",
        "Af", "Bf", "depth.normal", "adjusted.ratio")
    missing_cols <- setdiff(required_cols, names(seqz.tab))
    if (length(missing_cols) > 0) {
        stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
    }

    # Calculate weighted values if requested
    if (weighted.mean) {
        # Calculate weighted ratios and B-allele frequencies
        w.r <- sqrt(seqz.tab$depth.normal)  # Weights for depth ratio
        rw <- seqz.tab$adjusted.ratio * w.r  # Weighted depth ratios
        w.b <- sqrt(seqz.tab$good.reads)  # Weights for BAF
        bw <- seqz.tab$Bf * w.b  # Weighted BAF

        # Combine original data with weighted values
        seqz.tab <- cbind(seqz.tab[, c("chromosome", "position", "zygosity.normal",
            "good.reads", "Af", "Bf")], rw = rw, w.r = w.r, bw = bw, w.b = w.b)
    }

    # Get chromosome order for final result
    chr.order <- unique(seqz.tab$chromosome)

    # Split data by chromosome
    seqz.tab <- split(seqz.tab, f = seqz.tab$chromosome)
    segments <- list()

    # Helper function to ensure breaks are unique
    unique.breaks <- function(b, offset = 1) {
        while (any(diff(b) == 0)) {
            b[which(diff(b) == 0) + 1] <- b[diff(b) == 0] + offset
        }
        b
    }

    # Process each chromosome
    for (i in seq_along(seqz.tab)) {
        current_chrom <- names(seqz.tab)[i]

        # Filter heterozygous positions with sufficient reads
        seqz.b.i <- seqz.tab[[i]][seqz.tab[[i]]$zygosity.normal == "het", ]
        seqz.b.i <- seqz.b.i[seqz.b.i$good.reads >= min.reads.baf, ]

        # Get breaks for current chromosome
        breaks.i <- breaks[breaks$chrom == current_chrom, ]
        nb <- nrow(breaks.i)

        # Create break points vector
        breaks.vect <- do.call(cbind, split.data.frame(breaks.i[, c("start.pos",
            "end.pos")], f = 1:nb))
        breaks.vect <- unique.breaks(b = as.numeric(breaks.vect), offset = 1)

        # Cut data into segments
        fact.r.i <- cut(seqz.tab[[i]]$position, breaks.vect)
        fact.b.i <- cut(seqz.b.i$position, breaks.vect)

        # Count points in each segment
        seg.i.s.r <- sapply(split(seqz.tab[[i]]$chromosome, f = fact.r.i), length)
        seg.i.s.b <- sapply(split(seqz.b.i$chromosome, f = fact.b.i), length)

        # Calculate segment statistics based on weighting method
        if (weighted.mean) {
            segments.i <- calculate_weighted_segments(seqz_data = seqz.tab[[i]],
                baf_data = seqz.b.i, fact_r = fact.r.i, fact_b = fact.b.i, breaks_vect = breaks.vect,
                chrom_name = current_chrom, seg.i.s.r = seg.i.s.r, seg.i.s.b = seg.i.s.b)
        } else {
            segments.i <- calculate_unweighted_segments(seqz_data = seqz.tab[[i]],
                baf_data = seqz.b.i, fact_r = fact.r.i, fact_b = fact.b.i, breaks_vect = breaks.vect,
                chrom_name = current_chrom, seg.i.s.r = seg.i.s.r, seg.i.s.b = seg.i.s.b)
        }

        # Keep every other segment (odd-numbered segments)
        segments[[i]] <- segments.i[seq(1, nrow(segments.i), by = 2), ]
    }

    # Combine segments from all chromosomes
    segments <- do.call(rbind, segments[as.factor(chr.order)])
    row.names(segments) <- seq_len(nrow(segments))

    # Filter segments based on density (points per Mb)
    len.seg <- (segments$end.pos - segments$start.pos)/1e+06
    segments[(segments$N.ratio/len.seg) >= 2, ]
}

# Helper function for weighted segment calculations
calculate_weighted_segments <- function(seqz_data, baf_data, fact_r, fact_b, breaks_vect,
    chrom_name, seg.i.s.r, seg.i.s.b) {
    # Calculate weighted sums and standard deviations
    seg.i.rw <- sapply(split(seqz_data$rw, f = fact_r), sum, na.rm = TRUE)
    seg.i.w.r <- sapply(split(seqz_data$w.r, f = fact_r), sum, na.rm = TRUE)

    # Calculate standard deviations
    seg.i.r.sd <- sapply(split(seqz_data$rw/seqz_data$w.r, f = fact_r), sd, na.rm = TRUE)
    seg.i.b.sd <- sapply(split(baf_data$bw/baf_data$w.b, f = fact_b), sd, na.rm = TRUE)

    # Split BAF data
    A.split <- split(baf_data$Af, f = fact_b)
    B.split <- split(baf_data$Bf, f = fact_b)
    d.split <- split(baf_data$good.reads, f = fact_b)

    # Calculate B-allele frequencies
    window.quantiles <- mapply(b_allele_freq, Af = A.split, Bf = B.split, good.reads = d.split,
        conf = 0.95)

    # Create segments dataframe
    data.frame(chromosome = chrom_name, start.pos = as.numeric(breaks_vect[-length(breaks_vect)]),
        end.pos = as.numeric(breaks_vect[-1]), Bf = window.quantiles[2, ], N.BAF = seg.i.s.b,
        sd.BAF = seg.i.b.sd, depth.ratio = seg.i.rw/seg.i.w.r, N.ratio = seg.i.s.r,
        sd.ratio = seg.i.r.sd, stringsAsFactors = FALSE)
}

# Helper function for unweighted segment calculations
calculate_unweighted_segments <- function(seqz_data, baf_data, fact_r, fact_b, breaks_vect,
    chrom_name, seg.i.s.r, seg.i.s.b) {
    # Calculate means and standard deviations
    seg.i.r <- sapply(split(seqz_data$adjusted.ratio, f = fact_r), mean, na.rm = TRUE)
    seg.i.r.sd <- sapply(split(seqz_data$adjusted.ratio, f = fact_r), sd, na.rm = TRUE)
    seg.i.b.sd <- sapply(split(baf_data$Bf, f = fact_b), sd, na.rm = TRUE)

    # Split BAF data
    A.split <- split(baf_data$Af, f = fact_b)
    B.split <- split(baf_data$Bf, f = fact_b)
    d.split <- split(baf_data$good.reads, f = fact_b)

    # Calculate B-allele frequencies
    window.quantiles <- mapply(b_allele_freq, Af = A.split, Bf = B.split, good.reads = d.split,
        conf = 0.95)

    # Create segments dataframe
    data.frame(chromosome = chrom_name, start.pos = as.numeric(breaks_vect[-length(breaks_vect)]),
        end.pos = as.numeric(breaks_vect[-1]), Bf = window.quantiles[2, ], N.BAF = seg.i.s.b,
        sd.BAF = seg.i.b.sd, depth.ratio = seg.i.r, N.ratio = seg.i.s.r, sd.ratio = seg.i.r.sd,
        stringsAsFactors = FALSE)
}

compare_bins <- function(start, end, value, bins) {
    segs_vals <- data.frame(start, end, value)
    get_segs <- function(start, end, segs) {
        which(segs$start < end & segs$end > start)
    }
    is_similar <- apply(bins, 1, FUN = function(x, segs) {
        start <- x[1]
        end <- x[2]
        q0 <- x[4]
        q1 <- x[5]
        indexes <- get_segs(start, end, segs)
        segs_values <- segs[indexes, "value"]
        all(segs_values >= q0 & segs_values <= q1)
    }, segs = segs_vals)
    sum(is_similar, na.rm = TRUE)/length(na.exclude(is_similar))
}

merge_segments_clusters <- function(segs, clusters) {
    breaks <- list()
    last_clust <- NULL
    for (i in 1:nrow(segs)) {
        if (is.null(last_clust)) {
            breaks[[1]] <- c(segs$start.pos[i], segs$end.pos[i])
            last_clust <- clusters[i]
        } else {
            if (last_clust == clusters[i]) {
                breaks[[length(breaks)]][2] <- segs$end.pos[i]
                last_clust <- clusters[i]
            } else {
                breaks[[length(breaks) + 1]] <- c(segs$start.pos[i], segs$end.pos[i])
                last_clust <- clusters[i]
            }
        }
    }
    do.call(rbind, lapply(breaks, FUN = function(x) c(min(x), max(x))))
}

cluster_segments <- function(bf, depth_ratio, init_clust = 10, ...) {
    x <- cbind(bf, depth_ratio)
    sequenza:::gibbs(t(x), z_init = sample(1:init_clust, nrow(x), replace = T), ...)
}


process_segments_by_clusters <- function(sequenza_extract, seqz_file, out_path, file_out_prefix,
    init_n_clust = 10, dp_iter = 1000, pdf_out = FALSE, verbose = FALSE, ...) {
    segs_i <- do.call(rbind, sequenza_extract$segments)
    gc_stats <- sequenza_extract$gc
    segs_fitting_ratio <- sapply(sequenza_extract$chromosomes, FUN = function(x,
        extr) {
        compare_bins(extr$segments[[x]]$start.pos, extr$segments[[x]]$end.pos, extr$segments[[x]]$depth.ratio,
            extr$ratio[[x]])
    }, extr = sequenza_extract)
    segs_fitting_baf <- sapply(sequenza_extract$chromosomes, FUN = function(x, extr) {
        compare_bins(extr$segments[[x]]$start.pos, extr$segments[[x]]$end.pos, extr$segments[[x]]$Bf,
            extr$BAF[[x]])
    }, extr = sequenza_extract)
    if (pdf_out) {
        segs_i_clust <- file.path(out_path, paste(file_out_prefix, paste0("iter",
            i), "reclust.pdf", sep = "_"))
        pdf(segs_i_clust)
    }
    segs_clust <- cluster_segments(bf = segs_i$Bf, depth_ratio = segs_i$depth.ratio,
        init_clust = init_n_clust, progressbar = TRUE, iters = dp_iter, plots = pdf_out)
    if (pdf_out) {
        dev.off()
    }
    segs2 <- cbind(segs_i, cluster = segs_clust$cluster())
    if (verbose) {
        message("tot clusters: ", length(unique(segs_clust$cluster())))
    }
    segs_split <- split(segs2, f = segs2$chromosome)

    seg_res <- lapply(segs_split, FUN = function(x) {
        breaks <- merge_segments_clusters(x, x$cluster)
        chrom <- unique(x$chromosome)
        res <- data.frame(chrom, breaks)
        colnames(res) <- c("chrom", "start.pos", "end.pos")
        res
    })
    seg_res <- do.call(rbind, seg_res)
    if (verbose) {
        message("prev. N of segs ", nrow(segs_i))
        message("new N of segs ", nrow(seg_res))
    }
    temp_extract <- sequenza.extract(seqz_file, breaks = seg_res, gc.stats = gc_stats,
        verbose = verbose, chromosome.list = sequenza_extract$chromosomes, ...)
    segs_fitting_ratio_i <- sapply(temp_extract$chromosomes, FUN = function(x, extr) {
        compare_bins(extr$segments[[x]]$start.pos, extr$segments[[x]]$end.pos, extr$segments[[x]]$depth.ratio,
            extr$ratio[[x]])
    }, extr = temp_extract)
    segs_fitting_baf_i <- sapply(temp_extract$chromosomes, FUN = function(x, extr) {
        compare_bins(extr$segments[[x]]$start.pos, extr$segments[[x]]$end.pos, extr$segments[[x]]$Bf,
            extr$BAF[[x]])
    }, extr = temp_extract)

    fit_table <- data.frame(chromosome = sequenza_extract$chromosomes, ratio_fit = segs_fitting_ratio[sequenza_extract$chromosomes],
        baf_fit = segs_fitting_baf[sequenza_extract$chromosomes], N = sapply(sequenza_extract$segments[sequenza_extract$chromosomes],
            nrow))

    refit_table <- data.frame(chromosome = temp_extract$chromosomes, ratio_fit = segs_fitting_ratio_i[temp_extract$chromosomes],
        baf_fit = segs_fitting_baf_i[temp_extract$chromosomes], N = sapply(temp_extract$segments[temp_extract$chromosomes],
            nrow))
    list(extract = temp_extract, fit_table = fit_table, refit_table = refit_table)
}
