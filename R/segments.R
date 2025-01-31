# Core segmentation functions
find_breaks <- function(seqz.baf, slide_win, peak_win, arms, chr_name, verbose) {
  chromosome <- gsub(x = seqz.baf$chromosome, pattern = "chr", replacement = "")
  chromosome <- paste0("chr", chromosome)
  if (verbose) {
    message("Segmenting depth ratios")
  }
  ratio_diffs <- slide_matrix(seqz.baf$adjusted.ratio,
    w = slide_win, position = seqz.baf$position,
    verbose = verbose
  )
  if (verbose) {
    message("Segmenting allele frequencies")
  }
  bf_diffs <- slide_matrix(seqz.baf$Bf,
    w = slide_win, position = seqz.baf$position,
    verbose = verbose
  )

  peaks_both <- get_gaps_peaks(
    x = (ratio_diffs$y + bf_diffs$y) / 2, position = ratio_diffs$x,
    w = peak_win, arms = arms
  )
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

slide_tracks <- function(
  seqz.baf, slide_win, signal_out = c("both", "ratio", "baf"),
  verbose = TRUE
) {
  signal_out <- match.arg(arg = signal_out, choices = signal_out)
  chromosome <- gsub(x = seqz.baf$chromosome, pattern = "chr", replacement = "")
  chromosome <- paste0("chr", chromosome)
  if (signal_out %in% c("both", "ratio")) {
    if (verbose) {
      message("Segmenting depth ratios")
    }
    ratio_diffs <- slide_matrix(seqz.baf$adjusted.ratio,
      w = slide_win, position = seqz.baf$position,
      verbose = verbose
    )
  }

  if (signal_out %in% c("both", "baf")) {
    if (verbose) {
      message("Segmenting allele frequencies")
    }
    bf_diffs <- slide_matrix(seqz.baf$Bf,
      w = slide_win, position = seqz.baf$position,
      verbose = verbose
    )
  }

  if (signal_out == "both") {
    data.frame(y = (ratio_diffs$y + bf_diffs$y) / 2, x = ratio_diffs$x)
  } else if (signal_out == "baf") {
    bf_diffs
  } else {
    ratio_diffs
  }
}

peaks_tracks <- function(diff_track, peak_win, arms, chr_name, verbose) {
  peaks <- get_gaps_peaks(
    x = diff_track$y, position = diff_track$x, w = peak_win,
    arms = arms
  )

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


extract_breaks <- function(
  data, data_het, breaks, slide_win, peak_win, assembly,
  chromosome, verbose = TRUE
) {
  if (is.null(breaks)) {
    golden_path <- paste("http://hgdownload.cse.ucsc.edu", "goldenPath", assembly,
      "database", "cytoBandIdeo.txt.gz",
      sep = "/"
    )
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
      "database", "cytoBand.txt.gz",
      sep = "/"
    )
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
  required_cols <- c(
    "chromosome", "position", "zygosity.normal", "good.reads",
    "Af", "Bf", "depth.normal", "adjusted.ratio"
  )
  missing_cols <- setdiff(required_cols, names(seqz.tab))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  # Calculate weighted values if requested
  if (weighted.mean) {
    # Calculate weighted ratios and B-allele frequencies
    w.r <- sqrt(seqz.tab$depth.normal) # Weights for depth ratio
    rw <- seqz.tab$adjusted.ratio * w.r # Weighted depth ratios
    w.b <- sqrt(seqz.tab$good.reads) # Weights for BAF
    bw <- seqz.tab$Bf * w.b # Weighted BAF

    # Combine original data with weighted values
    seqz.tab <- cbind(seqz.tab[, c(
      "chromosome", "position", "zygosity.normal",
      "good.reads", "Af", "Bf"
    )], rw = rw, w.r = w.r, bw = bw, w.b = w.b)
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
    breaks.vect <- do.call(cbind, split.data.frame(breaks.i[, c(
      "start.pos",
      "end.pos"
    )], f = 1:nb))
    breaks.vect <- unique.breaks(b = as.numeric(breaks.vect), offset = 1)

    # Cut data into segments
    fact.r.i <- cut(seqz.tab[[i]]$position, breaks.vect)
    fact.b.i <- cut(seqz.b.i$position, breaks.vect)

    # Count points in each segment
    seg.i.s.r <- sapply(split(seqz.tab[[i]]$chromosome, f = fact.r.i), length)
    seg.i.s.b <- sapply(split(seqz.b.i$chromosome, f = fact.b.i), length)

    # Calculate segment statistics based on weighting method
    if (weighted.mean) {
      segments.i <- calculate_weighted_segments(
        seqz_data = seqz.tab[[i]],
        baf_data = seqz.b.i, fact_r = fact.r.i, fact_b = fact.b.i, breaks_vect = breaks.vect,
        chrom_name = current_chrom, seg.i.s.r = seg.i.s.r, seg.i.s.b = seg.i.s.b
      )
    } else {
      segments.i <- calculate_unweighted_segments(
        seqz_data = seqz.tab[[i]],
        baf_data = seqz.b.i, fact_r = fact.r.i, fact_b = fact.b.i, breaks_vect = breaks.vect,
        chrom_name = current_chrom, seg.i.s.r = seg.i.s.r, seg.i.s.b = seg.i.s.b
      )
    }

    # Keep every other segment (odd-numbered segments)
    segments[[i]] <- segments.i[seq(1, nrow(segments.i), by = 2), ]
  }

  # Combine segments from all chromosomes
  segments <- do.call(rbind, segments[as.factor(chr.order)])
  row.names(segments) <- seq_len(nrow(segments))

  # Filter segments based on density (points per Mb)
  len.seg <- (segments$end.pos - segments$start.pos) / 1e+06
  segments[(segments$N.ratio / len.seg) >= 2, ]
}

# Helper functions for segmentation
calculate_weighted_segments <- function(
  seqz_data, baf_data, fact_r, fact_b, breaks_vect,
  chrom_name, seg.i.s.r, seg.i.s.b
) {
  # Calculate weighted sums and standard deviations
  seg.i.rw <- sapply(split(seqz_data$rw, f = fact_r), sum, na.rm = TRUE)
  seg.i.w.r <- sapply(split(seqz_data$w.r, f = fact_r), sum, na.rm = TRUE)

  # Calculate standard deviations
  seg.i.r.sd <- sapply(split(seqz_data$rw / seqz_data$w.r, f = fact_r), sd, na.rm = TRUE)
  seg.i.b.sd <- sapply(split(baf_data$bw / baf_data$w.b, f = fact_b), sd, na.rm = TRUE)

  # Split BAF data
  A.split <- split(baf_data$Af, f = fact_b)
  B.split <- split(baf_data$Bf, f = fact_b)
  d.split <- split(baf_data$good.reads, f = fact_b)

  # Calculate B-allele frequencies
  window.quantiles <- mapply(b_allele_freq,
    Af = A.split, Bf = B.split, good.reads = d.split,
    conf = 0.95
  )

  # Create segments dataframe
  data.frame(
    chromosome = chrom_name, start.pos = as.numeric(breaks_vect[-length(breaks_vect)]),
    end.pos = as.numeric(breaks_vect[-1]), Bf = window.quantiles[2, ], N.BAF = seg.i.s.b,
    sd.BAF = seg.i.b.sd, depth.ratio = seg.i.rw / seg.i.w.r, N.ratio = seg.i.s.r,
    sd.ratio = seg.i.r.sd, stringsAsFactors = FALSE
  )
}

calculate_unweighted_segments <- function(
  seqz_data, baf_data, fact_r, fact_b, breaks_vect,
  chrom_name, seg.i.s.r, seg.i.s.b
) {
  # Calculate means and standard deviations
  seg.i.r <- sapply(split(seqz_data$adjusted.ratio, f = fact_r), mean, na.rm = TRUE)
  seg.i.r.sd <- sapply(split(seqz_data$adjusted.ratio, f = fact_r), sd, na.rm = TRUE)
  seg.i.b.sd <- sapply(split(baf_data$Bf, f = fact_b), sd, na.rm = TRUE)

  # Split BAF data
  A.split <- split(baf_data$Af, f = fact_b)
  B.split <- split(baf_data$Bf, f = fact_b)
  d.split <- split(baf_data$good.reads, f = fact_b)

  # Calculate B-allele frequencies
  window.quantiles <- mapply(b_allele_freq,
    Af = A.split, Bf = B.split, good.reads = d.split,
    conf = 0.95
  )

  # Create segments dataframe
  data.frame(
    chromosome = chrom_name, start.pos = as.numeric(breaks_vect[-length(breaks_vect)]),
    end.pos = as.numeric(breaks_vect[-1]), Bf = window.quantiles[2, ], N.BAF = seg.i.s.b,
    sd.BAF = seg.i.b.sd, depth.ratio = seg.i.r, N.ratio = seg.i.s.r, sd.ratio = seg.i.r.sd,
    stringsAsFactors = FALSE
  )
}

compare_bins <- function(start, end, value, bins) {
  # Input validation
  if (length(start) == 0 || length(end) == 0 || length(value) == 0 || nrow(bins) ==
    0) {
    warning("Empty input data in compare_bins")
    return(list(fit = 0, penalty = 1))
  }

  # Create data frame and validate
  segs_vals <- tryCatch(
    {
      df <- data.frame(start = as.numeric(start), end = as.numeric(end), value = as.numeric(value))
      if (any(is.na(df))) {
        warning("NA values found in segment data")
        return(list(fit = 0, penalty = 1))
      }
      df
    },
    error = function(e) {
      warning("Error creating segments dataframe: ", e$message)
      return(NULL)
    }
  )

  if (is.null(segs_vals)) {
    return(list(fit = 0, penalty = 1))
  }

  # Helper function to get overlapping segments
  get_segs <- function(start, end, segs) {
    idx <- which(segs$start < end & segs$end > start)
    if (length(idx) == 0) {
      return(integer(0))
    }
    return(idx)
  }

  # Process each bin with consecutive outlier detection
  bin_results <- tryCatch(
    {
      results <- apply(bins, 1, FUN = function(x, segs) {
        start <- as.numeric(x[1])
        end <- as.numeric(x[2])
        q0 <- as.numeric(x[4])
        q1 <- as.numeric(x[5])

        if (any(is.na(c(start, end, q0, q1)))) {
          return(c(fit = 0, penalty = 1))
        }

        indexes <- get_segs(start, end, segs)
        if (length(indexes) == 0) {
          return(c(fit = 0, penalty = 1))
        }

        segs_values <- segs$value[indexes]
        if (length(segs_values) == 0) {
          return(c(fit = 0, penalty = 1))
        }

        # Calculate fit and consecutive outlier penalties
        in_bounds <- segs_values >= q0 & segs_values <= q1
        fit_score <- mean(in_bounds, na.rm = TRUE)

        # Handle outlier runs with more severe penalties
        out_bounds <- !in_bounds
        if (length(out_bounds) > 1) {
          runs <- rle(out_bounds)
          outlier_runs <- runs$lengths[runs$values]
          if (length(outlier_runs) > 0) {
            # Apply quadratic penalty for longer runs
            max_run <- max(outlier_runs)
            run_penalty <- (max_run / length(out_bounds))^2

            # Additional penalty for multiple outlier runs
            if (length(outlier_runs) > 1) {
              run_penalty <- run_penalty * (1 + 0.1 * length(outlier_runs))
            }
          } else {
            run_penalty <- 0
          }
        } else {
          run_penalty <- if (length(out_bounds) == 1 && out_bounds) {
            1
          } else {
            0
          }
        }

        c(fit = fit_score, penalty = run_penalty)
      }, segs = segs_vals)

      if (is.null(results)) {
        return(NULL)
      }
      results
    },
    error = function(e) {
      warning("Error in bin comparison: ", e$message)
      return(NULL)
    }
  )

  if (is.null(bin_results)) {
    return(list(fit = 0, penalty = 1))
  }

  # Return both scores separately instead of combining them
  list(fit = mean(bin_results[1, ], na.rm = TRUE), penalty = mean(bin_results[2, ], na.rm = TRUE))
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


process_segments_by_clusters <- function(
  sequenza_extract, seqz_file, out_path, file_out_prefix,
  init_n_clust = 10, dp_iter = 1000, pdf_out = FALSE, verbose = FALSE, ...
) {
  segs_i <- do.call(rbind, sequenza_extract$segments)
  gc_stats <- sequenza_extract$gc
  segs_fitting_ratio <- sapply(sequenza_extract$chromosomes, FUN = function(x,
                                                                            extr) {
    compare_bins(
      extr$segments[[x]]$start.pos, extr$segments[[x]]$end.pos, extr$segments[[x]]$depth.ratio,
      extr$ratio[[x]]
    )
  }, extr = sequenza_extract)
  segs_fitting_baf <- sapply(sequenza_extract$chromosomes, FUN = function(x, extr) {
    compare_bins(
      extr$segments[[x]]$start.pos, extr$segments[[x]]$end.pos, extr$segments[[x]]$Bf,
      extr$BAF[[x]]
    )
  }, extr = sequenza_extract)
  if (pdf_out) {
    segs_i_clust <- file.path(out_path, paste(file_out_prefix, paste0(
      "iter",
      i
    ), "reclust.pdf", sep = "_"))
    pdf(segs_i_clust)
  }
  segs_clust <- cluster_segments(
    bf = segs_i$Bf, depth_ratio = segs_i$depth.ratio,
    init_clust = init_n_clust, progressbar = TRUE, iters = dp_iter, plots = pdf_out
  )
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
  temp_extract <- sequenza.extract(seqz_file,
    breaks = seg_res, gc.stats = gc_stats,
    verbose = verbose, chromosome.list = sequenza_extract$chromosomes, ...
  )
  segs_fitting_ratio_i <- sapply(temp_extract$chromosomes, FUN = function(x, extr) {
    compare_bins(
      extr$segments[[x]]$start.pos, extr$segments[[x]]$end.pos, extr$segments[[x]]$depth.ratio,
      extr$ratio[[x]]
    )
  }, extr = temp_extract)
  segs_fitting_baf_i <- sapply(temp_extract$chromosomes, FUN = function(x, extr) {
    compare_bins(
      extr$segments[[x]]$start.pos, extr$segments[[x]]$end.pos, extr$segments[[x]]$Bf,
      extr$BAF[[x]]
    )
  }, extr = temp_extract)

  fit_table <- data.frame(
    chromosome = sequenza_extract$chromosomes, ratio_fit = segs_fitting_ratio[sequenza_extract$chromosomes],
    baf_fit = segs_fitting_baf[sequenza_extract$chromosomes], N = sapply(
      sequenza_extract$segments[sequenza_extract$chromosomes],
      nrow
    )
  )

  refit_table <- data.frame(
    chromosome = temp_extract$chromosomes, ratio_fit = segs_fitting_ratio_i[temp_extract$chromosomes],
    baf_fit = segs_fitting_baf_i[temp_extract$chromosomes], N = sapply(
      temp_extract$segments[temp_extract$chromosomes],
      nrow
    )
  )
  list(extract = temp_extract, fit_table = fit_table, refit_table = refit_table)
}

rank_segments <- function(breaks_list, windows, params) {
  # Early return for single entry
  if (length(breaks_list) <= 1) {
    return(list(
      segs = breaks_list[[1]], selected_win = params$peak_wins[1],
      peak_win = data.frame(
        peak_win = params$peak_wins[1], baf_fit = 0, ratio_fit = 0,
        n_segs = nrow(breaks_list[[1]])
      )
    ))
  }

  # Validate windows input
  if (!all(c("baf", "ratio") %in% names(windows))) {
    stop("Windows must contain both 'baf' and 'ratio' components")
  }

  # Pre-extract window data for efficiency
  baf_win <- windows$baf[[1]]
  ratio_win <- windows$ratio[[1]]

  # Calculate segment comparisons
  compare_bins_list <- if (params$parallel > 1) {
    # ...parallel calculation code...
  } else {
    lapply(breaks_list, function(x) {
      baf_results <- safely_compare_bins(x$start.pos, x$end.pos, x$Bf, baf_win)
      ratio_results <- safely_compare_bins(
        x$start.pos, x$end.pos, x$depth.ratio,
        ratio_win
      )
      c(
        baf_fit = baf_results$fit, baf_penalty = baf_results$penalty, ratio_fit = ratio_results$fit,
        ratio_penalty = ratio_results$penalty
      )
    })
  }

  # Create comparison dataframe
  compare_bins_segs <- data.frame(
    peak_win = params$peak_wins, do.call(rbind, compare_bins_list),
    n_segs = vapply(breaks_list, nrow, numeric(1))
  )

  # Calculate scores using calculate_segment_scores
  scores <- tryCatch(
    {
      calculate_segment_scores(compare_bins_segs, params$segment_weights)
    },
    error = function(e) {
      warning("Error in score calculation: ", e$message)
      # Fallback scoring if calculation fails
      rep(0, nrow(compare_bins_segs))
    }
  )

  # Ensure we have valid scores
  if (all(is.na(scores)) || length(scores) == 0) {
    scores <- rep(0, nrow(compare_bins_segs))
    warning("No valid scores calculated, using default scoring")
  }

  # Select best window size based on maximum score
  best_idx <- which.max(scores)
  if (length(best_idx) == 0) {
    best_idx <- 1
  }
  select_win <- compare_bins_segs$peak_win[best_idx]

  # Create lookup table for breaks_list indices
  win_to_idx <- setNames(seq_along(params$peak_wins), as.character(params$peak_wins))
  breaks_idx <- win_to_idx[as.character(select_win)]

  # Add scores to output
  compare_bins_segs$composite_score <- scores

  if (params$verbose) {
    # ...verbose output code...
  }

  list(
    segs = breaks_list[[breaks_idx]], selected_win = select_win, peak_win = compare_bins_segs,
    weights = params$segment_weights,  # Store weights in output
    chromosome = params$chromosome     # Store chromosome name
  )
}

# Add error handling wrapper
safely_compare_bins <- function(start.pos, end.pos, values, windows) {
  tryCatch(
    {
      compare_bins(start.pos, end.pos, values, windows)
    },
    error = function(e) {
      message("Warning: Bin comparison failed: ", e$message)
      return(list(fit = 0, penalty = 1)) # Return neutral score on failure
    }
  )
}

# Improved rank_segments with better error handling and optimization
rank_segments <- function(breaks_list, windows, params) {
  if (length(breaks_list) <= 1) {
    return(list(
      segs = breaks_list[[1]], selected_win = params$peak_wins[1],
      peak_win = data.frame(
        peak_win = params$peak_wins[1], baf_fit = 0, ratio_fit = 0,
        n_segs = nrow(breaks_list[[1]])
      )
    ))
  }

  # Validate windows input
  if (!all(c("baf", "ratio") %in% names(windows))) {
    stop("Windows must contain both 'baf' and 'ratio' components")
  }

  # Pre-extract window data for efficiency
  baf_win <- windows$baf[[1]]
  ratio_win <- windows$ratio[[1]]

  # Calculate segment comparisons in parallel if possible
  compare_bins_list <- if (params$parallel > 1) {
    pbapply::pblapply(breaks_list, function(x) {
      baf_results <- safely_compare_bins(x$start.pos, x$end.pos, x$Bf, baf_win)
      ratio_results <- safely_compare_bins(
        x$start.pos, x$end.pos, x$depth.ratio,
        ratio_win
      )
      c(
        baf_fit = baf_results$fit, baf_penalty = baf_results$penalty, ratio_fit = ratio_results$fit,
        ratio_penalty = ratio_results$penalty
      )
    }, cl = params$parallel)
  } else {
    lapply(breaks_list, function(x) {
      baf_results <- safely_compare_bins(x$start.pos, x$end.pos, x$Bf, baf_win)
      ratio_results <- safely_compare_bins(
        x$start.pos, x$end.pos, x$depth.ratio,
        ratio_win
      )
      c(
        baf_fit = baf_results$fit, baf_penalty = baf_results$penalty, ratio_fit = ratio_results$fit,
        ratio_penalty = ratio_results$penalty
      )
    })
  }

  # Create comparison dataframe
  compare_bins_segs <- data.frame(
    peak_win = params$peak_wins, do.call(rbind, compare_bins_list),
    n_segs = vapply(breaks_list, nrow, numeric(1))
  )

  # Calculate metrics with focus on first significant drop
  calculate_segment_scores <- function(compare_bins_segs, weights) {
    # Normalize metrics to 0-1 scale
    normalize <- function(x) (x - min(x)) / (max(x) - min(x))

    # Calculate normalized scores
    baf_score <- normalize(compare_bins_segs$baf_fit)
    ratio_score <- normalize(compare_bins_segs$ratio_fit)

    # Calculate combined fit score
    combined_fit <- (baf_score + ratio_score) / 2

    # Find the elbow point using curvature
    find_elbow <- function(y) {
      # Handle edge cases
      if (length(y) < 3) {
        return(1)
      }
      if (all(is.na(y))) {
        return(1)
      }

      # Remove any NA values while preserving position information
      valid_idx <- which(!is.na(y))
      if (length(valid_idx) < 3) {
        return(1)
      }

      y_clean <- y[valid_idx]
      x <- seq_along(y_clean)

      # Normalize x and y to 0-1 scale
      x_norm <- (x - min(x)) / (diff(range(x)))
      y_norm <- (y_clean - min(y_clean)) / (diff(range(y_clean)))

      # Calculate first derivatives
      dy <- diff(y_norm)
      dx <- diff(x_norm)

      # Ensure we have enough points for second derivative
      if (length(dy) < 2) {
        return(1)
      }

      # Calculate second derivatives
      dy2 <- diff(dy)
      dx2 <- diff(dx)

      # Calculate curvature only where derivatives are defined
      curvature <- rep(0, length(y_clean))
      curve_idx <- 2:(length(y_clean) - 1)

      # Calculate curvature using finite differences
      for (i in curve_idx) {
        denom <- (1 + (dy[i] / dx[i])^2)^(3 / 2)
        if (!is.na(denom) && denom != 0) {
          curvature[i] <- abs(dy2[i - 1] / dx2[i - 1]) / denom
        }
      }

      # Only consider points where fit is improving (negative slope)
      valid_points <- dy < 0
      if (all(!valid_points, na.rm = TRUE)) {
        return(1)
      }

      curvature[!valid_points] <- 0

      # Map back to original indices
      elbow_idx_local <- which.max(curvature)
      return(valid_idx[elbow_idx_local])
    }

    # Find elbow point
    elbow_idx <- find_elbow(combined_fit)

    # Calculate distance score from elbow point
    distance_from_elbow <- abs(seq_along(combined_fit) - elbow_idx)
    elbow_score <- 1 - normalize(distance_from_elbow)

    # Segment count bonus (small preference for more segments up to elbow
    # point)
    n_segs <- compare_bins_segs$n_segs
    segment_bonus <- rep(0, length(n_segs))
    segment_bonus[1:elbow_idx] <- normalize(n_segs[1:elbow_idx]) * 0.1

    # Window size penalty (prefer smaller windows when fits are similar)
    window_sizes <- compare_bins_segs$peak_win
    window_penalty <- normalize(window_sizes) * 0.05

    # Combine scores with emphasis on elbow point
    final_scores <- weights$fit * combined_fit + weights$elbow * elbow_score +
      weights$segments * segment_bonus - weights$window * window_penalty

    if (params$verbose) {
      message("\nElbow point analysis:")
      message("Elbow detected at window size: ", compare_bins_segs$peak_win[elbow_idx])
      message("Fit score at elbow: ", round(combined_fit[elbow_idx], 4))
      message("Number of segments at elbow: ", n_segs[elbow_idx])
    }

    return(final_scores)
  }

  # Calculate comprehensive scores
  scores <- calculate_segment_scores(compare_bins_segs, params$segment_weights)

  # Select best window size based on maximum score
  best_idx <- which.max(scores)
  select_win <- compare_bins_segs$peak_win[best_idx]

  # Add scores to output for debugging
  compare_bins_segs$composite_score <- scores

  if (params$verbose) {
    message("Segment selection results:")
    message("Selected window size: ", select_win)
    message("Number of segments: ", compare_bins_segs$n_segs[best_idx])
    message("Composite score: ", round(scores[best_idx], 4))

    # Add more detailed diagnostics
    message("\nTop 3 solutions:")
    top3 <- head(order(scores, decreasing = TRUE), 3)
    for (i in top3) {
      message(sprintf(
        "Window: %d, Segments: %d, Score: %.4f", compare_bins_segs$peak_win[i],
        compare_bins_segs$n_segs[i], scores[i]
      ))
    }
  }

  list(
    segs = breaks_list[[as.character(select_win)]], selected_win = select_win,
    peak_win = compare_bins_segs,
    weights = params$segment_weights,  # Store weights in output
    chromosome = params$chromosome     # Store chromosome name
  )
}

# Update process_segments to use new rank_segments
process_segments <- function(seqz.data, breaks, chr, windows, params) {
  # Ensure weighted.mean has a default value if not in params
  weighted.mean <- if (!is.null(params$weighted.mean)) {
    params$weighted.mean
  } else {
    TRUE # Default value
  }

  # Handle segmentation
  if (is.null(breaks)) {
    diff_track <- slide_tracks(seqz.data, params$slide_win,
      signal_out = "both",
      verbose = params$verbose
    )

    breaks_chr_list <- lapply(params$peak_wins, function(x) {
      breaks_chr <- extract_breaks_tracks(
        track = diff_track, breaks = breaks,
        peak_win = x, assembly = params$assembly, chromosome = chr
      )

      if (inherits(breaks_chr, "try-error") || is.null(breaks_chr) || nrow(breaks_chr) ==
        0 || length(breaks_chr) == 0) {
        breaks_chr <- data.frame(chrom = chr, start.pos = min(seqz.data$position,
          na.rm = TRUE
        ), end.pos = max(seqz.data$position, na.rm = TRUE))
      }

      tryCatch(
        {
          segment.breaks(
            seqz.tab = seqz.data, breaks = breaks_chr, min.reads.baf = params$min.reads.baf,
            weighted.mean = weighted.mean # Use local variable
          )
        },
        error = function(e) {
          message("Warning: Segment calculation failed: ", e$message)
          data.frame(
            chrom = chr, start.pos = min(seqz.data$position, na.rm = TRUE),
            end.pos = max(seqz.data$position, na.rm = TRUE), Bf = 0, depth.ratio = mean(seqz.data$depth.ratio,
              na.rm = TRUE
            )
          )
        }
      )
    })

    names(breaks_chr_list) <- as.character(params$peak_wins)
    # Store chromosome name in params for rank_segments
    params$chromosome <- chr
    # Return segment rank and information
    segment_results <- rank_segments(breaks_chr_list, windows, params)
    return(list(
      seg = segment_results$segs, breaks_list = breaks_chr_list, selected_win = segment_results$selected_win,
      peak_win = segment_results$peak_win
    ))
  } else {
    segs <- segment.breaks(
      seqz.tab = seqz.data, breaks = breaks, min.reads.baf = params$min.reads.baf,
      weighted.mean = params$weighted.mean
    )
    select_win <- 0
    compare_bins_segs <- data.frame(peak_win = 0, baf_fit = compare_bins(
      segs$start.pos,
      segs$end.pos, segs$Bf, seqz.b.win[[chr]]
    ), ratio_fit = compare_bins(
      segs$start.pos,
      segs$end.pos, segs$depth.ratio, seqz.r.win[[chr]]
    ), n_segs = nrow(segs))

    list(segs = segs, selected_win = select_win, peak_win = compare_bins_segs)
  }
}
