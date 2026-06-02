#' @rdname sequenza
#' @export
sequenza.extract <- function(file, window = 1e+06, overlap = 1,
    slide_win = 100, peak_wins = 2^(1:10) * 10, support_threshold = 0.2,
    normalization.method = "mean", ignore.normal = FALSE, verbose = TRUE,
    chromosome.list = NULL, breaks = NULL, min.mut.freq = 0.1,
    min.reads = 40, min.reads.normal = 10, min.reads.baf = 1,
    max.mut.types = 1, min.type.freq = 0.9, min.fw.freq = 0,
    assembly = "hg38", female = TRUE, XY = c(X = "X", Y = "Y"),
    gc.stats = NULL, do_raster = FALSE, smooth_gc = FALSE, min_times_gc = 5,
    gc_grid = 250, parallel = 1, weighted.mean = TRUE, cytoband_file = NULL, ...) {
    # Track start time and memory
    start_time <- Sys.time()

    # Initialize parameters with all arguments
    params <- extract_initialize_parameters(file = file, window = window,
        overlap = overlap, slide_win = slide_win, peak_wins = peak_wins,
        support_threshold = support_threshold, normalization.method = normalization.method,
        ignore.normal = ignore.normal, verbose = verbose, chromosome.list = chromosome.list,
        breaks = breaks, assembly = assembly, female = female,
        XY = XY, gc.stats = gc.stats, do_raster = do_raster,
        smooth_gc = smooth_gc, min_times_gc = min_times_gc, gc_grid = gc_grid,
        parallel = parallel, weighted.mean = weighted.mean, cytoband_file = cytoband_file, ...)

    # Validate input parameters
    extract_validate_params(params)

    tryCatch({
        # Process GC content
        gc_data <- extract_process_gc_content(params$gc.stats, params$normalization.method)

        gc_splines <- list(
            normal = smooth.spline(data.frame(
                gc = as.numeric(names(gc_data$normal_vect)),
                depth = gc_data$normal_vect)),
            tumor = smooth.spline(data.frame(
                gc = as.numeric(names(gc_data$tumor_vect)),
                depth = gc_data$tumor_vect)))

        # Initialize containers
        containers <- initialize_extract_containers(params$chromosome.list)

        # Build per-chromosome file map for selective-read mode
        # Look for chr-specific seqz.gz files alongside the merged file
        chr_file_map <- NULL
        if (!dir.exists(file)) {
            file_dir <- dirname(file)
            # Derive sample id from the merged file name, then match
            # "<sample>_<chrom>.small.seqz.gz" for both hg38 (chr1, chrX)
            # and hg19/b37 (1, X, MT) naming. Anchoring on the sample id
            # also prevents matching the merged file itself.
            sample_id <- sub("\\.(small\\.)?seqz\\.gz$", "", basename(file))
            sample_re <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", sample_id)
            seqz_files <- sort(list.files(file_dir,
                pattern = paste0("^", sample_re, "_[^.]+\\.small\\.seqz\\.gz$"),
                full.names = TRUE))
            if (length(seqz_files) > 0) {
                chr_file_map <- setNames(seqz_files,
                    sub(paste0("^", sample_re, "_(.+)\\.small\\.seqz\\.gz$"), "\\1",
                        basename(seqz_files)))
                if (params$verbose)
                    message("Selective-read mode: found ", length(seqz_files),
                        " per-chromosome seqz files")
            }
        }

        # Resolve file path for a given chromosome
        get_chr_file <- function(chr) {
            if (!is.null(chr_file_map) && chr %in% names(chr_file_map))
                chr_file_map[[chr]]
            else
                file
        }

        # Process each chromosome with improved parallel
        # handling
        if (params$parallel > 1) {
            cl <- NULL
            tryCatch({
                cl <- manage_parallel_cluster(params$parallel)
                if (is.null(cl))
                  stop("Failed to create cluster")
                on.exit(if (!is.null(cl)) try(parallel::stopCluster(cl),
                  silent = TRUE))

                results <- pbapply::pblapply(seq_along(params$chromosome.list),
                  function(idx) {
                    chr <- params$chromosome.list[idx]
                    chr_f <- get_chr_file(chr)
                    extract_process_chromosome(chr, chr_f, params$gc.stats,
                      gc_splines, NULL, params)
                  }, cl = cl)

                # Merge results back maintaining original
                # data structure
                for (idx in seq_along(results)) {
                  chr <- params$chromosome.list[idx]
                  containers$windows.baf[[idx]] <- results[[idx]]$windows.baf[[idx]]
                  containers$windows.ratio[[idx]] <- results[[idx]]$windows.ratio[[idx]]
                  containers$windows.raw_ratio[[idx]] <- results[[idx]]$windows.raw_ratio[[idx]]
                  containers$windows.normal[[idx]] <- results[[idx]]$windows.normal[[idx]]
                  containers$windows.tumor[[idx]] <- results[[idx]]$windows.tumor[[idx]]
                  containers$windows.n_normal[[idx]] <- results[[idx]]$windows.n_normal[[idx]]
                  containers$windows.n_tumor[[idx]] <- results[[idx]]$windows.n_tumor[[idx]]
                  containers$segments.list[[idx]] <- results[[idx]]$segments.list[[idx]]
                  containers$mutation.list[[idx]] <- results[[idx]]$mutation.list[[idx]]
                  containers$norm.gc.list[[idx]] <- results[[idx]]$norm.gc.list[[idx]]
                  containers$rank_peaks.list[[idx]] <- results[[idx]]$rank_peaks.list[[idx]]
                  containers$all_segments[[idx]] <- results[[idx]]$all_segments[[idx]]
                }
            }, error = function(e) {
                message("Error in parallel processing: ", e$message)
                stop(e)
            })
        } else {
            for (chr in params$chromosome.list) {
                chr_f <- get_chr_file(chr)
                containers <- extract_process_chromosome(chr, chr_f,
                    params$gc.stats, gc_splines, containers, params)
            }
        }

        # Finalize and return results
        final_results <- finalize_extract_results(containers, params, params$gc.stats, gc_data)
        gc()  # Force garbage collection

        # Print performance summary if verbose
        if (params$verbose) {
            end_time <- Sys.time()
            end_mem <- gc(reset = FALSE)
            end_mem_used <- sum(end_mem[, 2])

            # Track memory usage of all processes
            if (params$parallel > 1) {
                # Get worker processes from parallel cluster
                n_workers <- if (!is.null(cl)) length(cl) else 0

                total_mem <- end_mem_used  # Start with main process memory

                if (n_workers > 0) {
                  worker_mems <- tryCatch({
                    if (.Platform$OS.type == "unix") {
                      # Use ps for Unix-like systems
                      pids <- unlist(parallel::clusterCall(cl, Sys.getpid))
                      mem_cmd <- sprintf("ps -p %s -o rss=", paste(pids, collapse = ","))
                      as.numeric(system(mem_cmd, intern = TRUE)) / 1024
                    } else {
                      # For Windows, just use main process memory
                      rep(end_mem_used / n_workers, n_workers)
                    }
                  }, error = function(e) {
                    message("Warning: Could not get worker memory usage")
                    rep(0, n_workers)
                  })

                  if (any(is.na(worker_mems)))
                      warning(sum(is.na(worker_mems)), " worker memory values are NA")
                  total_mem <- end_mem_used + sum(worker_mems, na.rm = TRUE)
                }

                message(sprintf("Memory usage (main, R gc peak): %.2f GB",
                    max(0, end_mem_used / 1024)))
                message(sprintf("Memory usage (all, main gc + worker RSS): %.2f GB",
                    max(0, total_mem / 1024)))
                message(sprintf("Number of worker processes: %d", n_workers))
            } else {
                total_mem <- end_mem_used
                message(sprintf("Peak memory usage: %.2f GB", max(0, total_mem / 1024)))
            }

            # Calculate total mutations and covered bases
            total_mutations <- sum(sapply(final_results$mutations, nrow))
            total_bases <- sum(sapply(final_results$segments, function(segs) {
                sum(segs$end.pos - segs$start.pos + 1)
            }))
            total_mb <- total_bases / 1e+06

            message("\nPerformance Summary:")
            message(sprintf("Total time: %.2f minutes",
                as.numeric(difftime(end_time, start_time, units = "mins"))))
            if (params$parallel > 1) {
                message(sprintf("Memory usage (main, R gc peak): %.2f GB",
                    max(0, end_mem_used / 1024)))
                message(sprintf("Memory usage (all, main gc + worker RSS): %.2f GB",
                    max(0, total_mem / 1024)))
                message(sprintf("Number of worker processes: %d", n_workers))
            } else {
                message(sprintf("Peak memory usage: %.2f GB", max(0, total_mem / 1024)))
            }
            message(sprintf("Number of chromosomes processed: %d", length(params$chromosome.list)))
            message(sprintf("Total segments identified: %d",
                sum(sapply(final_results$segments, nrow))))
            message(sprintf("Total mutations detected: %d", total_mutations))
            message(sprintf("Total megabases analyzed: %.1f", total_mb))
            message(sprintf("Mutation rate: %.2f mutations/Mb", total_mutations / total_mb))
        }

        return(final_results)
    }, error = function(e) {
        message("Error in sequenza.extract: ", e$message)
        stop(e)
    }, finally = {
        # Cleanup
        gc()
    })
}

# Helper function to validate input parameters
extract_validate_params <- function(params) {
    required <- c("window", "overlap", "normalization.method")
    missing <- required[!required %in% names(params)]
    if (length(missing) > 0) {
        stop("Missing required parameters: ", paste(missing,
            collapse = ", "))
    }

    if (!params$normalization.method %in% c("mean", "median")) {
        stop("normalization.method must be either 'mean' or 'median'")
    }
}

extract_process_gc_content <- function(gc.stats, normalization.method) {
    if (!is.list(gc.stats) || !all(c("normal", "tumor") %in%
        names(gc.stats))) {
        stop("Invalid gc.stats format")
    }

    tryCatch({
        if (normalization.method == "mean") {
            list(normal_vect = mean_gc(gc.stats$normal),
                tumor_vect = mean_gc(gc.stats$tumor),
                tum_depth = weighted.mean(x = gc.stats$tumor$depth,
                    w = colSums(gc.stats$tumor$n)),
                nor_depth = weighted.mean(x = gc.stats$normal$depth,
                    w = colSums(gc.stats$normal$n)))
        } else {
            list(normal_vect = median_gc(gc.stats$normal),
                tumor_vect = median_gc(gc.stats$tumor),
                tum_depth = weighted.median(x = gc.stats$tumor$depth,
                    w = colSums(gc.stats$tumor$n)),
                nor_depth = weighted.median(x = gc.stats$normal$depth,
                    w = colSums(gc.stats$normal$n)))
        }
    }, error = function(e) {
        stop("GC content processing failed: ", e$message)
    })
}

# check the XY parameter lists chromosomes included in the
# input file and if the difference is the 'chr' prefix, add
# a warning and modify the parameter. If the overal is no
# chromosome, add a waring that no chromosome overal where
# found

check_XY <- function(xy, chr_vector) {
    # Validate inputs
    if (!is.vector(xy) || !is.vector(chr_vector)) {
        stop("Both xy and chr_vector must be vectors")
    }

    xy_names <- names(xy)

    # Find direct matches
    in_input <- xy[xy %in% chr_vector]

    # If no direct matches, check for chr prefix
    # differences
    if (length(in_input) == 0) {
        # Check if chr_vector has 'chr' prefix while xy
        # doesn't
        chr_stripped <- gsub("^chr", "", chr_vector)
        xy_without_chr <- gsub("chr", "", xy)
        xy_with_chr <- paste0("chr", xy)

        if (any(xy %in% chr_stripped)) {
            warning("Adding 'chr' prefix to xy parameter to match input chromosomes")
            xy <- xy_with_chr
        } else if (any(xy_without_chr %in% chr_vector)) {
            warning("Removing 'chr' prefix to xy parameter to match input chromosomes")
            xy <- xy_without_chr
        } else {
            warning("No chromosome matches found between xy and input vector")
        }
    }
    names(xy) <- xy_names

    # Return modified xy parameter
    return(xy)
}

process_depths <- function(seqz.data, gc_splines, avg_depths,
    ignore.normal) {
    if (nrow(seqz.data) == 0) {
        return(list(tumor = numeric(0), normal = numeric(0),
            ratio = numeric(0), norm_gc_stats = NULL, seqz.data = seqz.data))
    }

    # Calculate normalized depths
    tumor_depth <- seqz.data$depth.tumor/predict(gc_splines$tumor,
        seqz.data$GC.percent)$y
    normal_depth <- seqz.data$depth.normal/predict(gc_splines$normal,
        seqz.data$GC.percent)$y

    norm_gc_stats <- depths_gc(depth_n = round(normal_depth *
        avg_depths$normal, 0), depth_t = round(tumor_depth *
        avg_depths$tumor, 0), gc = seqz.data$GC.percent)

    # Calculate ratio and add it to seqz.data
    ratio <- if (ignore.normal) {
        round(tumor_depth, 3)
    } else {
        round(tumor_depth/normal_depth, 3)
    }

    # Add calculated values to seqz.data
    seqz.data$adjusted.ratio <- ratio
    seqz.data$depth.ratio <- seqz.data$depth.tumor/seqz.data$depth.normal

    list(tumor = tumor_depth, normal = normal_depth, ratio = ratio,
        norm_gc_stats = norm_gc_stats, seqz.data = seqz.data  # Return modified seqz.data
)
}

# Improved window calculation with memory optimization
calculate_windows <- function(seqz.data, depths, window, overlap,
    avg_depths) {
    if (nrow(seqz.data) == 0) {
        return(list(ratio = data.frame(), normal = data.frame(),
            tumor = data.frame(), raw_ratio = data.frame(), n_normal = data.frame(),
            n_tumor = data.frame(), baf = list()))
    }

    required_cols <- c("adjusted.ratio", "depth.ratio", "position",
        "chromosome", "depth.normal")
    if (!all(required_cols %in% names(seqz.data))) {
        stop("Missing required columns in seqz.data: ", paste(setdiff(required_cols,
            names(seqz.data)), collapse = ", "))
    }

    # Precompute window indices once, reuse for all 6 variables
    precomp <- precompute_windows(seqz.data$position,
        seqz.data$chromosome, window = window, overlap = overlap)

    list(
        ratio = windowValues_fast(seqz.data$adjusted.ratio, precomp,
            weight = seqz.data$depth.normal),
        normal = windowValues_fast(seqz.data$depth.normal / avg_depths$normal, precomp),
        tumor = windowValues_fast(seqz.data$depth.tumor / avg_depths$tumor, precomp),
        raw_ratio = windowValues_fast(seqz.data$depth.ratio, precomp,
            weight = seqz.data$depth.normal),
        n_normal = windowValues_fast(depths$normal, precomp),
        n_tumor = windowValues_fast(depths$tumor, precomp),
        baf = list())
}

initialize_extract_containers <- function(chromosome.list) {
    n_chr <- length(chromosome.list)
    containers <- list(windows.baf = vector("list", n_chr), windows.ratio = vector("list",
        n_chr), windows.raw_ratio = vector("list", n_chr), windows.normal = vector("list",
        n_chr), windows.tumor = vector("list", n_chr), windows.n_normal = vector("list",
        n_chr), windows.n_tumor = vector("list", n_chr), mutation.list = vector("list",
        n_chr), segments.list = vector("list", n_chr), norm.gc.list = vector("list",
        n_chr), rank_peaks.list = vector("list", n_chr), all_segments = vector("list",
        n_chr)  # Add new field
)
    names(containers$windows.baf) <- chromosome.list
    names(containers$windows.ratio) <- chromosome.list
    names(containers$windows.raw_ratio) <- chromosome.list
    names(containers$windows.normal) <- chromosome.list
    names(containers$windows.tumor) <- chromosome.list
    names(containers$windows.n_normal) <- chromosome.list
    names(containers$windows.n_tumor) <- chromosome.list
    names(containers$mutation.list) <- chromosome.list
    names(containers$segments.list) <- chromosome.list
    names(containers$norm.gc.list) <- chromosome.list
    names(containers$rank_peaks.list) <- chromosome.list
    names(containers$all_segments) <- chromosome.list
    return(containers)
}

store_chromosome_results <- function(results, containers, chr,
    idx) {
    # Direct storage without list nesting
    containers$windows.ratio[[idx]] <- results$windows$ratio[[1]]
    containers$windows.raw_ratio[[idx]] <- results$windows$raw_ratio[[1]]
    containers$windows.normal[[idx]] <- results$windows$normal[[1]]
    containers$windows.tumor[[idx]] <- results$windows$tumor[[1]]
    containers$windows.n_normal[[idx]] <- results$windows$n_normal[[1]]
    containers$windows.n_tumor[[idx]] <- results$windows$n_tumor[[1]]
    containers$windows.baf[[idx]] <- results$windows$baf[[1]]
    containers$segments.list[[idx]] <- results$segments$seg
    containers$mutation.list[[idx]] <- results$mutations
    containers$norm.gc.list[[idx]] <- results$norm_gc_stats
    containers$rank_peaks.list[[idx]] <- list(selected_win = results$segments$selected_win,
        peak_win = results$segments$peak_win)
    containers$all_segments[[idx]] <- results$segments$breaks_list
    return(containers)
}

log_chromosome_results <- function(segments, seqz.data, mutations,
    num_het_positions) {
    message("Processed ", nrow(seqz.data), " data points.")
    message("Detected ", nrow(mutations), " mutations.")
    message("Detected ", num_het_positions, " heterozygous positions.")
}

extract_process_chromosome <- function(chr, file, gc_stats, gc_splines,
    containers, params) {
    # Memory-optimized version: reads file multiple times with
    # only needed columns, instead of holding all data in memory.
    # Works with per-chromosome seqz.gz files.
    # Prerequisite: input file must not be modified during execution.

    safely_execute({
        chr_file <- file
        # Selective-read mode: when file is a per-chromosome file
        # (not the merged file) and no custom breaks
        is_chr_file <- !is.null(params$original_file) && !identical(file, params$original_file)
        use_selective_read <- is_chr_file && is.null(params$breaks)

        if (!use_selective_read) {
            # Legacy mode: load everything at once
            file.lines <- gc_stats$file.metrics[which(params$chr.vect == chr), ]
            seqz.data <- read.seqz(file, n_lines = c(file.lines$start, file.lines$end),
                chr_name = chr)
            seqz.data <- validate_data_frame(seqz.data,
                c("depth.tumor", "depth.normal", "GC.percent"),
                "Chromosome data")
            if (is.null(seqz.data))
                return(NULL)

            depths_result <- process_depths(seqz.data, gc_splines, params$avg_depths,
                params$ignore.normal)
            seqz.data <- depths_result$seqz.data

            windows <- calculate_windows(seqz.data = seqz.data, depths = depths_result,
                window = params$window,
                overlap = params$overlap, avg_depths = params$avg_depths)

            seqz.het <- seqz.data[seqz.data$zygosity.normal == "het", ]
            num_het_positions <- nrow(seqz.het)
            if (num_het_positions > 0) {
                baf_precomp <- precompute_windows(seqz.het$position,
                    seqz.het$chromosome, window = params$window, overlap = params$overlap)
                windows$baf <- windowBf_fast(Af = seqz.het$Af, Bf = seqz.het$Bf,
                    good.reads = seqz.het$good.reads, precomp = baf_precomp, conf = 0.95)
            } else {
                windows$baf <- list(data.frame(
                    start = min(seqz.data$position, na.rm = TRUE),
                    end = max(seqz.data$position, na.rm = TRUE),
                    mean = 0, q0 = 0, q1 = 0, N = 1))
            }

            segments <- process_segments(seqz.data, params$breaks, chr, windows, params)

            mutations <- tryCatch({
                mutation.table(seqz.data, mufreq.threshold = params$min.mut.freq,
                    min.reads = params$min.reads,
                    min.reads.normal = params$min.reads.normal,
                    max.mut.types = params$max.mut.types,
                    min.type.freq = params$min.type.freq,
                    min.fw.freq = params$min.fw.freq,
                    segments = segments$seg)
            }, error = function(e) {
                message("Warning: Mutation table calculation failed: ", e$message)
                data.frame()
            })

            containers <- store_chromosome_results(
                list(windows = windows, segments = segments,
                     mutations = mutations,
                     norm_gc_stats = depths_result$norm_gc_stats),
                containers, chr, which(params$chromosome.list == chr))

            if (params$verbose) {
                log_chromosome_results(segments, seqz.data, mutations, nrow(seqz.het))
            }
            return(containers)
        }

        # ============================================
        # Selective-read mode: read only needed columns per step
        # ============================================

        if (params$verbose)
            message("\n[Selective-read] Processing chromosome ", chr)

        # --- Step 1: process_depths ---
        # Read only depth + GC columns, compute adjusted.ratio
        depth_data <- read.seqz.columns(chr_file,
            c("chromosome", "position", "depth.tumor", "depth.normal", "GC.percent"))
        depth_data <- validate_data_frame(depth_data,
            c("depth.tumor", "depth.normal", "GC.percent"),
            "Chromosome data")
        if (is.null(depth_data))
            return(NULL)

        depths_result <- process_depths(depth_data, gc_splines,
            params$avg_depths, params$ignore.normal)
        # Save computed vectors we need later
        n_rows <- nrow(depth_data)
        positions <- depth_data$position
        chromosomes <- depth_data$chromosome
        adjusted_ratio <- depths_result$seqz.data$adjusted.ratio
        depth_ratio <- depths_result$seqz.data$depth.ratio
        depth_normal <- depth_data$depth.normal
        depth_tumor <- depth_data$depth.tumor
        norm_gc_stats <- depths_result$norm_gc_stats
        norm_depths <- list(normal = depths_result$normal, tumor = depths_result$tumor)
        rm(depth_data, depths_result); gc(verbose = FALSE)

        if (params$verbose)
            message("  Step 1 (depths): done, ", n_rows, " rows")

        # --- Step 2: calculate_windows ---
        # Build minimal data.frame for calculate_windows
        win_data <- data.frame(
            chromosome = chromosomes,
            position = positions,
            adjusted.ratio = adjusted_ratio,
            depth.ratio = depth_ratio,
            depth.normal = depth_normal,
            depth.tumor = depth_tumor)

        win_depths <- list(normal = norm_depths$normal, tumor = norm_depths$tumor)

        windows <- calculate_windows(seqz.data = win_data, depths = win_depths,
            window = params$window,
            overlap = params$overlap, avg_depths = params$avg_depths)
        rm(win_data, win_depths, depth_tumor, depth_normal); gc(verbose = FALSE)

        if (params$verbose)
            message("  Step 2 (windows): done")

        # --- Step 3: BAF window ---
        # Read BAF-related columns, filter het only
        baf_data <- read.seqz.filtered(chr_file,
            c("chromosome", "position", "zygosity.normal", "Af", "Bf", "good.reads"),
            filter_fn = function(chunk) chunk$zygosity.normal == "het")

        num_het_positions <- nrow(baf_data)
        if (num_het_positions > 0) {
            baf_precomp <- precompute_windows(baf_data$position, baf_data$chromosome,
                window = params$window, overlap = params$overlap)
            windows$baf <- windowBf_fast(Af = baf_data$Af, Bf = baf_data$Bf,
                good.reads = baf_data$good.reads, precomp = baf_precomp, conf = 0.95)
        } else {
            windows$baf <- list(data.frame(
                start = min(positions, na.rm = TRUE),
                end = max(positions, na.rm = TRUE),
                mean = 0, q0 = 0, q1 = 0, N = 1))
        }
        rm(baf_data); gc(verbose = FALSE)

        if (params$verbose) {
            message("  Step 3 (BAF): done, ", num_het_positions, " het positions")
            message("  ratio entries: ", nrow(windows$ratio[[1]]))
            message("  BAF entries: ", nrow(windows$baf[[1]]))
        }

        # --- Step 4: segmentation ---
        # Read columns needed by slide_tracks + segment.breaks
        seg_data <- read.seqz.columns(chr_file,
            c("chromosome", "position", "zygosity.normal",
              "good.reads", "Af", "Bf", "depth.normal"))
        # Add adjusted.ratio (computed in Step 1, not in file)
        if (nrow(seg_data) != length(adjusted_ratio))
            stop("Row count mismatch between seg_data and adjusted_ratio")
        seg_data$adjusted.ratio <- adjusted_ratio

        # process_segments will call slide_tracks + segment.breaks
        segments <- process_segments(seg_data, params$breaks, chr, windows, params)
        rm(seg_data); gc(verbose = FALSE)

        if (params$verbose)
            message("  Step 4 (segmentation): done")

        # --- Step 5: mutation.table ---
        # Read with aggressive filtering
        mut_min_reads <- params$min.reads
        mut_min_reads_normal <- params$min.reads.normal
        mut_mufreq_threshold <- params$min.mut.freq

        mut_data <- read.seqz.filtered(chr_file,
            c("chromosome", "position", "zygosity.normal",
              "AB.tumor", "AB.normal", "good.reads", "depth.normal",
              "Af", "tumor.strand", "GC.percent"),
            filter_fn = function(chunk) {
                chunk$zygosity.normal == "hom" &
                chunk$AB.tumor != "." &
                chunk$good.reads >= mut_min_reads &
                chunk$depth.normal >= mut_min_reads_normal &
                chunk$Af <= (1 - mut_mufreq_threshold)
            })

        # Match adjusted.ratio from Step 1 by position
        if (nrow(mut_data) > 0) {
            idx <- match(mut_data$position, positions)
            valid <- !is.na(idx)
            if (!all(valid)) {
                warning(sum(!valid), " mutation positions not found in depth data")
                mut_data <- mut_data[valid, , drop = FALSE]
                idx <- idx[valid]
            }
            mut_data$adjusted.ratio <- adjusted_ratio[idx]
        }
        rm(adjusted_ratio); gc(verbose = FALSE)

        mutations <- tryCatch({
            mutation.table(mut_data,
                mufreq.threshold = params$min.mut.freq,
                min.reads = params$min.reads,
                min.reads.normal = params$min.reads.normal,
                max.mut.types = params$max.mut.types,
                min.type.freq = params$min.type.freq,
                min.fw.freq = params$min.fw.freq,
                segments = segments$seg)
        }, error = function(e) {
            message("Warning: Mutation table calculation failed: ", e$message)
            data.frame()
        })
        rm(mut_data); gc(verbose = FALSE)

        if (params$verbose)
            message("  Step 5 (mutations): done, ", nrow(mutations), " mutations")

        # --- Store results ---
        containers <- store_chromosome_results(
            list(windows = windows, segments = segments,
                 mutations = mutations, norm_gc_stats = norm_gc_stats),
            containers, chr, which(params$chromosome.list == chr))

        if (params$verbose) {
            log_chromosome_results(segments,
                data.frame(position = positions, chromosome = chromosomes),
                mutations, num_het_positions)
        }

        containers
    }, NULL, sprintf("Processing chromosome %s", chr))
}

extract_initialize_parameters <- function(file, window, overlap = 1,
    slide_win = 100, peak_wins = 2^(1:10) * 10, support_threshold = 0.2,
    normalization.method = "mean", ignore.normal = FALSE, verbose = TRUE,
    chromosome.list = NULL, breaks = NULL, min.mut.freq = 0.1,
    min.reads = 40, min.reads.normal = 10, min.reads.baf = 1,
    max.mut.types = 1, min.type.freq = 0.9, min.fw.freq = 0,
    assembly = "hg38", female = TRUE, XY = c(X = "X", Y = "Y"),
    gc.stats = NULL, do_raster = FALSE, smooth_gc = FALSE, min_times_gc = 5,
    gc_grid = 250, parallel = 1, weighted.mean = TRUE, cytoband_file = NULL, ...) {
    # Initialize GC stats if needed
    local_gc_stats <- if (is.null(gc.stats)) {
        gc.sample.stats(file, verbose = verbose, parallel = parallel,
            smooth = smooth_gc, min_times = min_times_gc, cl = parallel,
            n = gc_grid)
    } else {
        gc.stats
    }

    # Get chromosome vector
    chr.vect <- as.character(local_gc_stats$file.metrics$chr)

    # Initialize chromosome list if needed
    if (is.null(chromosome.list)) {
        chromosome.list <- select_chromosomes_with_centromere(assembly,
            chr.vect, cytoband_file = cytoband_file)
    } else {
        chromosome.list <- chromosome.list[chromosome.list %in%
            chr.vect]
    }

    # Check the XY parameter
    XY_adjusted <- check_XY(XY, chr.vect)

    # Return complete parameter list
    list(original_file = file, window = window, overlap = overlap,
        slide_win = slide_win,
        peak_wins = peak_wins, support_threshold = support_threshold,
        normalization.method = normalization.method, ignore.normal = ignore.normal,
        verbose = verbose, assembly = assembly, chromosome.list = chromosome.list,
        female = female, XY = XY_adjusted, breaks = if (is.null(dim(breaks))) NULL else breaks,
        chr.vect = chr.vect, min.mut.freq = min.mut.freq, min.reads = min.reads,
        min.reads.normal = min.reads.normal, min.reads.baf = min.reads.baf,
        max.mut.types = max.mut.types, min.type.freq = min.type.freq,
        min.fw.freq = min.fw.freq, gc.stats = local_gc_stats,
        do_raster = do_raster, smooth_gc = smooth_gc, min_times_gc = min_times_gc,
        gc_grid = gc_grid, parallel = parallel, ratio_baf_raster = data.frame(dr = NULL),
        avg_depths = list(normal = weighted.mean(local_gc_stats$normal$depth,
            colSums(local_gc_stats$normal$n)), tumor = weighted.mean(local_gc_stats$tumor$depth,
            colSums(local_gc_stats$tumor$n))), weighted.mean = weighted.mean,
            cytoband_file = cytoband_file)
}

finalize_extract_results <- function(containers, params, gc_stats,
    gc_data) {
    # Name lists
    for (list_name in names(containers)) {
        if (length(containers[[list_name]]) == length(params$chromosome.list)) {
            names(containers[[list_name]]) <- params$chromosome.list
        }
    }

    # Process GC normalization
    gc_norm <- unfold_gc(cbind(unique = 0, lines = 0, do.call(rbind,
        containers$norm.gc.list)), stats = FALSE, smooth = params$smooth_gc,
        min_times = params$min_times_gc, cl = params$parallel,
        grid_size = params$gc_grid)

    # Calculate final depths
    avg_tum_ndepth <- weighted.mean(x = gc_norm$tumor$depth,
        w = colSums(gc_norm$tumor$n))
    avg_nor_ndepth <- weighted.mean(x = gc_norm$normal$depth,
        w = colSums(gc_norm$normal$n))

    avg_depth_ratio <- if (params$ignore.normal) {
        avg_tum_ndepth/gc_data$tum_depth
    } else {
        (avg_tum_ndepth/gc_data$tum_depth)/(avg_nor_ndepth/gc_data$nor_depth)
    }

    depths <- list(avg_depth_ratio = avg_depth_ratio, avg_tum_depth = gc_data$tum_depth,
        avg_nor_depth = gc_data$nor_depth)

    # Calculate mutation statistics
    total_mutations <- sum(sapply(containers$mutation.list, nrow))
    total_bases <- sum(sapply(containers$segments.list, function(segs) {
        sum(segs$end.pos - segs$start.pos + 1)
    }))

    # Calculate per-chromosome mutation statistics
    chr_mutation_stats <- lapply(params$chromosome.list, function(chr) {
        mutations <- nrow(containers$mutation.list[[chr]])
        bases <- sum(containers$segments.list[[chr]]$end.pos -
            containers$segments.list[[chr]]$start.pos + 1)
        list(mutations = mutations, megabases = bases/1e+06,
            mutations_per_mb = mutations/(bases/1e+06))
    })
    names(chr_mutation_stats) <- params$chromosome.list

    # Compile complete mutation statistics
    mutation_stats <- list(total_mutations = total_mutations,
        total_megabases = total_bases/1e+06, mutations_per_mb = total_mutations/(total_bases/1e+06),
        chromosomes = chr_mutation_stats)

    # Construct and return complete results
    list(BAF = containers$windows.baf, ratio = containers$windows.ratio,
        raw_ratio = containers$windows.raw_ratio,
        depths = list(
            raw = list(normal = containers$windows.normal,
                tumor = containers$windows.tumor),
            norm = list(normal = containers$windows.n_normal,
                tumor = containers$windows.n_tumor)),
        mutations = containers$mutation.list,
        segments = containers$segments.list, win_peaks = containers$rank_peaks.list,
        chromosomes = params$chromosome.list, gc = gc_stats,
        ignore.normal = params$ignore.normal, gc_norm = gc_norm,
        avg.depth.ratio = depths$avg_depth_ratio, avg.depth.tumor = depths$avg_tum_depth,
        avg.depth.normal = depths$avg_nor_depth, mutation_stats = mutation_stats,
        all_segments = containers$all_segments, gender = if (params$female) "female" else "male",
        XY = params$XY)
}

# Add error handling wrapper
safely_compare_bins <- function(start.pos, end.pos, values, windows,
    sd.values) {
    tryCatch({
        compare_bins(start.pos, end.pos, values, windows, sd.values)
    }, error = function(e) {
        message("Warning: Bin comparison failed: ", e$message)
        return(0)  # Return neutral score on failure
    })
}

# Improved rank_segments with better error handling and
# optimization
rank_segments <- function(breaks_list, windows, params) {
    if (length(breaks_list) <= 1) {
        return(list(segs = breaks_list[[1]], selected_win = params$peak_wins[1],
            peak_win = data.frame(peak_win = params$peak_wins[1],
                baf_fit = 0, ratio_fit = 0, n_segs = nrow(breaks_list[[1]]))))
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
            baf_vs_bins <- safely_compare_bins(x$start.pos, x$end.pos,
                x$Bf, baf_win, x$sd.BAF)
            ratio_vs_bins <- safely_compare_bins(x$start.pos,
                x$end.pos, x$depth.ratio, ratio_win, x$sd.ratio)
            c(baf_fit = baf_vs_bins, ratio_fit = ratio_vs_bins)
        }, cl = params$parallel)
    } else {
        lapply(breaks_list, function(x) {
            baf_vs_bins <- safely_compare_bins(x$start.pos, x$end.pos,
                x$Bf, baf_win, x$sd.BAF)
            ratio_vs_bins <- safely_compare_bins(x$start.pos,
                x$end.pos, x$depth.ratio, ratio_win, x$sd.ratio)
            c(baf_fit = baf_vs_bins, ratio_fit = ratio_vs_bins)
        })
    }

    # Create comparison dataframe
    compare_bins_segs <- data.frame(peak_win = params$peak_wins,
        do.call(rbind, compare_bins_list), n_segs = vapply(breaks_list,
            nrow, numeric(1)))

    # Subsets the compare_bins_segs with both metrics above
    # 0.75
    compare_bins_segs_high <- compare_bins_segs[apply(compare_bins_segs[,
        c("baf_fit", "ratio_fit")], 1, function(x) {
        all(x > 0.75)
    }), ]
    # Fallback to all metrics if no high-quality fits
    if (nrow(compare_bins_segs_high) == 0) {
        compare_bins_segs_high <- compare_bins_segs
    }

    # Calculate ranks with weights for different metrics
    ranks_fits <- cbind(baf = rank(-compare_bins_segs_high$baf_fit,
        ties.method = "max"), ratio = rank(-compare_bins_segs_high$ratio_fit,
        ties.method = "max"), n_segs = rank(compare_bins_segs_high$n_segs,
        ties.method = "min"))

    # Select best fit considering all metrics
    total_ranks <- rowSums(ranks_fits)
    best_fits <- which(total_ranks == min(total_ranks))
    select_win <- max(compare_bins_segs_high$peak_win[best_fits])

    # Debug information if verbose
    if (params$verbose) {
        message("Segment ranking results:")
        message("Selected window size: ", select_win)
        message("Number of segments: ",
            compare_bins_segs_high$n_segs[compare_bins_segs_high$peak_win == select_win])
    }

    list(segs = breaks_list[[as.character(select_win)]], selected_win = select_win,
        peak_win = compare_bins_segs)
}

# Update process_segments to use new rank_segments
process_segments <- function(seqz.data, breaks, chr, windows,
    params) {
    # TODO: Add support for alternative segmentation
    # methods FIXME: Current segmentation can be memory
    # intensive for large chromosomes TODO: Consider adding
    # parallel processing for break detection

    # Ensure weighted.mean has a default value if not in
    # params
    weighted.mean <- if (!is.null(params$weighted.mean)) {
        params$weighted.mean
    } else {
        TRUE  # Default value
    }

    # Handle segmentation
    if (is.null(breaks)) {
        diff_track <- slide_tracks(seqz.data, params$slide_win,
            signal_out = "both", verbose = params$verbose)

        breaks_chr_list <- lapply(params$peak_wins, function(x) {
            breaks_chr <- extract_breaks_tracks(track = diff_track,
                breaks = breaks, peak_win = x, assembly = params$assembly,
                chromosome = chr, cytoband_file = params$cytoband_file)

            if (inherits(breaks_chr, "try-error") || is.null(breaks_chr) ||
                nrow(breaks_chr) == 0 || length(breaks_chr) ==
                0) {
                breaks_chr <- data.frame(chrom = chr, start.pos = min(seqz.data$position,
                  na.rm = TRUE), end.pos = max(seqz.data$position,
                  na.rm = TRUE))
            }

            tryCatch({
                segment.breaks(seqz.tab = seqz.data, breaks = breaks_chr,
                    min.reads.baf = params$min.reads.baf,
                    weighted.mean = weighted.mean)  # Use local variable
            }, error = function(e) {
                message("Warning: Segment calculation failed: ",
                  e$message)
                data.frame(chrom = chr, start.pos = min(seqz.data$position,
                  na.rm = TRUE), end.pos = max(seqz.data$position,
                  na.rm = TRUE), Bf = 0, depth.ratio = mean(seqz.data$depth.ratio,
                  na.rm = TRUE))
            })
        })

        names(breaks_chr_list) <- as.character(params$peak_wins)
        # Return segment rank and information with
        # breaks_list
        segment_results <- rank_segments(breaks_chr_list, windows,
            params)
        # breakpoints_consensus <-
        # debug_evaluate_consensus(breaks_chr_list,
        # windows, params$support_threshold)

        # newbreaks <- as.data.frame(do.call(rbind,
        # lapply(breakpoints_consensus$supported_breaks,
        # unlist))) newbreaks <-
        # newbreaks[order(newbreaks$position), ]

        # breaks_chr <- cbind(chrom=chr,
        # position_to_breaks(newbreaks$position)) segs <-
        # segment.breaks(seqz.tab = seqz.data, breaks =
        # breaks_chr, min.reads.baf = params$min.reads.baf,
        # weighted.mean = weighted.mean)
        return(list(seg = segment_results$segs, breaks_list = breaks_chr_list,
            selected_win = segment_results$selected_win, peak_win = segment_results$peak_win))
    } else {
        breaks_chr <- breaks[breaks$chrom == chr, ]
        segs <- segment.breaks(seqz.tab = seqz.data, breaks = breaks_chr,
            min.reads.baf = params$min.reads.baf, weighted.mean = weighted.mean)
        select_win <- 0
        compare_bins_segs <- data.frame(peak_win = 0, baf_fit = compare_bins(segs$start.pos,
            segs$end.pos, segs$Bf, windows$baf[[1]], segs$sd.BAF),
            ratio_fit = compare_bins(segs$start.pos, segs$end.pos,
                windows$ratio[[1]], segs$sd.ratio),
            n_segs = nrow(segs))

        # Create single-element breaks_list for
        # user-provided breaks
        breaks_chr_list <- list(user = segs)

        return(list(seg = segs, breaks_list = breaks_chr_list,
            selected_win = select_win, peak_win = compare_bins_segs))
    }
}

# Add safer parallel processing management function
manage_parallel_cluster <- function(n_cores, type = "FORK") {
    if (n_cores > 1) {
        cl <- NULL
        tryCatch({
            if (type == "PSOCK") {
                cl <- parallel::makeCluster(n_cores, type = "PSOCK")
                parallel::clusterEvalQ(cl, library(sequenza))
            } else {
                if (.Platform$OS.type != "unix")
                    stop("FORK clusters require Unix (use type='PSOCK' on Windows)")
                cl <- parallel::makeForkCluster(n_cores)
            }
            return(cl)
        }, error = function(e) {
            warning("Failed to create ", type, " cluster: ", e$message,
                "\nFalling back to sequential processing")
            if (!is.null(cl))
                try(parallel::stopCluster(cl), silent = TRUE)
            NULL
        })
    }
    return(NULL)
}
