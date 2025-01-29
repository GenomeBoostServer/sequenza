# Helper function to validate input parameters
validate_params <- function(params) {
    required <- c("window", "overlap", "normalization.method")
    missing <- required[!required %in% names(params)]
    if (length(missing) > 0) {
        stop("Missing required parameters: ", paste(missing, collapse = ", "))
    }

    if (!params$normalization.method %in% c("mean", "median")) {
        stop("normalization.method must be either 'mean' or 'median'")
    }
}

process_gc_content <- function(gc.stats, normalization.method) {
    if (!is.list(gc.stats) || !all(c("normal", "tumor") %in% names(gc.stats))) {
        stop("Invalid gc.stats format")
    }

    tryCatch({
        if (normalization.method == "mean") {
            list(normal_vect = mean_gc(gc.stats$normal), tumor_vect = mean_gc(gc.stats$tumor),
                tum_depth = weighted.mean(x = gc.stats$tumor$depth, w = colSums(gc.stats$tumor$n)),
                nor_depth = weighted.mean(x = gc.stats$normal$depth, w = colSums(gc.stats$normal$n)))
        } else {
            list(normal_vect = median_gc(gc.stats$normal), tumor_vect = median_gc(gc.stats$tumor),
                tum_depth = weighted.median(x = gc.stats$tumor$depth, w = colSums(gc.stats$tumor$n)),
                nor_depth = weighted.median(x = gc.stats$normal$depth, w = colSums(gc.stats$normal$n)))
        }
    }, error = function(e) {
        stop("GC content processing failed: ", e$message)
    })
}

process_depths <- function(seqz.data, gc_splines, avg_depths, ignore.normal) {
    if (nrow(seqz.data) == 0) {
        return(list(tumor = numeric(0), normal = numeric(0), ratio = numeric(0),
            norm_gc_stats = NULL, seqz.data = seqz.data))
    }

    # Calculate normalized depths
    tumor_depth <- seqz.data$depth.tumor/predict(gc_splines$tumor, seqz.data$GC.percent)$y
    normal_depth <- seqz.data$depth.normal/predict(gc_splines$normal, seqz.data$GC.percent)$y

    norm_gc_stats <- depths_gc(depth_n = round(normal_depth * avg_depths$normal,
        0), depth_t = round(tumor_depth * avg_depths$tumor, 0), gc = seqz.data$GC.percent)

    # Calculate ratio and add it to seqz.data
    ratio <- if (ignore.normal) {
        round(tumor_depth, 3)
    } else {
        round(tumor_depth/normal_depth, 3)
    }

    # Add calculated values to seqz.data
    seqz.data$adjusted.ratio <- ratio
    seqz.data$depth.ratio <- seqz.data$depth.tumor/seqz.data$depth.normal

    list(tumor = tumor_depth, normal = normal_depth, ratio = ratio, norm_gc_stats = norm_gc_stats,
        seqz.data = seqz.data  # Return modified seqz.data
)
}

# Improved window calculation with memory optimization
calculate_windows <- function(seqz.data, depths, window, overlap, avg_depths) {
    if (nrow(seqz.data) == 0) {
        return(list(ratio = data.frame(), normal = data.frame(), tumor = data.frame(),
            raw_ratio = data.frame(), n_normal = data.frame(), n_tumor = data.frame(),
            baf = list()))
    }

    required_cols <- c("adjusted.ratio", "depth.ratio", "position", "chromosome",
        "depth.normal")
    if (!all(required_cols %in% names(seqz.data))) {
        stop("Missing required columns in seqz.data: ", paste(setdiff(required_cols,
            names(seqz.data)), collapse = ", "))
    }

    # Calculate all window values at once
    list(ratio = windowValues(x = seqz.data$adjusted.ratio, positions = seqz.data$position,
        chromosomes = seqz.data$chromosome, window = window, overlap = overlap, weight = seqz.data$depth.normal),
        normal = windowValues(x = seqz.data$depth.normal/avg_depths$normal, positions = seqz.data$position,
            chromosomes = seqz.data$chromosome, window = window, overlap = overlap),
        tumor = windowValues(x = seqz.data$depth.tumor/avg_depths$tumor, positions = seqz.data$position,
            chromosomes = seqz.data$chromosome, window = window, overlap = overlap),
        raw_ratio = windowValues(x = seqz.data$depth.ratio, positions = seqz.data$position,
            chromosomes = seqz.data$chromosome, window = window, overlap = overlap,
            weight = seqz.data$depth.normal), n_normal = windowValues(x = depths$normal,
            positions = seqz.data$position, chromosomes = seqz.data$chromosome, window = window,
            overlap = overlap), n_tumor = windowValues(x = depths$tumor, positions = seqz.data$position,
            chromosomes = seqz.data$chromosome, window = window, overlap = overlap),
        baf = list())
}

initialize_extract_containers <- function(chromosome.list) {
    n_chr <- length(chromosome.list)
    containers <- list(windows.baf = vector("list", n_chr), windows.ratio = vector("list",
        n_chr), windows.raw_ratio = vector("list", n_chr), windows.normal = vector("list",
        n_chr), windows.tumor = vector("list", n_chr), windows.n_normal = vector("list",
        n_chr), windows.n_tumor = vector("list", n_chr), mutation.list = vector("list",
        n_chr), segments.list = vector("list", n_chr), norm.gc.list = vector("list",
        n_chr), rank_peaks.list = vector("list", n_chr)  # Add new field
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
    return(containers)
}

store_chromosome_results <- function(results, containers, chr, idx) {
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
    return(containers)
}

log_chromosome_results <- function(segments, seqz.data, mutations, num_het_positions) {
    message("Processed ", nrow(seqz.data), " data points.")
    message("Detected ", nrow(mutations), " mutations.")
    message("Detected ", num_het_positions, " heterozygous positions.")
}

process_single_chromosome <- function(chr, file, gc_stats, gc_splines, containers,
    params) {

    if (params$verbose) {
        message("Processing chromosome ", chr)
    }

    # Read chromosome data
    file.lines <- gc_stats$file.metrics[which(params$chr.vect == chr), ]
    seqz.data <- read.seqz(file, n_lines = c(file.lines$start, file.lines$end), chr_name = chr)

    # Process depths and get modified seqz.data
    depths_result <- process_depths(seqz.data, gc_splines, params$avg_depths, params$ignore.normal)
    seqz.data <- depths_result$seqz.data  # Use updated seqz.data

    # Calculate windows with initialized data
    windows <- calculate_windows(seqz.data = seqz.data, depths = depths_result, window = params$window,
        overlap = params$overlap, avg_depths = params$avg_depths)

    # Process BAF if heterozygous positions exist
    seqz.het <- seqz.data[seqz.data$zygosity.normal == "het", ]
    num_het_positions <- nrow(seqz.het)
    if (num_het_positions > 0) {
        windows$baf <- windowBf(Af = seqz.het$Af, Bf = seqz.het$Bf, good.reads = seqz.het$good.reads,
            chromosomes = seqz.het$chromosome, positions = seqz.het$position, conf = 0.95,
            window = params$window, overlap = params$overlap)
    } else {
        windows$baf <- list(data.frame(start = min(seqz.data$position, na.rm = TRUE),
            end = max(seqz.data$position, na.rm = TRUE), mean = 0, q0 = 0, q1 = 0,
            N = 1))
    }

    # Add debug message
    if (params$verbose) {
        message("Windows calculation results for chr ", chr, ":")
        message("  ratio entries: ", nrow(windows$ratio[[1]]))
        message("  raw_ratio entries: ", nrow(windows$raw_ratio[[1]]))
        message("  BAF entries: ", nrow(windows$baf[[1]]))
    }

    # Process segments
    segments <- process_segments(seqz.data, params$breaks, chr, windows, params)

    # Ensure seqz.data has necessary columns for mutation.table
    required_columns <- c("good.reads", "depth.normal")
    if (!all(required_columns %in% colnames(seqz.data))) {
        stop("seqz.data is missing required columns for mutation.table")
    }

    # Process mutations using mutation.table

    mutations <- tryCatch({
        mutation.table(seqz.data, mufreq.threshold = params$min.mut.freq, min.reads = params$min.reads,
            min.reads.normal = params$min.reads.normal, max.mut.types = params$max.mut.types,
            min.type.freq = params$min.type.freq, min.fw.freq = params$min.fw.freq,
            segments = segments$seg)
    }, error = function(e) {
        message("Warning: Mutation table calculation failed: ", e$message)
        data.frame()  # Return an empty data frame on error
    })


    # Store results
    containers <- store_chromosome_results(list(windows = windows, segments = segments,
        mutations = mutations, norm_gc_stats = depths_result$norm_gc_stats), containers,
        chr, which(params$chromosome.list == chr))

    if (params$verbose) {
        log_chromosome_results(segments, seqz.data, mutations, num_het_positions)
    }

    containers
}

# Add proper documentation
#' @title Extract Sequenza Data from Tumor/Normal Sequencing
#' @description Process and analyze genomic sequencing data from paired normal-tumor samples 
#' to extract copy number, mutation and cellularity information
#' @param file Path to the seqz format input file
#' @param window Size of windows for binning, in base pairs (default: 1e6)
#' @param overlap Fraction of overlap between adjacent windows (default: 1)
#' @param slide_win Window size for breakpoint detection (default: 100)
#' @param peak_wins Vector of window sizes for peak detection (default: seq(50,300,25))
#' @param normalization.method Method for GC normalization ('mean' or 'median', default: 'mean')
#' @param ignore.normal Ignore normal sample depth in ratio calculations (default: FALSE)
#' @param verbose Print detailed progress messages (default: TRUE)
#' @param chromosome.list Vector of chromosomes to analyze (default: NULL, all chromosomes)
#' @param breaks Pre-computed breakpoints data.frame (default: NULL)
#' @param min.mut.freq Minimum mutation frequency threshold (default: 0.1)
#' @param min.reads Minimum number of reads for mutation calling (default: 40)
#' @param min.reads.normal Minimum reads in normal sample (default: 10)
#' @param min.reads.baf Minimum reads for BAF calculation (default: 1)
#' @param max.mut.types Maximum mutation types per position (default: 1)
#' @param min.type.freq Minimum frequency for mutation type (default: 0.9)
#' @param min.fw.freq Minimum forward strand frequency (default: 0)
#' @param assembly Genome assembly version (default: 'hg38')
#' @param gc.stats Pre-computed GC statistics (default: NULL)
#' @param do_raster Generate raster plots (default: FALSE)
#' @param smooth_gc Smooth GC content calculations (default: FALSE)
#' @param min_times_gc Minimum times for GC normalization (default: 5)
#' @param gc_grid Grid size for GC binning (default: 250)
#' @param parallel Number of parallel processes (default: 1)
#' @param weighted.mean Use weighted means in calculations (default: TRUE)
#' @param ... Additional parameters passed to internal functions
#' @return A list containing:
#' \itemize{
#'   \item BAF - B-allele frequency windows for each chromosome
#'   \item ratio - Depth ratio windows normalized by GC content
#'   \item raw_ratio - Raw depth ratio windows without normalization
#'   \item depths - List of depth information:
#'     \itemize{
#'       \item raw - Raw depth windows for normal and tumor
#'       \item norm - GC-normalized depth windows
#'     }
#'   \item mutations - List of detected mutations per chromosome
#'   \item segments - List of segmentation results per chromosome
#'   \item win_peaks - Window peak analysis results
#'   \item chromosomes - Vector of processed chromosomes
#'   \item gc - GC content statistics
#'   \item gc_norm - GC normalization results
#'   \item avg.depth.ratio - Average depth ratio
#'   \item avg.depth.tumor - Average tumor depth
#'   \item avg.depth.normal - Average normal depth
#' }
#' @examples
#' \dontrun{
#' # Basic usage
#' results <- sequenza.extract('sample.seqz', window=1e6)
#' 
#' # With custom parameters
#' results <- sequenza.extract('sample.seqz',
#'                           window=5e5,
#'                           chromosome.list=c('chr1','chr2'),
#'                           parallel=4)
#' }
#' @references
#' Favero, F. et al. Sequenza: allele-specific copy number and mutation profiles 
#' from tumor sequencing data. Ann Oncol. 2015
#' @seealso
#' \code{\link{sequenza.fit}} for fitting copy number and cellularity models
#' @export
sequenza.extract <- function(file, window = 1e+06, overlap = 1, slide_win = 100,
    peak_wins = seq(from = 50, to = 300, by = 25), normalization.method = "mean",
    ignore.normal = FALSE, verbose = TRUE, chromosome.list = NULL, breaks = NULL,
    min.mut.freq = 0.1, min.reads = 40, min.reads.normal = 10, min.reads.baf = 1,
    max.mut.types = 1, min.type.freq = 0.9, min.fw.freq = 0, assembly = "hg38", gc.stats = NULL,
    do_raster = FALSE, smooth_gc = FALSE, min_times_gc = 5, gc_grid = 250, parallel = 1,
    weighted.mean = TRUE, ...) {

    # Initialize parameters with all arguments
    params <- initialize_extract_parameters(file = file, window = window, overlap = overlap,
        slide_win = slide_win, peak_wins = peak_wins, normalization.method = normalization.method,
        ignore.normal = ignore.normal, verbose = verbose, chromosome.list = chromosome.list,
        breaks = breaks, assembly = assembly, gc.stats = gc.stats, do_raster = do_raster,
        smooth_gc = smooth_gc, min_times_gc = min_times_gc, gc_grid = gc_grid, parallel = parallel,
        weighted.mean = weighted.mean, ...)

    # Validate input parameters
    validate_params(params)

    # Initialize logging
    log_file <- paste0(tempdir(), "/sequenza_", format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".log")
    log_conn <- file(log_file, "w")
    on.exit(close(log_conn))

    tryCatch({
        # Process GC content
        gc_data <- process_gc_content(params$gc.stats, params$normalization.method)

        gc_splines <- list(normal = smooth.spline(data.frame(gc = as.numeric(names(gc_data$normal_vect)),
            depth = gc_data$normal_vect)), tumor = smooth.spline(data.frame(gc = as.numeric(names(gc_data$tumor_vect)),
            depth = gc_data$tumor_vect)))

        # Initialize containers
        containers <- initialize_extract_containers(params$chromosome.list)

        # Process each chromosome
        if (params$parallel > 1) {
            containers <- initialize_extract_containers(params$chromosome.list)
            results <- parallel::mclapply(seq_along(params$chromosome.list), function(idx) {
                chr <- params$chromosome.list[idx]
                local_containers <- initialize_extract_containers(params$chromosome.list)
                result <- process_single_chromosome(chr, file, params$gc.stats, gc_splines,
                  local_containers, params)
                result  # Return the complete result
            }, mc.cores = params$parallel)

            # Merge results back into main containers
            for (idx in seq_along(results)) {
                chr <- params$chromosome.list[idx]
                containers <- store_chromosome_results(list(windows = list(ratio = results[[idx]]$windows.ratio[[idx]],
                  raw_ratio = results[[idx]]$windows.raw_ratio[[idx]], normal = results[[idx]]$windows.normal[[idx]],
                  tumor = results[[idx]]$windows.tumor[[idx]], n_normal = results[[idx]]$windows.n_normal[[idx]],
                  n_tumor = results[[idx]]$windows.n_tumor[[idx]], baf = results[[idx]]$windows.baf[[idx]]),
                  segments = list(seg = results[[idx]]$segments.list[[idx]]), mutations = results[[idx]]$mutation.list[[idx]],
                  norm_gc_stats = results[[idx]]$norm.gc.list[[idx]]), containers,
                  chr, idx)
            }
        } else {
            for (chr in params$chromosome.list) {
                containers <- process_single_chromosome(chr, file, params$gc.stats,
                  gc_splines, containers, params)
            }
        }

        # Finalize and return results
        final_results <- finalize_extract_results(containers, params, params$gc.stats,
            gc_data)
        gc()  # Force garbage collection

        return(final_results)

    }, error = function(e) {
        message("Error in sequenza.extract: ", e$message)
        write(paste("Error:", e$message), log_conn)
        stop(e)
    }, finally = {
        # Cleanup
        gc()
    })
}

initialize_extract_parameters <- function(file, window, overlap = 1, slide_win = 100,
    peak_wins = seq(from = 50, to = 300, by = 25), normalization.method = "mean",
    ignore.normal = FALSE, verbose = TRUE, chromosome.list = NULL, breaks = NULL,
    min.mut.freq = 0.1, min.reads = 40, min.reads.normal = 10, min.reads.baf = 1,
    max.mut.types = 1, min.type.freq = 0.9, min.fw.freq = 0, assembly = "hg38", gc.stats = NULL,
    do_raster = FALSE, smooth_gc = FALSE, min_times_gc = 5, gc_grid = 250, parallel = 1,
    weighted.mean = TRUE, ...) {
    # Initialize GC stats if needed
    local_gc_stats <- if (is.null(gc.stats)) {
        gc.sample.stats(file, verbose = verbose, parallel = parallel, smooth = smooth_gc,
            min_times = min_times_gc, cl = parallel, n = gc_grid)
    } else {
        gc.stats
    }

    # Get chromosome vector
    chr.vect <- as.character(local_gc_stats$file.metrics$chr)

    # Initialize chromosome list if needed
    if (is.null(chromosome.list)) {
        chromosome.list <- select_chromosomes_with_centromere(assembly, chr.vect)
    } else {
        chromosome.list <- chromosome.list[chromosome.list %in% chr.vect]
    }

    # Return complete parameter list
    list(window = window, overlap = overlap, slide_win = slide_win, peak_wins = peak_wins,
        normalization.method = normalization.method, ignore.normal = ignore.normal,
        verbose = verbose, assembly = assembly, chromosome.list = chromosome.list,
        breaks = if (is.null(dim(breaks))) NULL else breaks, chr.vect = chr.vect,
        min.mut.freq = min.mut.freq, min.reads = min.reads, min.reads.normal = min.reads.normal,
        min.reads.baf = min.reads.baf, max.mut.types = max.mut.types, min.type.freq = min.type.freq,
        min.fw.freq = min.fw.freq, gc.stats = local_gc_stats, do_raster = do_raster,
        smooth_gc = smooth_gc, min_times_gc = min_times_gc, gc_grid = gc_grid, parallel = parallel,
        ratio_baf_raster = data.frame(dr = NULL), avg_depths = list(normal = weighted.mean(local_gc_stats$normal$depth,
            colSums(local_gc_stats$normal$n)), tumor = weighted.mean(local_gc_stats$tumor$depth,
            colSums(local_gc_stats$tumor$n))), weighted.mean = weighted.mean)
}

finalize_extract_results <- function(containers, params, gc_stats, gc_data) {
    # Name lists
    for (list_name in names(containers)) {
        if (length(containers[[list_name]]) == length(params$chromosome.list)) {
            names(containers[[list_name]]) <- params$chromosome.list
        }
    }

    # Process GC normalization
    gc_norm <- unfold_gc(cbind(unique = 0, lines = 0, do.call(rbind, containers$norm.gc.list)),
        stats = FALSE, smooth = params$smooth_gc, min_times = params$min_times_gc,
        cl = params$parallel, grid_size = params$gc_grid)

    # Calculate final depths
    avg_tum_ndepth <- weighted.mean(x = gc_norm$tumor$depth, w = colSums(gc_norm$tumor$n))
    avg_nor_ndepth <- weighted.mean(x = gc_norm$normal$depth, w = colSums(gc_norm$normal$n))

    avg_depth_ratio <- if (params$ignore.normal) {
        avg_tum_ndepth/gc_data$tum_depth
    } else {
        (avg_tum_ndepth/gc_data$tum_depth)/(avg_nor_ndepth/gc_data$nor_depth)
    }

    depths <- list(avg_depth_ratio = avg_depth_ratio, avg_tum_depth = gc_data$tum_depth,
        avg_nor_depth = gc_data$nor_depth)

    # Construct return value
    list(BAF = containers$windows.baf, ratio = containers$windows.ratio, raw_ratio = containers$windows.raw_ratio,
        depths = list(raw = list(normal = containers$windows.normal, tumor = containers$windows.tumor),
            norm = list(normal = containers$windows.n_normal, tumor = containers$windows.n_tumor)),
        mutations = containers$mutation.list, segments = containers$segments.list,
        win_peaks = containers$rank_peaks.list, chromosomes = params$chromosome.list,
        gc = gc_stats, gc_norm = gc_norm, avg.depth.ratio = depths$avg_depth_ratio,
        avg.depth.tumor = depths$avg_tum_depth, avg.depth.normal = depths$avg_nor_depth)
}

# Add error handling wrapper
safely_compare_bins <- function(start.pos, end.pos, values, windows) {
    tryCatch({
        compare_bins(start.pos, end.pos, values, windows)
    }, error = function(e) {
        message("Warning: Bin comparison failed: ", e$message)
        return(0)  # Return neutral score on failure
    })
}

# Improved rank_segments with better error handling and optimization
rank_segments <- function(breaks_list, windows, params) {
    if (length(breaks_list) <= 1) {
        return(list(segs = breaks_list[[1]], selected_win = params$peak_wins[1],
            peak_win = data.frame(peak_win = params$peak_wins[1], baf_fit = 0, ratio_fit = 0,
                n_segs = nrow(breaks_list[[1]]))))
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
        parallel::mclapply(breaks_list, function(x) {
            baf_vs_bins <- safely_compare_bins(x$start.pos, x$end.pos, x$Bf, baf_win)
            ratio_vs_bins <- safely_compare_bins(x$start.pos, x$end.pos, x$depth.ratio,
                ratio_win)
            c(baf_fit = baf_vs_bins, ratio_fit = ratio_vs_bins)
        }, mc.cores = params$parallel)
    } else {
        lapply(breaks_list, function(x) {
            baf_vs_bins <- safely_compare_bins(x$start.pos, x$end.pos, x$Bf, baf_win)
            ratio_vs_bins <- safely_compare_bins(x$start.pos, x$end.pos, x$depth.ratio,
                ratio_win)
            c(baf_fit = baf_vs_bins, ratio_fit = ratio_vs_bins)
        })
    }

    # Create comparison dataframe
    compare_bins_segs <- data.frame(peak_win = params$peak_wins, do.call(rbind, compare_bins_list),
        n_segs = vapply(breaks_list, nrow, numeric(1)))

    # Calculate ranks with weights for different metrics
    ranks_fits <- cbind(baf = rank(-compare_bins_segs$baf_fit, ties.method = "max"),
        ratio = rank(-compare_bins_segs$ratio_fit, ties.method = "max"), n_segs = rank(compare_bins_segs$n_segs,
            ties.method = "min"))

    # Select best fit considering all metrics
    total_ranks <- rowSums(ranks_fits)
    best_fits <- which(total_ranks == min(total_ranks))
    select_win <- max(compare_bins_segs$peak_win[best_fits])

    # Debug information if verbose
    if (params$verbose) {
        message("Segment ranking results:")
        message("Selected window size: ", select_win)
        message("Number of segments: ", compare_bins_segs$n_segs[compare_bins_segs$peak_win ==
            select_win])
    }

    list(segs = breaks_list[[as.character(select_win)]], selected_win = select_win,
        peak_win = compare_bins_segs)
}

# Update process_segments to use new rank_segments
process_segments <- function(seqz.data, breaks, chr, windows, params) {
    # Ensure weighted.mean has a default value if not in params
    weighted.mean <- if (!is.null(params$weighted.mean)) {
        params$weighted.mean
    } else {
        TRUE  # Default value
    }

    # Handle segmentation
    if (is.null(breaks)) {
        diff_track <- slide_tracks(seqz.data, params$slide_win, signal_out = "both",
            verbose = params$verbose)

        breaks_chr_list <- lapply(params$peak_wins, function(x) {
            breaks_chr <- extract_breaks_tracks(track = diff_track, breaks = breaks,
                peak_win = x, assembly = params$assembly, chromosome = chr)

            if (inherits(breaks_chr, "try-error") || is.null(breaks_chr) || nrow(breaks_chr) ==
                0 || length(breaks_chr) == 0) {
                breaks_chr <- data.frame(chrom = chr, start.pos = min(seqz.data$position,
                  na.rm = TRUE), end.pos = max(seqz.data$position, na.rm = TRUE))
            }

            tryCatch({
                segment.breaks(seqz.tab = seqz.data, breaks = breaks_chr, min.reads.baf = params$min.reads.baf,
                  weighted.mean = weighted.mean  # Use local variable
)
            }, error = function(e) {
                message("Warning: Segment calculation failed: ", e$message)
                data.frame(chrom = chr, start.pos = min(seqz.data$position, na.rm = TRUE),
                  end.pos = max(seqz.data$position, na.rm = TRUE), Bf = 0, depth.ratio = mean(seqz.data$depth.ratio,
                    na.rm = TRUE))
            })
        })

        names(breaks_chr_list) <- as.character(params$peak_wins)
        # Return segment rank and information
        segment_results <- rank_segments(breaks_chr_list, windows, params)
        return(list(seg = segment_results$segs, breaks_list = breaks_chr_list, selected_win = segment_results$selected_win,
            peak_win = segment_results$peak_win))

    } else {
        segs <- segment.breaks(seqz.tab = seqz.data, breaks = breaks, min.reads.baf = params$min.reads.baf,
            weighted.mean = params$weighted.mean)
        select_win <- 0
        compare_bins_segs <- data.frame(peak_win = 0, baf_fit = compare_bins(segs$start.pos,
            segs$end.pos, segs$Bf, seqz.b.win[[chr]]), ratio_fit = compare_bins(segs$start.pos,
            segs$end.pos, segs$depth.ratio, seqz.r.win[[chr]]), n_segs = nrow(segs))

        list(segs = segs, selected_win = select_win, peak_win = compare_bins_segs)
    }
}
