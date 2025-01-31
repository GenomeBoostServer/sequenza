#' @rdname sequenza
#' @export
#' @import segments
sequenza.extract <- function(file, window = 1e+06, overlap = 1, slide_win = 100,
                             peak_wins = seq(from = 50, to = 1000, by = 75), normalization.method = "mean",
                             ignore.normal = FALSE, verbose = TRUE, chromosome.list = NULL, breaks = NULL,
                             min.mut.freq = 0.1, min.reads = 40, min.reads.normal = 10, min.reads.baf = 1,
                             max.mut.types = 1, min.type.freq = 0.9, min.fw.freq = 0, assembly = "hg38", gc.stats = NULL,
                             do_raster = FALSE, smooth_gc = FALSE, min_times_gc = 5, gc_grid = 250, parallel = 1,
                             weighted.mean = TRUE, segment_weights = list(
                               fit = 0.8, # Weight for actual fit quality
                               penalty = 0.4, # Weight for consecutive outliers penalty
                               elbow = 0.1, # Weight for proximity to elbow
                               segments = 0.1, # Weight for segment count before elbow
                               window = 0.05 # Window size penalty weight
                             ), ...) {
  # Track start time and memory
  start_time <- Sys.time()
  start_mem <- gc(reset = TRUE)
  start_mem_used <- sum(start_mem[, 2])

  # Initialize parameters with all arguments
  params <- initialize_extract_parameters(
    file = file, window = window, overlap = overlap,
    slide_win = slide_win, peak_wins = peak_wins, normalization.method = normalization.method,
    ignore.normal = ignore.normal, verbose = verbose, chromosome.list = chromosome.list,
    breaks = breaks, assembly = assembly, gc.stats = gc.stats, do_raster = do_raster,
    smooth_gc = smooth_gc, min_times_gc = min_times_gc, gc_grid = gc_grid, parallel = parallel,
    weighted.mean = weighted.mean, segment_weights = segment_weights, ...
  )

  # Validate input parameters
  validate_params(params)

  tryCatch({
    # Process GC content
    gc_data <- process_gc_content(params$gc.stats, params$normalization.method)

    gc_splines <- list(normal = smooth.spline(data.frame(
      gc = as.numeric(names(gc_data$normal_vect)),
      depth = gc_data$normal_vect
    )), tumor = smooth.spline(data.frame(
      gc = as.numeric(names(gc_data$tumor_vect)),
      depth = gc_data$tumor_vect
    )))

    # Initialize containers
    containers <- initialize_extract_containers(params$chromosome.list)

    # Process each chromosome with improved parallel handling
    if (params$parallel > 1) {
      cl <- NULL
      tryCatch(
        {
          cl <- manage_parallel_cluster(params$parallel)
          if (is.null(cl)) stop("Failed to create cluster")
          on.exit(if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE))

          # Export necessary objects
          parallel::clusterExport(cl,
            c("file", "params", "gc_splines"),
            envir = environment()
          )

          results <- pbapply::pblapply(
            seq_along(params$chromosome.list),
            function(idx) {
              chr <- params$chromosome.list[idx]
              process_single_chromosome(
                chr, file, params$gc.stats,
                gc_splines, NULL, params
              )
            },
            cl = cl
          )

          # Merge results back maintaining original data structure
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
          }
        },
        error = function(e) {
          message("Error in parallel processing: ", e$message)
          stop(e)
        }
      )
    } else {
      for (chr in params$chromosome.list) {
        containers <- process_single_chromosome(
          chr, file, params$gc.stats,
          gc_splines, containers, params
        )
      }
    }

    # Finalize and return results
    final_results <- finalize_extract_results(
      containers, params, params$gc.stats,
      gc_data
    )
    gc() # Force garbage collection

    # Print performance summary if verbose
    if (params$verbose) {
      end_time <- Sys.time()
      end_mem <- gc(reset = FALSE)
      end_mem_used <- sum(end_mem[, 2])

      # Track memory usage of all processes
      if (params$parallel > 1) {
        # Get all child process IDs (on Unix-like systems)
        child_pids <- tryCatch(
          {
            suppressWarnings(
              system(sprintf("pgrep -P %d", Sys.getpid()), intern = TRUE)
            )
          },
          error = function(e) character(0)
        )

        total_mem <- end_mem_used # Start with main process memory

        if (length(child_pids) > 0) {
          # Use ps command to get memory usage for each child process
          mem_cmd <- sprintf("ps -o rss= %s", paste(child_pids, collapse = " "))
          child_mems <- try(as.numeric(system(mem_cmd, intern = TRUE)) / 1024, silent = TRUE)

          if (!inherits(child_mems, "try-error")) {
            total_mem <- total_mem + sum(child_mems, na.rm = TRUE)
          }
        }
      } else {
        total_mem <- end_mem_used
      }

      # Calculate total mutations and covered bases
      total_mutations <- sum(sapply(final_results$mutations, nrow))
      total_bases <- sum(sapply(final_results$segments, function(segs) {
        sum(segs$end.pos - segs$start.pos + 1)
      }))
      total_mb <- total_bases / 1e6

      message("\nPerformance Summary:")
      message(sprintf(
        "Total time: %.2f minutes",
        as.numeric(difftime(end_time, start_time, units = "mins"))
      ))
      if (params$parallel > 1) {
        message(sprintf(
          "Peak memory usage (main process): %.2f GB",
          max(0, (end_mem_used) / 1024)
        ))
        message(sprintf(
          "Peak memory usage (all processes): %.2f GB",
          max(0, total_mem / 1024)
        ))
        message(sprintf(
          "Number of worker processes: %d",
          length(child_pids)
        ))
      } else {
        message(sprintf(
          "Peak memory usage: %.2f GB",
          max(0, total_mem / 1024)
        ))
      }
      message(sprintf(
        "Number of chromosomes processed: %d",
        length(params$chromosome.list)
      ))
      message(sprintf(
        "Total segments identified: %d",
        sum(sapply(final_results$segments, nrow))
      ))
      message(sprintf("Total mutations detected: %d", total_mutations))
      message(sprintf("Total megabases analyzed: %.1f", total_mb))
      message(sprintf(
        "Mutation rate: %.2f mutations/Mb",
        total_mutations / total_mb
      ))
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

# Extraction helper functions
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

  tryCatch(
    {
      if (normalization.method == "mean") {
        list(
          normal_vect = mean_gc(gc.stats$normal), tumor_vect = mean_gc(gc.stats$tumor),
          tum_depth = weighted.mean(x = gc.stats$tumor$depth, w = colSums(gc.stats$tumor$n)),
          nor_depth = weighted.mean(x = gc.stats$normal$depth, w = colSums(gc.stats$normal$n))
        )
      } else {
        list(
          normal_vect = median_gc(gc.stats$normal), tumor_vect = median_gc(gc.stats$tumor),
          tum_depth = weighted.median(x = gc.stats$tumor$depth, w = colSums(gc.stats$tumor$n)),
          nor_depth = weighted.median(x = gc.stats$normal$depth, w = colSums(gc.stats$normal$n))
        )
      }
    },
    error = function(e) {
      stop("GC content processing failed: ", e$message)
    }
  )
}

process_depths <- function(seqz.data, gc_splines, avg_depths, ignore.normal) {
  if (nrow(seqz.data) == 0) {
    return(list(
      tumor = numeric(0), normal = numeric(0), ratio = numeric(0),
      norm_gc_stats = NULL, seqz.data = seqz.data
    ))
  }

  # Calculate normalized depths
  tumor_depth <- seqz.data$depth.tumor / predict(gc_splines$tumor, seqz.data$GC.percent)$y
  normal_depth <- seqz.data$depth.normal / predict(gc_splines$normal, seqz.data$GC.percent)$y

  norm_gc_stats <- depths_gc(depth_n = round(
    normal_depth * avg_depths$normal,
    0
  ), depth_t = round(tumor_depth * avg_depths$tumor, 0), gc = seqz.data$GC.percent)

  # Calculate ratio and add it to seqz.data
  ratio <- if (ignore.normal) {
    round(tumor_depth, 3)
  } else {
    round(tumor_depth / normal_depth, 3)
  }

  # Add calculated values to seqz.data
  seqz.data$adjusted.ratio <- ratio
  seqz.data$depth.ratio <- seqz.data$depth.tumor / seqz.data$depth.normal

  list(
    tumor = tumor_depth, normal = normal_depth, ratio = ratio, norm_gc_stats = norm_gc_stats,
    seqz.data = seqz.data # Return modified seqz.data
  )
}

# Improved window calculation with memory optimization
calculate_windows <- function(seqz.data, depths, window, overlap, avg_depths) {
  if (nrow(seqz.data) == 0) {
    return(list(
      ratio = data.frame(), normal = data.frame(), tumor = data.frame(),
      raw_ratio = data.frame(), n_normal = data.frame(), n_tumor = data.frame(),
      baf = list()
    ))
  }

  required_cols <- c(
    "adjusted.ratio", "depth.ratio", "position", "chromosome",
    "depth.normal"
  )
  if (!all(required_cols %in% names(seqz.data))) {
    stop("Missing required columns in seqz.data: ", paste(setdiff(
      required_cols,
      names(seqz.data)
    ), collapse = ", "))
  }

  # Calculate all window values at once
  list(
    ratio = windowValues(
      x = seqz.data$adjusted.ratio, positions = seqz.data$position,
      chromosomes = seqz.data$chromosome, window = window, overlap = overlap, weight = seqz.data$depth.normal
    ),
    normal = windowValues(
      x = seqz.data$depth.normal / avg_depths$normal, positions = seqz.data$position,
      chromosomes = seqz.data$chromosome, window = window, overlap = overlap
    ),
    tumor = windowValues(
      x = seqz.data$depth.tumor / avg_depths$tumor, positions = seqz.data$position,
      chromosomes = seqz.data$chromosome, window = window, overlap = overlap
    ),
    raw_ratio = windowValues(
      x = seqz.data$depth.ratio, positions = seqz.data$position,
      chromosomes = seqz.data$chromosome, window = window, overlap = overlap,
      weight = seqz.data$depth.normal
    ), n_normal = windowValues(
      x = depths$normal,
      positions = seqz.data$position, chromosomes = seqz.data$chromosome, window = window,
      overlap = overlap
    ), n_tumor = windowValues(
      x = depths$tumor, positions = seqz.data$position,
      chromosomes = seqz.data$chromosome, window = window, overlap = overlap
    ),
    baf = list()
  )
}

initialize_extract_containers <- function(chromosome.list) {
  n_chr <- length(chromosome.list)
  containers <- list(
    windows.baf = vector("list", n_chr), windows.ratio = vector(
      "list",
      n_chr
    ), windows.raw_ratio = vector("list", n_chr), windows.normal = vector(
      "list",
      n_chr
    ), windows.tumor = vector("list", n_chr), windows.n_normal = vector(
      "list",
      n_chr
    ), windows.n_tumor = vector("list", n_chr), mutation.list = vector(
      "list",
      n_chr
    ), segments.list = vector("list", n_chr), norm.gc.list = vector(
      "list",
      n_chr
    ), rank_peaks.list = vector("list", n_chr) # Add new field
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
  containers$rank_peaks.list[[idx]] <- list(
    selected_win = results$segments$selected_win,
    peak_win = results$segments$peak_win
  )
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
    message("\nProcessing chromosome ", chr)
  }

  # Read chromosome data
  file.lines <- gc_stats$file.metrics[which(params$chr.vect == chr), ]
  seqz.data <- read.seqz(file, n_lines = c(file.lines$start, file.lines$end), chr_name = chr)

  # Process depths and get modified seqz.data
  depths_result <- process_depths(seqz.data, gc_splines, params$avg_depths, params$ignore.normal)
  seqz.data <- depths_result$seqz.data # Use updated seqz.data

  # Calculate windows with initialized data
  windows <- calculate_windows(
    seqz.data = seqz.data, depths = depths_result, window = params$window,
    overlap = params$overlap, avg_depths = params$avg_depths
  )

  # Process BAF if heterozygous positions exist
  seqz.het <- seqz.data[seqz.data$zygosity.normal == "het", ]
  num_het_positions <- nrow(seqz.het)
  if (num_het_positions > 0) {
    windows$baf <- windowBf(
      Af = seqz.het$Af, Bf = seqz.het$Bf, good.reads = seqz.het$good.reads,
      chromosomes = seqz.het$chromosome, positions = seqz.het$position, conf = 0.95,
      window = params$window, overlap = params$overlap
    )
  } else {
    windows$baf <- list(data.frame(
      start = min(seqz.data$position, na.rm = TRUE),
      end = max(seqz.data$position, na.rm = TRUE), mean = 0, q0 = 0, q1 = 0,
      N = 1
    ))
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

  mutations <- tryCatch(
    {
      mutation.table(seqz.data,
        mufreq.threshold = params$min.mut.freq, min.reads = params$min.reads,
        min.reads.normal = params$min.reads.normal, max.mut.types = params$max.mut.types,
        min.type.freq = params$min.type.freq, min.fw.freq = params$min.fw.freq,
        segments = segments$seg
      )
    },
    error = function(e) {
      message("Warning: Mutation table calculation failed: ", e$message)
      data.frame() # Return an empty data frame on error
    }
  )


  # Store results
  containers <- store_chromosome_results(
    list(
      windows = windows, segments = segments,
      mutations = mutations, norm_gc_stats = depths_result$norm_gc_stats
    ), containers,
    chr, which(params$chromosome.list == chr)
  )

  if (params$verbose) {
    log_chromosome_results(segments, seqz.data, mutations, num_het_positions)
  }

  containers
}

initialize_extract_parameters <- function(file, window, overlap = 1, slide_win = 100,
                                          peak_wins = seq(from = 50, to = 1000, by = 75), normalization.method = "mean",
                                          ignore.normal = FALSE, verbose = TRUE, chromosome.list = NULL, breaks = NULL,
                                          min.mut.freq = 0.1, min.reads = 40, min.reads.normal = 10, min.reads.baf = 1,
                                          max.mut.types = 1, min.type.freq = 0.9, min.fw.freq = 0, assembly = "hg38", gc.stats = NULL,
                                          do_raster = FALSE, smooth_gc = FALSE, min_times_gc = 5, gc_grid = 250, parallel = 1,
                                          weighted.mean = TRUE, segment_weights = list(
                                            fit = 0.8,
                                            penalty = 0.4, # Weight for consecutive outliers penalty
                                            elbow = 0.1,
                                            segments = 0.1,
                                            window = 0.05
                                          ), ...) {
  # Initialize GC stats if needed
  local_gc_stats <- if (is.null(gc.stats)) {
    gc.sample.stats(file,
      verbose = verbose, parallel = parallel, smooth = smooth_gc,
      min_times = min_times_gc, cl = parallel, n = gc_grid
    )
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
  list(
    window = window, overlap = overlap, slide_win = slide_win, peak_wins = peak_wins,
    normalization.method = normalization.method, ignore.normal = ignore.normal,
    verbose = verbose, assembly = assembly, chromosome.list = chromosome.list,
    breaks = if (is.null(dim(breaks))) NULL else breaks, chr.vect = chr.vect,
    min.mut.freq = min.mut.freq, min.reads = min.reads, min.reads.normal = min.reads.normal,
    min.reads.baf = min.reads.baf, max.mut.types = max.mut.types, min.type.freq = min.type.freq,
    min.fw.freq = min.fw.freq, gc.stats = local_gc_stats, do_raster = do_raster,
    smooth_gc = smooth_gc, min_times_gc = min_times_gc, gc_grid = gc_grid, parallel = parallel,
    ratio_baf_raster = data.frame(dr = NULL), avg_depths = list(normal = weighted.mean(
      local_gc_stats$normal$depth,
      colSums(local_gc_stats$normal$n)
    ), tumor = weighted.mean(
      local_gc_stats$tumor$depth,
      colSums(local_gc_stats$tumor$n)
    )), weighted.mean = weighted.mean, segment_weights = segment_weights
  )
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
    cl = params$parallel, grid_size = params$gc_grid
  )

  # Calculate final depths
  avg_tum_ndepth <- weighted.mean(x = gc_norm$tumor$depth, w = colSums(gc_norm$tumor$n))
  avg_nor_ndepth <- weighted.mean(x = gc_norm$normal$depth, w = colSums(gc_norm$normal$n))

  avg_depth_ratio <- if (params$ignore.normal) {
    avg_tum_ndepth / gc_data$tum_depth
  } else {
    (avg_tum_ndepth / gc_data$tum_depth) / (avg_nor_ndepth / gc_data$nor_depth)
  }

  depths <- list(
    avg_depth_ratio = avg_depth_ratio,
    avg_tum_depth = gc_data$tum_depth,
    avg_nor_depth = gc_data$nor_depth
  )

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
    list(
      mutations = mutations,
      megabases = bases / 1e6,
      mutations_per_mb = mutations / (bases / 1e6)
    )
  })
  names(chr_mutation_stats) <- params$chromosome.list

  # Compile complete mutation statistics
  mutation_stats <- list(
    total_mutations = total_mutations,
    total_megabases = total_bases / 1e6,
    mutations_per_mb = total_mutations / (total_bases / 1e6),
    chromosomes = chr_mutation_stats
  )

  # Construct and return complete results
  list(
    BAF = containers$windows.baf,
    ratio = containers$windows.ratio,
    raw_ratio = containers$windows.raw_ratio,
    depths = list(
      raw = list(
        normal = containers$windows.normal,
        tumor = containers$windows.tumor
      ),
      norm = list(
        normal = containers$windows.n_normal,
        tumor = containers$windows.n_tumor
      )
    ),
    mutations = containers$mutation.list,
    segments = containers$segments.list,
    win_peaks = containers$rank_peaks.list,
    chromosomes = params$chromosome.list,
    gc = gc_stats,
    gc_norm = gc_norm,
    avg.depth.ratio = depths$avg_depth_ratio,
    avg.depth.tumor = depths$avg_tum_depth,
    avg.depth.normal = depths$avg_nor_depth,
    mutation_stats = mutation_stats
  )
}

# Add safer parallel processing management function
manage_parallel_cluster <- function(n_cores) {
  if (n_cores > 1) {
    cl <- NULL
    tryCatch(
      {
        cl <- parallel::makeCluster(n_cores)
        # Load required packages on worker nodes
        parallel::clusterEvalQ(cl, {
          library(pbapply)
          library(stringr)
        })
        return(cl)
      },
      error = function(e) {
        if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE)
        stop("Failed to create cluster: ", e$message)
      }
    )
  }
  return(NULL)
}
