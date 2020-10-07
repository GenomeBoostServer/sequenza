median_gc <- function(gc_list) {
    if (any(gc_list$n < 0)) {
        gc_list$n <- gc_list$n + abs(min(gc_list$n, na.rm = TRUE))
    }
    apply(gc_list$n, 1, FUN = function(x, w) {

            weighted.median(x = w, w = x, na.rm = TRUE)
        },
        w = gc_list$depth)
}

mean_gc <- function(gc_list) {
    if (any(gc_list$n < 0)) {
        gc_list$n <- gc_list$n + abs(min(gc_list$n, na.rm = TRUE))
    }
    apply(gc_list$n, 1, FUN = function(x, w) {
            weighted.mean(x = w, w = x, na.rm = TRUE)
        },
        w = gc_list$depth)
}

depths_gc <- function(depth_n, depth_t, gc) {
    gc_nor <- lapply(split(depth_n, gc), table)
    gc_tum <- lapply(split(depth_t, gc), table)
    list(normal = gc_nor, tumor = gc_tum)
}

gc_data_smooth <- function(gc_list, min_times = 20, n = 100,
    scale.subset = 1.5, ...) {

    mengc <- mean_gc(gc_list)
    medgc <- median_gc(gc_list)
    max_depth <- round(max(c(mengc, medgc)) * scale.subset, 0)

    regrid <- smooth_matrix(x = gc_list$gc,
        y = gc_list$depth[gc_list$depth <= max_depth],
        z = gc_list$n[, gc_list$depth <= max_depth])
    n_tab <- regrid$estimate
    colnames(n_tab) <- as.character(regrid$eval.points[[2]])
    rownames(n_tab) <- as.character(regrid$eval.points[[1]])
    list(gc = regrid$eval.points[[1]],
        depth = regrid$eval.points[[2]], n = n_tab)
}

get_gc <- function(gc_col, smooth = TRUE,
    min_times = 20, grid_size = 100, scale.subset = 1.5, ...) {
    sort_char <- function(x) {
        as.character(sort(as.numeric(x)))
    }
    all_depths <- splash_table(gc_col)
    all_depths <- lapply(all_depths, FUN = function(x) {
        sapply(splash_table(x), sum)
    })
    names_gc <- sort_char(names(all_depths))
    all_depths <- all_depths[names_gc]
    names_depths <- sort_char(unique(Reduce("c", lapply(all_depths, names))))
    n <- do.call(rbind, lapply(all_depths, FUN = function(x, names_depths) {
            res <- x[names_depths]
            names(res) <- names_depths
            res
        },
        names_depths = names_depths))
    n[is.na(n)] <- 0
    if (smooth == TRUE) {
        gc_data_smooth(list(
            gc = as.numeric(names_gc), depth = as.numeric(names_depths), n = n),
            min_times = min_times, n = grid_size,
            scale.subset = scale.subset, ...)
    } else {
        list(gc = as.numeric(names_gc), depth = as.numeric(names_depths), n = n)
    }
}

gc_table <- function(x, y, args) {
    lapply(split(y[, x[1]],
        y[, x[2]]), table)
}

do_get_gc <- function(x, args) {
    smooth <- args[["smooth"]]
    min_times <- args[["min_times"]]
    grid_size <- args[["grid_size"]]
    scale_subset <- args[["scale_subset"]]
    get_gc(x, smooth = smooth, min_times = min_times,
        grid_size = grid_size, scale.subset = scale_subset)
}

gc.sample.stats <- function(file, col_types = "c--dd----d----",
    buffer = 33554432, parallel = 2L, stats = TRUE, smooth = TRUE,
    min_times = 20, n = 100, scale.subset = 1.5, verbose = TRUE, ...) {

    get_gc_args <- list("smooth" = smooth,
        "min_times" = min_times, "grid_size" = n,
        "scale_subset" = scale.subset)

    data_fast_stats(file, col_types = col_types,
        buffer = buffer, parallel = parallel, stats = stats,
        verbose = verbose,
        col_sets = list("normal" = c(2, 4), "tumor" = c(3, 4)),
        f1 = gc_table, f2 = do_get_gc,
        args_f1 = list(), args_f2 = get_gc_args,
        msg = "Collecting GC information ")
}

unfold_gc <- function(x, stats = FALSE, smooth = TRUE, min_times = 20,
        cl = 2L, grid_size = 100) {
    get_gc_args <- list("smooth" = smooth,
        "min_times" = min_times, "grid_size" = grid_size,
        "scale_subset" = 1.5)
    unfold_data(x = x, f = do_get_gc, args_f = get_gc_args, stats = stats)
}

unfold_gc_old <- function(x, stats = TRUE, smooth = TRUE,
    min_times = 20, grid_size = 100, scale.subset = 1.5, ...) {
    gc_norm <- get_gc(x[, "gc_nor"], smooth, min_times,
        grid_size, scale.subset, ...)
    gc_tum <- get_gc(x[, "gc_tum"], smooth, min_times,
        grid_size, scale.subset, ...)
    if (stats) {
        ord_chrom <- unique(Reduce("c", Reduce("c", x[, "unique"])))
        stats_chrom <- Reduce("c", x[, "lines"])
        stats_chrom <- sapply(splash_table(x[, "lines"]), sum)
        stats_chrom <- stats_chrom[ord_chrom]
        stats_start <- cumsum(c(1, stats_chrom[-length(stats_chrom)]))
        stats_end   <- stats_start + stats_chrom - 1
        stats_chrom <- data.frame(chr = ord_chrom, n_lines = stats_chrom,
            start = stats_start, end = stats_end)

        list(file.metrics = stats_chrom, normal = gc_norm, tumor = gc_tum)
    } else {
        list(normal = gc_norm, tumor = gc_tum)
    }
}