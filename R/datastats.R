splash_table <- function(lis_obj) {
    lis_obj <- Reduce("c", lis_obj)
    split(lis_obj, names(lis_obj))
}

unfold_data <- function(
    x, f = do_get_gc, args_f = get_gc_defaults, stats = TRUE) {

    data_items <- colnames(x)[3:ncol(x)]
    xyz_mats <- lapply(data_items, FUN = function(x, stat_obj,
        f, args) {
        f(stat_obj[, x], args)
    }, stat_obj = x, f = f, args = args_f)

    names(xyz_mats) <- data_items

    if (stats) {
        ord_chrom <- unique(Reduce("c", Reduce("c", x[, "unique"])))
        stats_chrom <- Reduce("c", x[, "lines"])
        stats_chrom <- sapply(splash_table(x[, "lines"]), sum)
        stats_chrom <- stats_chrom[ord_chrom]
        stats_start <- cumsum(c(1, stats_chrom[-length(stats_chrom)]))
        stats_end   <- stats_start + stats_chrom - 1
        stats_chrom <- data.frame(chr = ord_chrom, n_lines = stats_chrom,
            start = stats_start, end = stats_end)

        c(list(file.metrics = stats_chrom), xyz_mats)
    } else {
        xyz_mats
    }
}

data_fast_stats <- function(file, col_types = "c--dd----d----",
    buffer = 33554432, parallel = 2L, stats = TRUE, verbose = TRUE,
    col_sets = list("normal" = c(2, 4), "tumor" = c(3, 4)),
    f1 = gc_table, f2 = do_get_gc, args_f1 = list(), args_f2 = get_gc_defaults,
    msg = NULL) {
    con <- gzfile(file, "rb")

    suppressWarnings(skip_line <- readLines(con, n = 1))
    remove(skip_line)
    parse_chunck <- function(x, col_types,
        sets = col_sets, f = f1, args_f = args_f1) {
        x <- read_tsv(file = paste(mstrsplit(x), collapse = "\n"),
            col_types = col_types, col_names = FALSE,
            skip = 0, n_max = Inf, progress = FALSE)
        u_chr <- unique(x[, 1])
        n_chr <- table(x[, 1])
        set_lists <- lapply(sets, f, y = x, args = args_f)
        if (verbose) {
            message(".", appendLF = FALSE)
        }
        c(list(unique = u_chr, lines = n_chr), set_lists)
    }
    if (verbose) {
        message(msg, appendLF = FALSE)
    }
    res <- chunk.apply(
        input = con, FUN = function(x, col_types, sets, f, args) {
            parse_chunck(
                x, col_types = col_types, sets = sets, f = f, args)
        }, col_types = col_types, sets = col_sets, f = f1, args = args_f1,
        CH.MAX.SIZE = buffer, parallel = parallel)
    close(con)
    if (verbose) {
        message(" done")
    }
    #res 
    unfold_data(x = res, f = f2, args_f = args_f2, stats = stats)
}

get_baf_ratio <- function(baf_col, smooth = TRUE,
    min_times = 20, grid_size = 100, scale.subset = 1.5, ...) {
    sort_char <- function(x) {
        as.character(sort(as.numeric(x)))
    }
    all_ratios <- sequenza:::splash_table(baf_col)
    all_ratios <- lapply(all_ratios, FUN = function(x) {
        sapply(sequenza:::splash_table(x), sum)
    })
    names_baf <- sort_char(names(all_ratios))
    all_ratios <- all_ratios[names_baf]
    names_ratios <- sort_char(unique(Reduce("c", lapply(all_ratios, names))))
    n <- do.call(rbind, lapply(all_ratios, FUN = function(x, names_ratios) {
            res <- x[names_ratios]
            names(res) <- names_ratios
            res
        },
        names_ratios = names_ratios))
    n[is.na(n)] <- 0
    if (smooth == TRUE) {
        part <- gc_data_smooth(list(
            gc = as.numeric(names_baf), depth = as.numeric(
                names_ratios), n = n),
            min_times = min_times, n = grid_size,
            scale.subset = scale.subset, ...)
        names(part) <- c("x", "y", "z")
        part
    } else {
        list(baf = as.numeric(names_baf), ratio = as.numeric(
            names_ratios), n = n)
    }
}

baf_ratio_raster <- function(file_name, gc_normal, gc_tumor, verbose = TRUE,
    min_times = 20, smooth = TRUE, grid_size = 100, scale.subset = 4,
    round_baf = 2, round_ratio = 1) {

    ratio_args <- list(
        "model_normal" = gc_normal,
        "model_tumor" = gc_tumor,
        "round_dr" = round_ratio,
        "round_bf" = round_baf
    )

    get_gc_defaults <- list("smooth" = smooth,
        "min_times" = min_times, "grid_size" = grid_size,
        "scale_subset" = scale.subset)

    do_ratio <- function(x, y, args) {
        model_normal <- args[["model_normal"]]
        model_tumor <- args[["model_tumor"]]
        round_r <- args[["round_dr"]]
        round_b <- args[["round_bf"]]
        colnames(y)[x] <- names(x)
        y <- y[y$zyg == "het", ]
        norm_tumor_depth <- y$tumor /
            predict(model_tumor, y$gc)$y
        norm_normal_depth <- y$normal /
            predict(model_normal, y$gc)$y
        d_r <- round(norm_tumor_depth / norm_normal_depth, round_r)
        baf <- round(y$baf, round_b)

        lapply(split(c(d_r, d_r), c(baf, 1 - baf)), table)
    }


    do_get_ratio_baf <- function(x, args) {
        smooth <- args[["smooth"]]
        min_times <- args[["min_times"]]
        grid_size <- args[["grid_size"]]
        scale_subset <- args[["scale_subset"]]
        get_baf_ratio(x, smooth = smooth, min_times = min_times,
                grid_size = grid_size, scale.subset = scale_subset) 
    }

    data_fast_stats(file_name, col_types = "---dd--dcd----",
        stats = FALSE, verbose = verbose, col_sets = list("dr" = c(
            "normal" = 1, "tumor" = 2, "baf" = 3, "zyg" = 4, "gc" = 5)),
        f1 = do_ratio, f2 = do_get_ratio_baf, args_f1 = ratio_args,
        args_f2 = get_gc_defaults,
        msg = "Collecting BAF/ratio information ")
}

rs_baf_ratio <- function(dens, n = 100, min_prob = 1) {
    df_r <- approxfun(density(rep(dens$y, times = colSums(dens$z))))
    df_b <- approxfun(density(rep(dens$x, times = rowSums(dens$z))))
    #df_r <- splinefun(x = dens$y, y = colSums(dens$z))
    #df_b <- splinefun(x = dens$x, y = rowSums(dens$z))
    bs <- 0:1000 / 2000
    rs <- 0:2000 / 100
    p_bs <- df_b(bs)
    p_rs <- df_r(rs)
    p_bs[is.na(p_bs)] <- 0
    p_rs[is.na(p_rs)] <- 0
    p_bs[p_bs <= 0] <- 0
    p_rs[p_rs <= 0] <- 0

    sample_b <- sample(bs, n, TRUE, p_bs)
    sample_r <- sample(rs, n, TRUE, p_rs)
    res <- cbind(
        baf = sample_b,
        ratio = sample_r,
        l = df_b(
            sample_b) / max(p_bs) * df_r(sample_r) / max(p_rs))
    res[res[, 3] >= min_prob, 1:2]
}
