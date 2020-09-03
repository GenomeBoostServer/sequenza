# file = data.file; window = 1e6; overlap = 1;
#     slide_win = 100; peak_wins = seq(from = 50, to = 300, by = 25);
#     mufreq.treshold = 0.10; min.reads = 40; min.reads.normal = 10;
#     min.reads.baf = 1; max.mut.types = 1; min.type.freq = 0.9;
#     min.fw.freq = 0; verbose = TRUE; chromosome.list = NULL;
#     breaks = NULL; assembly = "hg19"; weighted.mean = TRUE;
#     normalization.method = "mean"; ignore.normal = FALSE;
#     parallel = 1; gc.stats = NULL; segments.samples = FALSE

sequenza.extract <- function(file, window = 1e6, overlap = 1,
    slide_win = 100, peak_wins = seq(from = 0.0005, to = 0.05, by = 0.005),
    mufreq.treshold = 0.10, min.reads = 40, min.reads.normal = 10,
    min.reads.baf = 1, max.mut.types = 1, min.type.freq = 0.9,
    min.fw.freq = 0, verbose = TRUE, chromosome.list = NULL,
    breaks = NULL, assembly = "hg19", weighted.mean = TRUE,
    normalization.method = "mean", ignore.normal = FALSE,
    parallel = 1, gc.stats = NULL, segments.samples = FALSE,
    smooth_gc = TRUE, min_times_gc = 20, gc_grid = 250) {

    pbo <- pboptions()

    if (verbose == FALSE) {
        pboptions(type = "none")
    }

    if (is.null(gc.stats)) {
        gc.stats <- gc.sample.stats(file, verbose = verbose,
            parallel = parallel, smooth = smooth_gc, min_times = min_times_gc,
                cl = parallel, n = gc_grid)
    }
    if (normalization.method == "mean") {
        gc.normal.vect <- mean_gc(gc.stats$normal)
        gc.tumor.vect  <- mean_gc(gc.stats$tumor)
        avg_tum_depth <- weighted.mean(x = gc.stats$tumor$depth,
            w = colSums(gc.stats$tumor$n))
        avg_nor_depth <- weighted.mean(x = gc.stats$normal$depth,
            w = colSums(gc.stats$normal$n))
    } else {
        gc.normal.vect <- median_gc(gc.stats$normal)
        gc.tumor.vect  <- median_gc(gc.stats$tumor)
        avg_tum_depth <- weighted.median(x = gc.stats$tumor$depth,
            w = colSums(gc.stats$tumor$n))
        avg_nor_depth <- weighted.median(x = gc.stats$normal$depth,
            w = colSums(gc.stats$normal$n))
    }

    gc_spline_normal <- smooth.spline(
        data.frame(
            gc = as.numeric(names(gc.normal.vect)),
            depth = gc.normal.vect))
    gc_spline_tumor <- smooth.spline(
        data.frame(
            gc = as.numeric(names(gc.tumor.vect)),
            depth = gc.tumor.vect))


    windows.baf   <- list()
    windows.ratio <- list()
    windows.raw_ratio <- list()
    windows.normal <- list()
    windows.tumor <- list()
    windows.n_normal <- list()
    windows.n_tumor <- list()
    rank_peaks.list <- list()
    mutation.list <- list()
    segments.list <- list()
    segments_samples.list <- list()
    norm.gc.list <- list()
    if (is.null(dim(breaks))) {
        breaks <- NULL
    }
    chr.vect <- as.character(gc.stats$file.metrics$chr)
    if (is.null(chromosome.list)) {
        chromosome.list <- chr.vect
    } else {
        chromosome.list <- chromosome.list[chromosome.list %in% chr.vect]
    }
    for (chr in chromosome.list) {
        if (verbose) {
            message("Processing ", chr, ":", appendLF = TRUE)
        }
        tbi <- file.exists(paste0(file, ".tbi"))
        if (tbi) {
            seqz.data   <- read.seqz(file, chr_name = chr)
        } else {
            file.lines <- gc.stats$file.metrics[which(chr.vect == chr), ]
            seqz.data   <- read.seqz(file, n_lines = c(file.lines$start,
                file.lines$end))
        }

        norm_tumor_depth <- seqz.data$depth.tumor /
            predict(gc_spline_tumor, seqz.data$GC.percent)$y
        norm_normal_depth <- seqz.data$depth.normal /
            predict(gc_spline_normal, seqz.data$GC.percent)$y

        norm.gc.stats <- depths_gc(
            depth_n = round(norm_normal_depth * avg_nor_depth, 0),
            depth_t = round(norm_tumor_depth * avg_tum_depth, 0),
            gc = seqz.data$GC.percent)
        if (ignore.normal) {
            seqz.data$adjusted.ratio <- round(norm_tumor_depth, 3)
        } else {
            seqz.data$adjusted.ratio <- round(
                norm_tumor_depth / norm_normal_depth, 3)
        }
        if (segments.samples == TRUE) {
            breaks_normal_chr <- NULL
            breaks_tumor_chr <- NULL
            # breaks_normal_chr <- breaks_full(
            #     data = data.frame(chromosome = seqz.data$chromosome,
            #                       position = seqz.data$position,
            #                       adjusted.ratio = norm_normal_depth,
            #                       singsAsFactors = FALSE),
            #     gamma = gamma.pcf, kmin = kmin.pcf, assembly = assembly,
            #     breaks.het = NULL)
            # breaks_tumor_chr <- breaks_full(
            #    data = data.frame(chromosome = seqz.data$chromosome,
            #                      position = seqz.data$position,
            #                      adjusted.ratio = norm_tumor_depth,
            #                      singsAsFactors = FALSE),
            #    gamma = gamma.pcf, kmin = kmin.pcf, assembly = assembly,
            #    breaks.het = NULL)
        } else {
           breaks_normal_chr <- NULL
           breaks_tumor_chr <- NULL
        }
        seqz.r.win <- windowValues(x = seqz.data$adjusted.ratio,
            positions = seqz.data$position,
            chromosomes = seqz.data$chromosome,
            window = window, overlap = overlap,
            weight = seqz.data$depth.normal)
        seqz.n.win <- windowValues(x = seqz.data$depth.normal / avg_nor_depth,
            positions = seqz.data$position,
            chromosomes = seqz.data$chromosome,
            window = window, overlap = overlap)
        seqz.t.win <- windowValues(x = seqz.data$depth.tumor / avg_tum_depth,
            positions = seqz.data$position,
            chromosomes = seqz.data$chromosome,
            window = window, overlap = overlap)
        seqz.r_r.win <- windowValues(x = seqz.data$depth.ratio / (
                avg_tum_depth / avg_nor_depth),
            positions = seqz.data$position,
            chromosomes = seqz.data$chromosome,
            window = window, overlap = overlap,
            weight = seqz.data$depth.normal)
        seqz.n_n.win <- windowValues(x = norm_normal_depth,
            positions = seqz.data$position,
            chromosomes = seqz.data$chromosome,
            window = window, overlap = overlap)
        seqz.n_t.win <- windowValues(x = norm_tumor_depth,
            positions = seqz.data$position,
            chromosomes = seqz.data$chromosome,
            window = window, overlap = overlap)


        seqz.hom <- seqz.data$zygosity.normal == "hom"
        seqz.het <- seqz.data[!seqz.hom, ]
        het.filt <- seqz.het$good.reads >= min.reads.baf
        seqz.het <- seqz.het[het.filt, ]
        het_ok <- nrow(seqz.het) > 0
        if (is.null(breaks)) {
            breaks_chr <- NULL
        } else {
            breaks_chr <- breaks[breaks$chrom == chr, ]
        }
        if (het_ok) {
            seqz.b.win <- windowBf(Af = seqz.het$Af, Bf = seqz.het$Bf,
                good.reads = seqz.het$good.reads,
                chromosomes = seqz.het$chromosome,
                positions = seqz.het$position, conf = 0.95,
                window = window, overlap = overlap)
        } else {
            seqz.b.win <- list()
            seqz.b.win[[1]] <- data.frame(start = min(seqz.data$position,
                na.rm = TRUE), end = max(seqz.data$position, na.rm = TRUE),
                mean = 0, q0 = 0,  q1 = 0, N = 1)
        }
        #diff_track <- slide_tracks(seqz.het, slide_win,
        #    signal_out = "both", verbose = verbose)
        breaks_chr_list <- lapply(peak_wins, FUN = function(
            x, data, data_het, breaks, slide_win,
            assembly, chromosome, verbose,
            min.reads.baf, weighted.mean) {
            # breaks_chr <- extract_breaks_tracks(
            #     track = diff_track, breaks = breaks,
            breaks_chr <- extract_breaks_tracks(
                track = data_het, breaks = breaks,
                peak_win = x, assembly = assembly,
                chromosome = chr)
            if (class(breaks_chr) == "try-error") {
                breaks_chr <- NULL
            }
            if (is.null(breaks_chr) || nrow(breaks_chr) == 0 ||
                length(breaks_chr) == 0) {
                breaks_chr <- data.frame(chrom = chr,
                    start.pos = min(seqz.data$position, na.rm = TRUE),
                    end.pos = max(seqz.data$position, na.rm = TRUE))
            }
            segment.breaks(seqz.tab = data, breaks = breaks_chr,
                min.reads.baf = min.reads.baf,
                weighted.mean = weighted.mean)

        }, data = seqz.data, data_het = seqz.het,
            breaks = breaks_chr, slide_win = slide_win,
            assembly = assembly, chromosome = chr, verbose = verbose,
            min.reads.baf = min.reads.baf, weighted.mean = weighted.mean)

        names(breaks_chr_list) <- as.character(peak_wins)

        compare_bins_list <- lapply(breaks_chr_list, FUN = function(
            x, baf_win, ratio_win) {
            baf_vs_bins <- compare_bins(
                x$start.pos, x$end.pos, x$Bf, baf_win)
            ratio_vs_bins <- compare_bins(
                x$start.pos, x$end.pos, x$depth.ratio, ratio_win)
            cbind(baf_fit = baf_vs_bins, ratio_fit = ratio_vs_bins)
        }, baf_win = seqz.b.win[[chr]], ratio_win = seqz.r.win[[chr]])

        compare_bins_segs <- data.frame(peak_win = peak_wins,
            do.call(rbind, compare_bins_list))
        compare_bins_segs$n_segs <- sapply(breaks_chr_list, nrow)

        ranks_fits <- cbind(
            apply(-compare_bins_segs[,
                c("baf_fit", "ratio_fit")], 2, rank, ties.method = "max"),
            n_segs = rank(compare_bins_segs$n_segs, ties.method = "min"))

        best_fits <- which(rowSums(ranks_fits) %in%  min(rowSums(ranks_fits)))
        select_win <- max(compare_bins_segs$peak_win[best_fits])

        seg.s1 <- breaks_chr_list[[as.character(select_win)]]

        mut.tab <- mutation.table(
            seqz.data, mufreq.treshold = mufreq.treshold,
            min.reads = min.reads, min.reads.normal = min.reads.normal,
            max.mut.types = max.mut.types, min.type.freq = min.type.freq,
            min.fw.freq = min.fw.freq, segments = seg.s1)

        windows.ratio[[which(chromosome.list == chr)]] <- seqz.r.win[[1]]
        windows.raw_ratio[[which(chromosome.list == chr)]] <- seqz.r_r.win[[1]]
        windows.normal[[which(chromosome.list == chr)]] <- seqz.n.win[[1]]
        windows.tumor[[which(chromosome.list == chr)]] <- seqz.t.win[[1]]
        windows.n_normal[[which(chromosome.list == chr)]] <- seqz.n_n.win[[1]]
        windows.n_tumor[[which(chromosome.list == chr)]] <- seqz.n_t.win[[1]]
        windows.baf[[which(chromosome.list == chr)]]   <- seqz.b.win[[1]]
        segments.list[[which(chromosome.list == chr)]] <- seg.s1
        mutation.list[[which(chromosome.list == chr)]] <- mut.tab
        norm.gc.list[[which(chromosome.list == chr)]] <- norm.gc.stats
        segments_samples.list[[which(chromosome.list == chr)]] <- list(
            normal = breaks_normal_chr, tumor = breaks_tumor_chr)
        rank_peaks.list[[which(chromosome.list == chr)]] <- list(
            selected_win = select_win, peak_win = compare_bins_segs)

        if (verbose) {
            message("   ", nrow(mut.tab), " variant calls.", appendLF = TRUE)
            message("   ", nrow(seg.s1), " copy-number segments.",
                appendLF = TRUE)
            message("   ", nrow(seqz.het), " heterozygous positions.",
                appendLF = TRUE)
            message("   ", sum(seqz.hom), " homozygous positions.",
                appendLF = TRUE)
        }
    }
    names(windows.baf)   <- chromosome.list
    names(windows.ratio) <- chromosome.list
    names(windows.raw_ratio) <- chromosome.list
    names(windows.normal) <- chromosome.list
    names(windows.tumor) <- chromosome.list
    names(windows.n_normal) <- chromosome.list
    names(windows.n_tumor) <- chromosome.list
    names(mutation.list) <- chromosome.list
    names(segments.list) <- chromosome.list
    names(segments_samples.list) <- chromosome.list
    names(rank_peaks.list) <- chromosome.list

    gc_norm <- unfold_gc(do.call(rbind, norm.gc.list), stats = FALSE,
        smooth = smooth_gc, min_times = min_times_gc,
        cl = parallel, grid_size = gc_grid)

    if (normalization.method == "mean") {
        avg_tum_ndepth <- weighted.mean(x = gc_norm$tumor$depth,
            w = colSums(gc_norm$tumor$n))
        avg_nor_ndepth <- weighted.mean(x = gc_norm$normal$depth,
            w = colSums(gc_norm$normal$n))
    } else {
        avg_tum_ndepth <- weighted.median(x = gc_norm$tumor$depth,
            w = colSums(gc_norm$tumor$n))
        avg_nor_ndepth <- weighted.median(x = gc_norm$normal$depth,
            w = colSums(gc_norm$normal$n))
    }
    if (ignore.normal) {
        avg_depth_ratio <- avg_tum_ndepth / avg_tum_depth
    } else {
        avg_depth_ratio <- (avg_tum_ndepth / avg_tum_depth) /
            (avg_nor_ndepth / avg_nor_depth)
    }
    pboptions(pbo)

    list(BAF = windows.baf, ratio = windows.ratio,
        raw_ratio = windows.raw_ratio,
        depths =  list(
            raw = list(normal = windows.normal, tumor = windows.tumor),
            norm = list(normal = windows.n_normal, tumor = windows.n_tumor)),
        mutations = mutation.list, segments = segments.list,
        chromosomes = chromosome.list, gc = gc.stats,
        gc_norm = gc_norm, avg.depth.ratio = avg_depth_ratio,
        avg.depth.tumor = avg_tum_depth, avg.depth.normal = avg_nor_depth,
        segments_samples = segments_samples.list, win_peaks = rank_peaks.list)
}
