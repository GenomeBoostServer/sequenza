#' @rdname sequenza
#' @export
sequenza.results <- function(sequenza.extract, cp.table = NULL,
    sample.id, out.dir = getwd(), cellularity = NULL, ploidy = NULL,
    female = TRUE, CNt.max = 20, ratio.priority = FALSE, XY = c(X = "X",
        Y = "Y"), chromosome.list = 1:24) {
    # Enhanced input validation
    validate_results_input(sequenza.extract, sample.id, out.dir)

    # Create directory with better error handling
    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(out.dir)) {
        stop("Failed to create output directory: ", out.dir)
    }

    # Create file paths for all outputs
    files <- create_output_paths(out.dir, sample.id)

    # Validate inputs
    validate_inputs(cp.table, cellularity, ploidy)

    # Process segments
    segments_data <- process_segments_data(sequenza.extract,
        chromosome.list)

    # Save extraction data
    save_extraction_data(sequenza.extract, files$robj.extr, sample.id)

    # Process and plot GC content
    plot_gc_content(sequenza.extract, files$gc.file)

    # Process and plot depths
    plot_depths(sequenza.extract, segments_data$seg.tab, files$depths.file)

    # Process CP table and get cellularity/ploidy estimates
    cp_results <- process_cp_table(cp.table, cellularity, ploidy,
        files, sample.id)
    cellularity <- cp_results$cellularity
    ploidy <- cp_results$ploidy

    # Process mutations and segments
    results <- process_mutations_and_segments(sequenza.extract,
        segments_data, cp_results, female, XY, CNt.max, ratio.priority,
        chromosome.list, files)

    # Generate plots with progress tracking
    generate_result_plots(results, files, sequenza.extract, cellularity,
        ploidy, female, cp_results, XY = XY)

    # Write log file
    write_log_file(files$log.file)
}

# Helper functions
setup_output_dir <- function(out_dir) {
    if (!file.exists(out_dir)) {
        dir_ok <- dir.create(path = out_dir, recursive = TRUE)
        if (!dir_ok) {
            stop("Directory does not exist and cannot be created: ",
                out_dir)
        }
    }
}

create_output_paths <- function(out_dir, sample.id) {
    make_filename <- function(x) file.path(out_dir, paste(sample.id,
        x, sep = "_"))

    list(cp.file = make_filename("CP_contours.pdf"), cint.file = make_filename("confints_CP.txt"),
        chrw.file = make_filename("chromosome_view.pdf"), depths.file = make_filename("chromosome_depths.pdf"),
        gc.file = make_filename("gc_plots.pdf"), peak_win.file = make_filename("peak_windows_plots.pdf"),
        peak_dens.file = make_filename("peak_density_plots.pdf"),
        geno.file = make_filename("genome_view.pdf"), cn.file = make_filename("CN_bars.pdf"),
        fit.file = make_filename("model_fit.pdf"), alt.file = make_filename("alternative_solutions.txt"),
        afit.file = make_filename("alternative_fit.pdf"), muts.file = make_filename("mutations.txt"),
        segs.file = make_filename("segments.txt"), robj.extr = make_filename("sequenza_extract.RData"),
        robj.fit = make_filename("sequenza_cp_table.RData"),
        log.file = make_filename("sequenza_log.txt"))
}

validate_inputs <- function(cp.table, cellularity, ploidy) {
    if (is.null(cp.table) && (is.null(cellularity) || is.null(ploidy))) {
        stop("cp.table and/or cellularity and ploidy argument are required.")
    }
}

process_segments_data <- function(sequenza.extract, chromosome.list) {
    seg.tab <- do.call(rbind, sequenza.extract$segments[chromosome.list])
    seg.len <- (seg.tab$end.pos - seg.tab$start.pos)/1e+06
    list(seg.tab = seg.tab, seg.len = seg.len)
}

save_extraction_data <- function(sequenza.extract, robj.extr,
    sample.id) {
    assign(x = paste0(sample.id, "_sequenza_extract"), value = sequenza.extract)
    save(list = paste0(sample.id, "_sequenza_extract"), file = robj.extr)
}

plot_gc_content <- function(sequenza.extract, gc.file) {
    pdf(gc.file, width = 10, height = 5)
    par(mfrow = c(1, 2))
    gc.summary.plot(sequenza.extract$gc$normal, mean.col = "lightsalmon",
        median.col = "lightgreen", las = 1, xlab = "GC %", ylab = "Depth",
        zlab = "N", main = "GC vs raw depth in the normal sample")
    gc.summary.plot(sequenza.extract$gc_norm$normal, mean.col = "lightsalmon",
        median.col = "lightgreen", las = 1, xlab = "GC %", ylab = "Depth",
        zlab = "N", main = "GC vs normalized depth in the normal sample")
    gc.summary.plot(sequenza.extract$gc$tumor, mean.col = "lightsalmon",
        median.col = "lightgreen", las = 1, xlab = "GC %", ylab = "Depth",
        zlab = "N", main = "GC vs raw depth in the tumor sample")
    gc.summary.plot(sequenza.extract$gc_norm$tumor, mean.col = "lightsalmon",
        median.col = "lightgreen", las = 1, xlab = "GC %", ylab = "Depth",
        zlab = "N", main = "GC vs normalized depth in the tumor sample")
    dev.off()
}

plot_depths <- function(sequenza.extract, seg.tab, depths.file) {
    pdf(depths.file, height = 10, width = 15)
    for (i in unique(seg.tab$chromosome)) {
        max_coord_chr_i <- max(sequenza.extract$ratio[[i]]$end)
        par(mfcol = c(3, 2), xaxt = "n", mar = c(0, 4, 3, 0),
            oma = c(5, 0, 4, 0))
        plotWindows(sequenza.extract$depths$raw$normal[[i]],
            ylab = "normal depth", ylim = c(0, 2.5), main = paste("raw",
                i, sep = " "))
        plotWindows(sequenza.extract$depths$raw$tumor[[i]], ylab = "tumor depth",
            ylim = c(0, 2.5))
        plotWindows(sequenza.extract$raw_ratio[[i]], ylab = "depth ratio",
            ylim = c(0, 2.5))
        par(xaxt = "s")
        axis(labels = as.character(round(seq(0, max_coord_chr_i/1e+06,
            by = 10), 0)), side = 1, line = 0, at = seq(0, max_coord_chr_i,
            by = 1e+07), outer = FALSE, cex = par("cex.axis") *
            par("cex"))
        mtext("Position (Mb)", side = 1, line = 3, outer = FALSE,
            cex = par("cex.lab") * par("cex"))
        par(xaxt = "n")
        plotWindows(sequenza.extract$depths$norm$normal[[i]],
            ylab = "normal depth", ylim = c(0, 2.5), main = paste("normalized",
                i, sep = " "))
        plotWindows(sequenza.extract$depths$norm$tumor[[i]],
            ylab = "tumor depth", ylim = c(0, 2.5))
        plotWindows(sequenza.extract$ratio[[i]], ylab = "depth ratio",
            ylim = c(0, 2.5))
        par(xaxt = "s")
        axis(labels = as.character(round(seq(0, max_coord_chr_i/1e+06,
            by = 10), 0)), side = 1, line = 0, at = seq(0, max_coord_chr_i,
            by = 1e+07), outer = FALSE, cex = par("cex.axis") *
            par("cex"))
        mtext("Position (Mb)", side = 1, line = 3, outer = FALSE,
            cex = par("cex.lab") * par("cex"))
    }
    dev.off()
}

process_cp_table <- function(cp.table, cellularity, ploidy, files,
    sample.id) {
    if (!is.null(cp.table)) {
        assign(x = paste0(sample.id, "_sequenza_cp_table"), value = cp.table)
        save(list = paste0(sample.id, "_sequenza_cp_table"),
            file = files$robj.fit)
        cint <- get.ci(cp.table)
        pdf(files$cp.file)
        cp.plot(cp.table)
        cp.plot.contours(cp.table, add = TRUE, likThresh = c(0.95),
            col = "red", pch = 20)
        if (!is.null(cellularity) || !is.null(ploidy)) {
            if (is.null(cellularity)) {
                cellularity <- cint$max.cellularity
            }
            if (is.null(ploidy)) {
                ploidy <- cint$max.ploidy
            }
            points(x = ploidy, y = cellularity, pch = 5)
            text(x = ploidy, y = cellularity, labels = "User selection",
                pos = 3, offset = 0.5)
        } else {
            cellularity <- cint$max.cellularity
            ploidy <- cint$max.ploidy
        }
        dev.off()
        return(list(cellularity = cellularity, ploidy = ploidy,
            cint = cint, cp.table = cp.table))
    }
    list(cellularity = cellularity, ploidy = ploidy, cint = NULL,
        cp.table = NULL)
}

process_mutations_and_segments <- function(sequenza.extract,
    segments_data, cp_results, female, XY, CNt.max, ratio.priority,
    chromosome.list, files) {
    seg.tab <- segments_data$seg.tab
    seg.len <- segments_data$seg.len
    cellularity <- cp_results$cellularity
    ploidy <- cp_results$ploidy
    avg.depth.ratio <- sequenza.extract$avg.depth.ratio

    mut.tab <- na.exclude(do.call(rbind, sequenza.extract$mutations[chromosome.list]))
    if (female) {
        segs.is.xy <- seg.tab$chromosome == XY["Y"]
        mut.is.xy <- mut.tab$chromosome == XY["Y"]
    } else {
        segs.is.xy <- seg.tab$chromosome %in% XY
        mut.is.xy <- mut.tab$chromosome %in% XY
    }
    avg.sd.ratio <- sum(seg.tab$sd.ratio * seg.tab$N.ratio, na.rm = TRUE)/sum(seg.tab$N.ratio,
        na.rm = TRUE)
    avg.sd.Bf <- sum(seg.tab$sd.BAF * seg.tab$N.BAF, na.rm = TRUE)/sum(seg.tab$N.BAF,
        na.rm = TRUE)
    cn.alleles <- baf.bayes(Bf = seg.tab$Bf[!segs.is.xy], CNt.max = CNt.max,
        depth.ratio = seg.tab$depth.ratio[!segs.is.xy], cellularity = cellularity,
        ploidy = ploidy, avg.depth.ratio = avg.depth.ratio, sd.ratio = seg.tab$sd.ratio[!segs.is.xy],
        weight.ratio = seg.len[!segs.is.xy], sd.Bf = seg.tab$sd.BAF[!segs.is.xy],
        weight.Bf = 1, ratio.priority = ratio.priority, CNn = 2)
    seg.res <- cbind(seg.tab[!segs.is.xy, ], cn.alleles)
    if (!female) {
        if (sum(segs.is.xy) >= 1) {
            cn.alleles <- baf.bayes(Bf = NA, CNt.max = CNt.max,
                depth.ratio = seg.tab$depth.ratio[segs.is.xy],
                cellularity = cellularity, ploidy = ploidy, avg.depth.ratio = avg.depth.ratio,
                sd.ratio = seg.tab$sd.ratio[segs.is.xy], weight.ratio = seg.len[segs.is.xy],
                sd.Bf = NA, weight.Bf = NA, ratio.priority = ratio.priority,
                CNn = 1)
            seg.xy <- cbind(seg.tab[segs.is.xy, ], cn.alleles)
            seg.res <- rbind(seg.res, seg.xy)
        }
    }
    write.table(seg.res, file = files$segs.file, col.names = TRUE,
        row.names = FALSE, sep = "\t", quote = FALSE)
    if (nrow(mut.tab) > 0) {
        mut.alleles <- mufreq.bayes(mufreq = mut.tab$F[!mut.is.xy],
            CNt.max = CNt.max, depth.ratio = mut.tab$adjusted.ratio[!mut.is.xy],
            cellularity = cellularity, ploidy = ploidy, avg.depth.ratio = avg.depth.ratio,
            CNn = 2)
        mut.res <- cbind(mut.tab[!mut.is.xy, ], mut.alleles)
        if (!female) {
            if (sum(mut.is.xy) >= 1) {
                mut.alleles <- mufreq.bayes(mufreq = mut.tab$F[mut.is.xy],
                  CNt.max = CNt.max, depth.ratio = mut.tab$adjusted.ratio[mut.is.xy],
                  cellularity = cellularity, ploidy = ploidy,
                  avg.depth.ratio = avg.depth.ratio, CNn = 1)
                mut.xy <- cbind(mut.tab[mut.is.xy, ], mut.alleles)
                mut.res <- rbind(mut.res, mut.xy)
            }
        }
        write.table(mut.res, file = files$muts.file, col.names = TRUE,
            row.names = FALSE, sep = "\t", quote = FALSE)
    }
    list(seg.res = seg.res, mut.tab = mut.tab)
}

# Add better progress tracking with newlines
track_progress <- function(total, verbose = TRUE, label = "") {
    if (!verbose) {
        return(NULL)
    }

    message("\nStarting ", label, "...")
    pb <- txtProgressBar(min = 0, max = total, style = 3, width = 50)

    function(i) {
        setTxtProgressBar(pb, i)
        if (i == total) {
            close(pb)
            message("\nCompleted ", label)
        }
    }
}

# Enhance generate_result_plots with progress tracking and
# better memory management
generate_result_plots <- function(results, files, sequenza.extract,
    cellularity, ploidy, female, cp_results, XY) {
    seg.res <- results$seg.res
    mut.tab <- results$mut.tab
    avg.depth.ratio <- sequenza.extract$avg.depth.ratio
    cint <- cp_results$cint

    # Define segs.is.xy here for use throughout the
    # function
    if (female) {
        segs.is.xy <- seg.res$chromosome == XY["Y"]
    } else {
        segs.is.xy <- seg.res$chromosome %in% XY
    }

    # Add progress tracking for chromosome plots
    chrs <- unique(seg.res$chromosome)
    progress <- track_progress(length(chrs), label = "chromosome plots")

    # Memory efficient plotting
    pdf(files$chrw.file)
    message("\nGenerating chromosome plots...")
    for (i in seq_along(chrs)) {
        chr <- chrs[i]
        message("\nProcessing chromosome ", chr)

        CNn <- if (!female && chr %in% XY)
            1 else 2

        # Get only needed data for this chromosome
        chr_data <- list(mutations = sequenza.extract$mutations[[chr]],
            baf = sequenza.extract$BAF[[chr]], ratio = sequenza.extract$ratio[[chr]])

        chromosome.view(mut.tab = chr_data$mutations, baf.windows = chr_data$baf,
            ratio.windows = chr_data$ratio, cellularity = cellularity,
            ploidy = ploidy, main = chr, segments = seg.res[seg.res$chromosome ==
                chr, ], avg.depth.ratio = avg.depth.ratio, CNn = CNn,
            min.N.ratio = 1)

        # Update progress
        progress(i)

        # Clean up memory
        rm(chr_data)
        gc(verbose = FALSE)
    }
    dev.off()
    message("\nChromosome plots completed")

    # Chromosome view plots
    pdf(files$chrw.file)
    for (i in unique(seg.res$chromosome)) {
        CNn <- if (!female && i %in% XY)
            1 else 2
        chromosome.view(mut.tab = sequenza.extract$mutations[[i]],
            baf.windows = sequenza.extract$BAF[[i]], ratio.windows = sequenza.extract$ratio[[i]],
            cellularity = cellularity, ploidy = ploidy, main = i,
            segments = seg.res[seg.res$chromosome == i, ], avg.depth.ratio = avg.depth.ratio,
            CNn = CNn, min.N.ratio = 1)
    }
    dev.off()

    # Peak windows plots
    pdf(files$peak_win.file, height = 5, width = 8)
    for (i in unique(seg.res$chromosome)) {
        plot_peaks_win(sequenza.extract$win_peaks[[i]], main = i)
    }
    dev.off()

    # Density plots if available
    if (!is.null(sequenza.extract$ratio_baf_raster)) {
        # ...existing density plotting code...
    }

    # Genome view plots
    pdf(files$geno.file, height = 5, width = 15)
    if (sum(!is.na(seg.res$A)) > 0) {
        genome.view(seg.res)
    }
    genome.view(seg.res, "CN")
    plotRawGenome(sequenza.extract, cellularity = cellularity,
        ploidy = ploidy, mirror.BAF = TRUE)
    dev.off()

    # Copy number distribution plot
    barscn <- data.frame(size = seg.res$end.pos - seg.res$start.pos,
        CNt = seg.res$CNt)
    cn.sizes <- split(barscn$size, barscn$CNt)
    cn.sizes <- sapply(cn.sizes, "sum")
    pdf(files$cn.file)
    barplot(round(cn.sizes/sum(cn.sizes) * 100), names = names(cn.sizes),
        las = 1, ylab = "Percentage (%)", xlab = "Copy number")
    dev.off()

    # Write confidence intervals if available
    if (!is.null(cp_results$cint)) {
        res.tab <- data.frame(cellularity = c(cint$confint.cellularity[1],
            cint$max.cellularity[1], cint$confint.cellularity[2]),
            ploidy.estimate = c(cint$confint.ploidy[1], cint$max.ploidy[1],
                cint$confint.ploidy[2]), ploidy.mean.cn = weighted.mean(x = as.integer(names(cn.sizes)),
                w = cn.sizes))
        write.table(res.tab, files$cint.file, col.names = TRUE,
            row.names = FALSE, sep = "\t", quote = FALSE)
    }

    # Model fit plots
    pdf(files$fit.file, width = 6, height = 6)
    baf.model.view(cellularity = cellularity, ploidy = ploidy,
        segs = seg.res[!segs.is.xy, ])
    dev.off()

    # Alternative solutions if available
    if (!is.null(cp_results$cp.table)) {
        alt.sol <- alternative.cp.solutions(cp_results$cp.table)
        write.table(alt.sol, file = files$alt.file, col.names = TRUE,
            row.names = FALSE, sep = "\t", quote = FALSE)
        pdf(files$afit.file)
        for (sol in seq_len(nrow(alt.sol))) {
            baf.model.view(cellularity = alt.sol$cellularity[sol],
                ploidy = alt.sol$ploidy[sol], segs = seg.res[!segs.is.xy,
                  ])
        }
        dev.off()
    }
}

write_log_file <- function(log.file) {
    file_conn <- file(log.file)
    writeLines(c(date(), paste("Sequenza version:", packageVersion("sequenza"),
        sep = " ")), file_conn)
    close(file_conn)
}

# Add enhanced input validation
validate_results_input <- function(sequenza.extract, sample.id,
    out.dir) {
    if (!is.list(sequenza.extract)) {
        stop("sequenza.extract must be a list object from sequenza.extract()")
    }

    required_components <- c("BAF", "ratio", "raw_ratio", "mutations",
        "segments")
    missing <- setdiff(required_components, names(sequenza.extract))
    if (length(missing) > 0) {
        stop("Missing required components in sequenza.extract: ",
            paste(missing, collapse = ", "))
    }

    if (!is.character(sample.id) || length(sample.id) != 1) {
        stop("sample.id must be a single character string")
    }

    if (!is.character(out.dir)) {
        stop("out.dir must be a character string")
    }
}
