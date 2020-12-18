find_breaks <- function(seqz.baf, slide_win, peak_win, arms,
    chr_name, verbose) {
    chromosome <- gsub(x = seqz.baf$chromosome,
        pattern = "chr", replacement = "")
    chromosome <- paste0("chr", chromosome)
    if (verbose) {
      message("Segmenting depth ratios")
    }
    ratio_diffs <- slide_matrix(seqz.baf$adjusted.ratio, w = slide_win,
                                position = seqz.baf$position, verbose = verbose)
    if (verbose) {
      message("Segmenting allele frequencies")
    }
    bf_diffs <- slide_matrix(seqz.baf$Bf, w = slide_win,
                             position = seqz.baf$position, verbose = verbose)

    peaks_both <- get_gaps_peaks(x = (ratio_diffs$y + bf_diffs$y) / 2,
                                 position = ratio_diffs$x, w = peak_win,
                                 arms = arms)
    breaks <- lapply(peaks_both, FUN = function(peaks) {
         coords <- peaks
         if (is.null(coords)) {
           NULL
         } else {
           pos_start <- coords[-length(coords)]
           pos_end <- coords[-1]
           data.frame(
             chrom = chr_name, start.pos = pos_start,
             end.pos = pos_end)
         }
    })
    breaks <- do.call(rbind, breaks)
    not.uniq <- which(breaks$end.pos == c(breaks$start.pos[-1], 0))
    breaks$end.pos[not.uniq] <- breaks$end.pos[not.uniq] - 1
    breaks
}

slide_tracks <- function(seqz.baf, slide_win,
    signal_out = c("both", "ratio", "baf"), verbose = TRUE) {
    signal_out <- match.arg(arg = signal_out,
        choices = signal_out)
    chromosome <- gsub(x = seqz.baf$chromosome,
        pattern = "chr", replacement = "")
    chromosome <- paste0("chr", chromosome)
    if (signal_out %in% c("both", "ratio")) {
        if (verbose) {
            message("Segmenting depth ratios")
        }
        ratio_diffs <- slide_matrix(
            seqz.baf$adjusted.ratio, w = slide_win,
            position = seqz.baf$position, verbose = verbose)
    }

    if (signal_out %in% c("both", "baf")) {
        if (verbose) {
            message("Segmenting allele frequencies")
        }
        bf_diffs <- slide_matrix(
            seqz.baf$Bf, w = slide_win,
        position = seqz.baf$position, verbose = verbose)
    }

    if (signal_out == "both") {
        data.frame(y = (ratio_diffs$y + bf_diffs$y) / 2,
            x = ratio_diffs$x)
    } else if (signal_out == "baf") {
        bf_diffs
    } else {
        ratio_diffs
    }
}

peaks_tracks <- function(diff_track, peak_win, arms,
    chr_name, verbose) {
    peaks <- get_gaps_peaks(
        x = diff_track$y, position = diff_track$x,
        w = peak_win, arms = arms)

    breaks <- lapply(peaks, FUN = function(peaks) {
         coords <- peaks
         if (is.null(coords)) {
           NULL
         } else {
           pos_start <- coords[-length(coords)]
           pos_end <- coords[-1]
           data.frame(
             chrom = chr_name, start.pos = pos_start,
             end.pos = pos_end)
         }
    })
    breaks <- do.call(rbind, breaks)
    not_uniq <- which(breaks$end.pos == c(breaks$start.pos[-1], 0))
    breaks$end.pos[not_uniq] <- breaks$end.pos[not_uniq] - 1
    breaks
}


extract_breaks <- function(data, data_het, breaks,
                           slide_win, peak_win, assembly, chromosome,
                           verbose = TRUE) {
  if (is.null(breaks)) {
      golden_path <- paste("http://hgdownload.cse.ucsc.edu",
                           "goldenPath", assembly, "database",
                           "cytoBandIdeo.txt.gz", sep = "/")
      arms <- get_assembly(url = golden_path, prefix = "chr")
      chr_arm <- gsub(x = chromosome,
                      pattern = "chr", replacement = "")
      chr_arm <- paste0("chr", chr_arm)

      arms_i <- arms[arms$chromosome == chr_arm, ]
      data_het <- data_het[data_het$chromosome == chromosome, ]
      find_breaks(data_het, slide_win, peak_win, arms_i, chromosome, verbose)
  } else {
      breaks
  }
}

extract_breaks_tracks <- function(track, breaks,
    slide_win, peak_win, assembly, chromosome) {
    if (is.null(breaks)) {
        golden_path <- paste("http://hgdownload.cse.ucsc.edu",
            "goldenPath", assembly, "database",
            "cytoBand.txt.gz", sep = "/")
        arms <- get_assembly(url = golden_path, prefix = "chr")
        chr_arm <- gsub(x = chromosome,
             pattern = "chr", replacement = "")
        chr_arm <- paste0("chr", chr_arm)

        arms_i <- arms[arms$chromosome == chr_arm, ]
        peaks_tracks(track, peak_win, arms_i, chromosome)
  } else {
      breaks
  }
}


segment.breaks <- function(seqz.tab, breaks, min.reads.baf = 1,
    weighted.mean = TRUE) {
    if (weighted.mean){
        w.r <- sqrt(seqz.tab$depth.normal)
        rw <- seqz.tab$adjusted.ratio * w.r
        w.b <- sqrt(seqz.tab$good.reads)
        bw <- seqz.tab$Bf * w.b
        seqz.tab <- cbind(seqz.tab[, c("chromosome", "position",
            "zygosity.normal", "good.reads", "Af", "Bf")],
            rw = rw, w.r = w.r, bw = bw, w.b = w.b)
    }
    chr.order <- unique(seqz.tab$chromosome)
    seqz.tab <- split(seqz.tab, f = seqz.tab$chromosome)
    segments <- list()
    for (i in 1:length(seqz.tab)) {
        seqz.b.i <- seqz.tab[[i]][seqz.tab[[i]]$zygosity.normal == "het", ]
        seqz.b.i <- seqz.b.i[seqz.b.i$good.reads >= min.reads.baf, ]
        breaks.i <- breaks[breaks$chrom == names(seqz.tab)[i], ]
        nb <- nrow(breaks.i)
        breaks.vect <- do.call(cbind, split.data.frame(breaks.i[,
            c("start.pos", "end.pos")], f = 1:nb))
        unique.breaks <- function(b, offset = 1) {
            while(any(diff(b) == 0)) {
                b[which(diff(b) == 0) + 1] <- b[diff(b) == 0] + offset
            }
            b
        }
        breaks.vect <- unique.breaks(b = as.numeric(breaks.vect), offset = 1)
        fact.r.i <- cut(seqz.tab[[i]]$position, breaks.vect)
        fact.b.i <- cut(seqz.b.i$position, breaks.vect)
        seg.i.s.r <- sapply(X = split(seqz.tab[[i]]$chromosome,
            f = fact.r.i), FUN = length)
        seg.i.s.b <- sapply(X = split(seqz.b.i$chromosome,
            f = fact.b.i), FUN = length)

        if (weighted.mean) {
            seg.i.rw    <- sapply(X = split(seqz.tab[[i]]$rw, f = fact.r.i),
                FUN = function(a) sum(a, na.rm = TRUE))
            seg.i.w.r   <- sapply(X = split(seqz.tab[[i]]$w.r, f = fact.r.i),
                FUN = function(a) sum(a, na.rm = TRUE))
            seg.i.r.sd  <- sapply(X = split(seqz.tab[[i]]$rw /
                seqz.tab[[i]]$w.r, f = fact.r.i),
                FUN = function(a) sd(a, na.rm = TRUE))
            seg.i.b.sd  <- sapply(X = split(seqz.b.i$bw /
                seqz.b.i$w.b, f = fact.b.i),
                FUN = function(a) sd(a, na.rm = TRUE))
            A.split <- split(seqz.b.i$Af, f = fact.b.i)
            B.split <- split(seqz.b.i$Bf, f = fact.b.i)
            d.split <- split(seqz.b.i$good.reads, f = fact.b.i)
            window.quantiles <- mapply(b_allele_freq, Af = A.split,
                Bf = B.split, good.reads = d.split, conf = 0.95)
            segments.i <- data.frame(chromosome  = names(seqz.tab)[i],
                start.pos = as.numeric(breaks.vect[-length(breaks.vect)]),
                end.pos = as.numeric(breaks.vect[-1]),
                Bf = window.quantiles[2, ], N.BAF = seg.i.s.b,
                sd.BAF = seg.i.b.sd, depth.ratio = seg.i.rw / seg.i.w.r,
                N.ratio = seg.i.s.r, sd.ratio = seg.i.r.sd,
                stringsAsFactors = FALSE)
        } else {
            seg.i.r    <- sapply(X = split(seqz.tab[[i]]$adjusted.ratio,
                f = fact.r.i), FUN = function(a) mean(a, na.rm = TRUE))
            A.split <- split(seqz.b.i$Af, f = fact.b.i)
            B.split <- split(seqz.b.i$Bf, f = fact.b.i)
            d.split <- split(seqz.b.i$good.reads, f = fact.b.i)
            window.quantiles <- mapply(b_allele_freq, Af = A.split,
                Bf = B.split, good.reads = d.split, conf = 0.95)
            seg.i.r.sd <- sapply(X = split(seqz.tab[[i]]$adjusted.ratio,
                f = fact.r.i), FUN = function(a) sd(a, na.rm = TRUE))
            seg.i.b.sd <- sapply(X = split(seqz.b.i$Bf, f = fact.b.i),
                FUN = function(a) sd(a, na.rm = TRUE))
            segments.i <- data.frame(chromosome  = names(seqz.tab)[i],
                start.pos = as.numeric(breaks.vect[-length(breaks.vect)]),
                end.pos = as.numeric(breaks.vect[-1]),
                Bf = window.quantiles[2, ], N.BAF = seg.i.s.b,
                sd.BAF = seg.i.b.sd, depth.ratio = seg.i.r,
                N.ratio = seg.i.s.r, sd.ratio = seg.i.r.sd,
                stringsAsFactors = FALSE)
        }
        segments[[i]] <- segments.i[seq(from = 1,
            to = nrow(segments.i), by = 2),]
    }
    segments <- do.call(rbind, segments[as.factor(chr.order)])
    row.names(segments) <- 1:nrow(segments)
    len.seg <- (segments$end.pos - segments$start.pos) / 1e6
    segments[(segments$N.ratio / len.seg) >= 2, ]
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
    sum(is_similar, na.rm = TRUE) / length(na.exclude(is_similar))
}