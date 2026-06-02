# Precompute window indices from positions and chromosomes.
# Returns an index structure that can be reused across multiple variables,
# avoiding repeated split(data.frame(...)) calls.
precompute_windows <- function(positions, chromosomes, window = 1e6,
    overlap = 0, start.coord = 1) {
    overlap <- as.integer(overlap)
    window.offset <- as.integer(window - round(window * (overlap / (overlap + 1))))
    chr.ordered <- unique(chromosomes)

    # Split only integer indices by chromosome (not data.frames)
    global_idx <- split(seq_along(positions), factor(chromosomes, levels = chr.ordered))

    chr_windows <- lapply(chr.ordered, function(chr) {
        g_idx <- global_idx[[chr]]
        pos <- positions[g_idx]

        range.pos <- range(pos, na.rm = TRUE)
        if (!is.null(start.coord)) {
            range.pos[1] <- as.integer(start.coord)
        }
        beam.coords <- seq.int(range.pos[1], range.pos[2], by = window.offset)
        if (max(beam.coords) != range.pos[2]) {
            beam.coords <- c(beam.coords, range.pos[2])
        }
        nWindows <- length(beam.coords) - overlap - 1
        if (nWindows < 1) return(list(starts = integer(0), ends = integer(0), window_idx = list()))
        pos.cut <- cut(pos, breaks = beam.coords)

        # Split global indices by window bin
        bin_idx <- split(g_idx, pos.cut)

        # Merge bins for overlap
        idx.list <- lapply(1:nWindows, function(ii) ii + (0:overlap))
        window_idx <- lapply(idx.list, function(idx) {
            unlist(bin_idx[idx], use.names = FALSE)
        })

        list(starts = beam.coords[1:nWindows],
             ends = beam.coords[(1:nWindows) + 1 + overlap],
             window_idx = window_idx)
    })
    names(chr_windows) <- chr.ordered
    chr_windows
}

# Compute weighted window statistics using precomputed indices.
# x: numeric vector (same length as positions used in precompute)
# weight: optional weight vector (sqrt applied internally)
# precomp: result of precompute_windows()
windowValues_fast <- function(x, precomp, weight = NULL) {
    if (!is.null(weight)) {
        w <- sqrt(weight)
    }
    lapply(precomp, function(cw) {
        if (is.null(weight)) {
            window.means <- vapply(cw$window_idx, function(idx) {
                mean(x[idx], na.rm = TRUE)
            }, numeric(1))
        } else {
            window.means <- vapply(cw$window_idx, function(idx) {
                weighted.mean(x[idx], w[idx], na.rm = TRUE)
            }, numeric(1))
        }
        window.quantiles <- vapply(cw$window_idx, function(idx) {
            quantile(x[idx], probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE)
        }, numeric(2))
        window.counts <- vapply(cw$window_idx, function(idx) {
            length(idx)
        }, integer(1))
        data.frame(start = cw$starts, end = cw$ends,
            mean = window.means, q0 = window.quantiles[1, ],
            q1 = window.quantiles[2, ], N = window.counts)
    })
}

# Compute BAF window statistics using precomputed indices.
windowBf_fast <- function(Af, Bf, good.reads, precomp, conf = 0.95) {
    lapply(precomp, function(cw) {
        window.quantiles <- vapply(cw$window_idx, function(idx) {
            b_allele_freq(Af[idx], Bf[idx], good.reads[idx], conf = conf)
        }, numeric(3))
        window.counts <- vapply(cw$window_idx, function(idx) {
            length(idx)
        }, integer(1))
        data.frame(start = cw$starts, end = cw$ends, mean = window.quantiles[2, ],
            q0 = window.quantiles[1, ], q1 = window.quantiles[3, ], N = window.counts)
    })
}

# Original windowValues (kept for backward compatibility)
windowValues <- function(x, positions, chromosomes, window = 1e+06,
    overlap = 0, weight = rep.int(x = 1, times = length(x)),
    start.coord = 1) {
    weight <- sqrt(weight)
    overlap <- as.integer(overlap)
    window.offset <- as.integer(window - round(window * (overlap/(overlap +
        1))))
    chr.ordered <- unique(chromosomes)
    data.splitByChr <- split(data.frame(pos = positions, x = x,
        weight = weight), f = factor(chromosomes, levels = chr.ordered))
    lapply(data.splitByChr, function(x) {
        range.pos <- range(x$pos, na.rm = TRUE)
        if (!is.null(start.coord)) {
            range.pos[1] <- as.integer(start.coord)
        }
        beam.coords <- seq.int(range.pos[1], range.pos[2], by = window.offset)
        if (max(beam.coords) != range.pos[2]) {
            beam.coords <- c(beam.coords, range.pos[2])
        }
        nWindows <- length(beam.coords) - overlap - 1
        pos.cut <- cut(x$pos, breaks = beam.coords)
        x.split <- split(x$x, f = pos.cut)
        weight.split <- split(x$weight, f = pos.cut)
        window.starts <- beam.coords[1:nWindows]
        window.ends <- beam.coords[(1:nWindows) + 1 + overlap]
        idx.list <- lapply(1:nWindows, function(ii) ii + (0:overlap))
        x.window <- lapply(idx.list, function(idx) {
            unlist(x.split[idx], use.names = FALSE)
        })
        weight.window <- lapply(idx.list, function(idx) {
            unlist(weight.split[idx], use.names = FALSE)
        })
        window.means <- mapply(weighted.mean, x = x.window, w = weight.window)
        window.quantiles <- sapply(x.window, quantile, probs = c(0.25,
            0.75), na.rm = TRUE, names = FALSE)
        window.counts <- sapply(x.window, length)
        data.frame(start = window.starts, end = window.ends,
            mean = window.means, q0 = window.quantiles[1, ],
            q1 = window.quantiles[2, ], N = window.counts)
    })
}

windowBf <- function(Af, Bf, good.reads, positions, chromosomes,
    window = 1e+06, overlap = 0, start.coord = 1, conf = 0.95) {
    overlap <- as.integer(overlap)
    window.offset <- as.integer(window - round(window * (overlap/(overlap +
        1))))
    chr.ordered <- unique(chromosomes)
    data.splitByChr <- split(data.frame(pos = positions, Bf,
        Af, good.reads), f = factor(chromosomes, levels = chr.ordered))
    lapply(data.splitByChr, function(data.oneChr) {
        range.pos <- range(data.oneChr$pos, na.rm = TRUE)
        if (!is.null(start.coord)) {
            range.pos[1] <- as.integer(start.coord)
        }
        beam.coords <- seq.int(range.pos[1], range.pos[2], by = window.offset)
        if (max(beam.coords) != range.pos[2]) {
            beam.coords <- c(beam.coords, range.pos[2])
        }
        nWindows <- length(beam.coords) - overlap - 1
        pos.cut <- cut(data.oneChr$pos, breaks = beam.coords)
        A.split <- split(data.oneChr$Af, f = pos.cut)
        B.split <- split(data.oneChr$Bf, f = pos.cut)
        d.split <- split(data.oneChr$good.reads, f = pos.cut)
        window.starts <- beam.coords[1:nWindows]
        window.ends <- beam.coords[(1:nWindows) + 1 + overlap]
        idx.list <- lapply(1:nWindows, function(ii) ii + (0:overlap))
        A.window <- lapply(idx.list, function(idx) {
            unlist(A.split[idx], use.names = FALSE)
        })
        B.window <- lapply(idx.list, function(idx) {
            unlist(B.split[idx], use.names = FALSE)
        })
        d.window <- lapply(idx.list, function(idx) {
            unlist(d.split[idx], use.names = FALSE)
        })
        window.quantiles <- mapply(b_allele_freq, Af = A.window,
            Bf = B.window, good.reads = d.window, conf = conf)
        window.counts <- sapply(B.window, length)
        data.frame(start = window.starts, end = window.ends,
            mean = window.quantiles[2, ], q0 = window.quantiles[1,
                ], q1 = window.quantiles[3, ], N = window.counts)
    })
}
