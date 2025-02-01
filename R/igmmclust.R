ClusterTraces <- setRefClass("ClusterTraces", fields = c(traces = "matrix",
    traces.z = "matrix", tot_labels = "numeric"), methods = c(cluster = function(discard = 0.1) {
    .Call(`_sequenza_get_clusters`, traces.z, tot_labels, discard)
}, cluster_stats = function(discard = 0.1) {
    st <- .Call(`_sequenza_get_clusters_stats`, traces.z, tot_labels,
        discard)
    colnames(st) <- 1:ncol(st)
    return(st)
}))

gibbs <- function(A, z_init, iters = 100, progressbar = TRUE,
    plots = FALSE, plots.freq = 100) {
    D <- dim(A)[1]
    N <- dim(A)[2]
    init_clust <- length(unique(z_init))
    # traces
    traces.z <- matrix(NA, iters, N)
    traces <- matrix(NA, iters, 6)
    colnames(traces) <- c("n_clusters", "alpha", "mean_mu_a0",
        "det(R_a0)", "det(W_a0)", "beta_a0")

    # data mean and covariance
    Sigma_a <- cov(t(A))
    Lambda_a <- solve(Sigma_a)
    mu_a <- rowMeans(A)

    # init state
    alpha <- 2
    z <- z_init
    mu_a0 <- mu_a
    R_a0 <- solve(Sigma_a)
    W_a0 <- Sigma_a
    beta_a0 <- dim(A)[1]
    S_ar <- array(NA, dim = c(D, D, 100))
    mu_ar <- matrix(NA, D, 100)
    for (k in 1:length(unique(z))) {
        mu_ar[, k] <- rowMeans(A[, z == k])
    }
    if (progressbar) {
        pb <- timerProgressBar(min = 1, max = iters, initial = 1,
            style = 3, char = "+", width = 50)
    }
    for (i in 1:iters) {
        if (progressbar) {
            setTxtProgressBar(pb, i)
        }
        # Active components
        K <- length(unique(z))

        # Sample component parameters
        active_clusts <- sort(unique(z))
        for (k in active_clusts) {
            mask <- (z == k)
            S_ar[, , k] <- sample_S_ar(A[, mask, drop = FALSE],
                mu_ar[, k, drop = FALSE], beta_a0, W_a0)
            mu_ar[, k] <- sample_mu_ar(A[, mask, drop = FALSE],
                S_ar[, , k], mu_a0, R_a0)
        }

        # Sample hyperpriors
        mu_a0 <- sample_mu_a0(Lambda_a, mu_a, mu_ar[, active_clusts,
            drop = FALSE], R_a0)
        R_a0 <- sample_R_a0(Sigma_a, mu_ar[, active_clusts, drop = FALSE],
            mu_a0)
        W_a0 <- sample_W_a0(Lambda_a, S_ar[, , active_clusts,
            drop = FALSE], beta_a0)
        beta_a0 <- sample_beta_a0(S_ar[, , active_clusts, drop = FALSE],
            W_a0)

        res <- sample_z_loop(N, A, alpha, z, mu_ar, S_ar, mu_a0,
            R_a0, beta_a0, W_a0, init_clust)
        z <- res$z
        mu_ar <- res$mu_ar
        S_ar <- res$S_ar
        init_clust <- res$nclust


        # Sample concentration parameter
        K <- length(unique(z))
        alpha <- sample_alpha(K = K, U = N)

        traces[i, ] <- c(length(unique(z)), alpha, mean(mu_a0),
            det(R_a0), det(W_a0), beta_a0)
        traces.z[i, ] <- z
        if (plots && (i%%plots.freq == 0)) {
            pairs_ellipses(A, z, i, mu_a0, W_a0, mu_ar, active_clusts,
                S_ar)
        }
    }
    if (progressbar) {
        closepb(pb)
    }

    return(ClusterTraces$new(traces = traces, traces.z = traces.z,
        tot_labels = init_clust))
}

ellipse <- function(mu, sigma, alpha = 0.05, npoints = 250, newplot = FALSE,
    draw = TRUE, ...) {
    es <- eigen(sigma)
    e1 <- es$vec %*% diag(sqrt(es$val))
    r1 <- sqrt(qchisq(1 - alpha, 2))
    theta <- seq(0, 2 * pi, len = npoints)
    v1 <- cbind(r1 * cos(theta), r1 * sin(theta))
    pts <- t(mu - (e1 %*% t(v1)))
    if (newplot && draw) {
        plot(pts, ...)
    } else if (!newplot && draw) {
        lines(pts, ...)
    }
    invisible(pts)
}

pairs_ellipses <- function(A, z, iter, mu_a0, W_a0, mu_ar, active_clusts,
    S_ar) {
    n_size <- nrow(A)
    par(mfrow = c(n_size, n_size), mar = c(0.5, 0.5, 0.5, 0.5),
        oma = c(2, 2, 6, 2))

    for (i in 1:n_size) {
        for (j in 1:n_size) {
            if (j > i) {
                plot(t(A[c(j, i), ]), col = z, xlab = "", ylab = "",
                  xaxt = "n", yaxt = "n")
                if (i == 1) {
                  axis(3, las = 1)
                }
                if (j == n_size) {
                  axis(4, las = 1)
                }
                for (k in active_clusts) {
                  ellipse(mu = mu_a0[c(j, i)], sigma = W_a0[c(i,
                    j), c(i, j)], alpha = 0.25, lwd = 5, npoints = 250,
                    col = "black")
                  ellipse(mu = mu_ar[c(j, i), k], sigma = solve(S_ar[c(i,
                    j), c(i, j), k]), alpha = 0.25, lwd = 3,
                    npoints = 250, col = k)
                }
            } else if (i == j) {
                plot(t(A[c(j, i), ]), type = "n", xlab = "",
                  ylab = "", xaxt = "n", yaxt = "n")
                text(mean(A[i, ]), mean(A[j, ]), rownames(A)[i])
                if (i == 1) {
                  axis(2, las = 1)
                }
                if (i == n_size) {
                  axis(1, las = 1)
                }
            } else {
                plot.new()
            }
        }
    }
    mtext(paste("iter #", iter, sep = ""), side = 3, line = 3,
        outer = TRUE)
}


cluster_segments <- function(bf, depth_ratio, init_clust = 100,
    ...) {
    x <- cbind(bf, depth_ratio)
    gibbs(t(x), z_init = sample(1:init_clust, nrow(x), replace = T),
        ...)
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
                breaks[[length(breaks) + 1]] <- c(segs$start.pos[i],
                  segs$end.pos[i])
                last_clust <- clusters[i]
            }
        }
    }
    do.call(rbind, lapply(breaks, FUN = function(x) c(min(x),
        max(x))))
}
