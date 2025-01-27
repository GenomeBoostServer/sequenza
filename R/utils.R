rmvnorm <- function(n, mean, sigma) {
    .Call(`_sequenza_rmvnorm`, n, mean, sigma)
}

rwishart <- function(df, S) {
    .Call(`_sequenza_rwishart`, df, S)
}

ginv <- function(m, tol = NULL) {
    .Call(`_sequenza_ginv`, m, tol)
}

sample_z_loop <- function(N, A, alpha, z, mu_ar, S_ar, mu_a0, R_a0, beta_a0, W_a0,
    nclust) {
    res <- .Call(`_sequenza_sample_z_loop`, N, A, alpha, z, mu_ar, S_ar, mu_a0, R_a0,
        beta_a0, W_a0, nclust)
    res$z <- as.vector(res$z)
    return(res)
}

ClusterTraces <- setRefClass("ClusterTraces", fields = c(traces = "matrix", traces.z = "matrix",
    tot_labels = "numeric"), methods = c(cluster = function(discard = 0.1) {
    .Call(`_sequenza_get_clusters`, traces.z, tot_labels, discard)
}, cluster_stats = function(discard = 0.1) {
    st <- .Call(`_sequenza_get_clusters_stats`, traces.z, tot_labels, discard)
    colnames(st) <- 1:ncol(st)
    return(st)
}))

ellipse <- function(mu, sigma, alpha = 0.05, npoints = 250, newplot = FALSE, draw = TRUE,
    ...) {
    es <- eigen(sigma)
    e1 <- es$vec %*% diag(sqrt(es$val))
    r1 <- sqrt(qchisq(1 - alpha, 2))
    theta <- seq(0, 2 * pi, len = npoints)
    v1 <- cbind(r1 * cos(theta), r1 * sin(theta))
    pts = t(mu - (e1 %*% t(v1)))
    if (newplot && draw) {
        plot(pts, ...)
    } else if (!newplot && draw) {
        lines(pts, ...)
    }
    invisible(pts)
}

pairs_ellipses <- function(A, z, iter, mu_a0, W_a0, mu_ar, active_clusts, S_ar) {
    n_size <- nrow(A)
    par(mfrow = c(n_size, n_size), mar = c(0.5, 0.5, 0.5, 0.5), oma = c(2, 2, 6,
        2))

    for (i in 1:n_size) {
        for (j in 1:n_size) {
            if (j > i) {
                plot(t(A[c(j, i), ]), col = z, xlab = "", ylab = "", xaxt = "n",
                  yaxt = "n")
                if (i == 1) {
                  axis(3, las = 1)
                }
                if (j == n_size) {
                  axis(4, las = 1)
                }
                for (k in active_clusts) {
                  ellipse(mu = mu_a0[c(j, i)], sigma = W_a0[c(i, j), c(i, j)], alpha = 0.25,
                    lwd = 5, npoints = 250, col = "black")
                  ellipse(mu = mu_ar[c(j, i), k], sigma = solve(S_ar[c(i, j), c(i,
                    j), k]), alpha = 0.25, lwd = 3, npoints = 250, col = k)
                }
            } else if (i == j) {
                plot(t(A[c(j, i), ]), type = "n", xlab = "", ylab = "", xaxt = "n",
                  yaxt = "n")
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
    mtext(paste("iter #", iter, sep = ""), side = 3, line = 3, outer = TRUE)
}
