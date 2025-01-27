
#' Simulate subclonal segments in genomic data
#' 
#' @param CNn Copy number in normal cells
#' @param CNt Copy number in tumor cells
#' @param CNst Copy number in subclonal tumor cells
#' @param N Number of samples to simulate
#' @param cellularity Tumor purity (0-1)
#' @param ploidy Base ploidy of the genome
#' @param CCF Cancer cell fraction
#' @param depth Sequencing depth
#' @param ... Additional parameters passed to density plot
#' @return Mean normalized coverage
sim_subclonal_segs <- function(CNn, CNt, CNst, N, cellularity, ploidy, CCF, depth,
    ...) {
    # Generate distribution for primary tumor population
    dist_1 <- rpois(round(N * cellularity * CCF, 0), CNt * depth)

    # Generate distribution for normal cells if not pure tumor
    if (cellularity < 1) {
        normal_dist <- rpois(round(N * (1 - cellularity), 0), CNn * depth)
    } else {
        normal_dist <- NA
    }

    # Generate distribution for subclonal population if present
    if (CCF < 1) {
        dist_2 <- rpois(round(N * cellularity * (1 - CCF), 0), CNst * depth)
    } else {
        dist_2 <- NA
    }

    # Plot density distributions
    plot(density(dist_1, adjust = CNt), ...)

    if (sum(is.na(normal_dist) == 0)) {
        lines(density(normal_dist, adjust = CNn))
    }
    if (sum(is.na(dist_2)) == 0) {
        lines(density(dist_2, adjust = CNst))
    }

    # Return normalized mean coverage
    mean(c(dist_1, dist_2, normal_dist)/(ploidy * depth), na.rm = TRUE)
}

#' Calculate Poisson probability for depth ratio given CCF
#' 
#' @param x Observed depth ratio
#' @param CNn Copy number in normal cells
#' @param CNt Copy number in tumor cells
#' @param CNst Copy number in subclonal tumor cells
#' @param cellularity Tumor purity (0-1)
#' @param ploidy Base ploidy
#' @param CCF Cancer cell fraction
#' @param depth Sequencing depth
#' @param avg.depth.ratio Average depth ratio (default: 1)
#' @param normal.ploidy Ploidy of normal cells (default: 2)
#' @param ... Additional parameters passed to dpois
#' @return Poisson probability
depth_ratio_dpois_ccf <- function(x, CNn, CNt, CNst, cellularity, ploidy, CCF, depth,
    avg.depth.ratio = 1, normal.ploidy = 2, ...) {
    # Calculate theoretical depth ratios
    dr_t <- theoretical.depth.ratio(CNt = CNt, cellularity = cellularity, ploidy = ploidy,
        CNn = CNn, normal.ploidy = normal.ploidy, avg.depth.ratio = avg.depth.ratio)

    dr_s <- theoretical.depth.ratio(CNt = CNst, cellularity = cellularity, ploidy = ploidy,
        CNn = CNn, normal.ploidy = normal.ploidy, avg.depth.ratio = avg.depth.ratio)

    # Calculate weighted mean of depth ratios
    dr_x <- weighted.mean(c(dr_t, dr_s), c(CCF, 1 - CCF))
    depth_x <- round(x * depth, 0)

    dpois(depth_x, round(dr_x * depth, 0), ...)
}

#' Calculate likelihood for different CCF values
#' 
#' @param x Observed depth ratio
#' @param CCFs Vector of CCF values to test
#' @param CNn Copy number in normal cells
#' @param CNt Copy number in tumor cells
#' @param CNst Copy number in subclonal tumor cells
#' @param cellularity Tumor purity (0-1)
#' @param ploidy Base ploidy
#' @param depth Sequencing depth
#' @param avg.depth.ratio Average depth ratio (default: 1)
#' @param normal.ploidy Ploidy of normal cells (default: 2)
#' @return Vector of likelihood values
depth_ratio_ccf_likelihood <- function(x, CCFs, CNn, CNt, CNst, cellularity, ploidy,
    depth, avg.depth.ratio = 1, normal.ploidy = 2) {
    sapply(CCFs, FUN = function(x, r, CNn, CNt, CNst, cellularity, ploidy, depth,
        avg.depth.ratio, normal.ploidy) {
        depth_ratio_dpois_ccf(r, CNn, CNt, CNst, cellularity, ploidy, x, depth, avg.depth.ratio,
            normal.ploidy)
    }, r = x, CNn = CNn, CNt = CNt, CNst = CNst, cellularity = cellularity, ploidy = ploidy,
        depth = depth, avg.depth.ratio = avg.depth.ratio, normal.ploidy = normal.ploidy)
}

# Example usage for CCF likelihood plot: plot( seq(0, 1, 0.001),
# depth_ratio_ccf_likelihood( x = 2.2, CCFs = seq(0, 1, 0.001), CNn = 2, CNt =
# 4, CNst = 5, cellularity = 1, ploidy = 2, depth = 100 ), type = 'l' )

# Example usage for subclonal segment simulation: par(mfrow = c(4, 1)) for (ccf
# in seq(1, 0.1, -0.3)) { m <- sim_subclonal_segs( CNn = 2, CNt = 3, CNst = 4,
# N = 1000, cellularity = 1, ploidy = 2, CCF = ccf, depth = 60, xlim = c(0,
# 350), ylim = c(0, 0.04), main = paste('simulated ccf', ccf), xlab =
# 'coverage' ) message(m) }
