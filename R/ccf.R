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
sim_subclonal_segs <- function(CNn, CNt, CNst, N, cellularity,
    ploidy, CCF, depth, ...) {
    # Generate distribution for primary tumor population
    dist_1 <- rpois(round(N * cellularity * CCF, 0), CNt * depth)

    # Generate distribution for normal cells if not pure
    # tumor
    if (cellularity < 1) {
        normal_dist <- rpois(round(N * (1 - cellularity), 0),
            CNn * depth)
    } else {
        normal_dist <- NA
    }

    # Generate distribution for subclonal population if
    # present
    if (CCF < 1) {
        dist_2 <- rpois(round(N * cellularity * (1 - CCF), 0),
            CNst * depth)
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
depth_ratio_dpois_ccf <- function(x, CNn, CNt, CNst, cellularity,
    ploidy, CCF, depth, avg.depth.ratio = 1, normal.ploidy = 2,
    ...) {
    # Calculate theoretical depth ratios
    dr_t <- theoretical.depth.ratio(CNt = CNt, cellularity = cellularity,
        ploidy = ploidy, CNn = CNn, normal.ploidy = normal.ploidy,
        avg.depth.ratio = avg.depth.ratio)

    dr_s <- theoretical.depth.ratio(CNt = CNst, cellularity = cellularity,
        ploidy = ploidy, CNn = CNn, normal.ploidy = normal.ploidy,
        avg.depth.ratio = avg.depth.ratio)

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
depth_ratio_ccf_likelihood <- function(x, CCFs, CNn, CNt, CNst,
    cellularity, ploidy, depth, avg.depth.ratio = 1, normal.ploidy = 2) {
    sapply(CCFs, FUN = function(x, r, CNn, CNt, CNst, cellularity,
        ploidy, depth, avg.depth.ratio, normal.ploidy) {
        depth_ratio_dpois_ccf(r, CNn, CNt, CNst, cellularity,
            ploidy, x, depth, avg.depth.ratio, normal.ploidy)
    }, r = x, CNn = CNn, CNt = CNt, CNst = CNst, cellularity = cellularity,
        ploidy = ploidy, depth = depth, avg.depth.ratio = avg.depth.ratio,
        normal.ploidy = normal.ploidy)
}

#' Calculate Cancer Cell Fraction (CCF) per Segment
#'
#' This function calculates the Cancer Cell Fraction (CCF) for each segment in the given segments table.
#'
#' @param segs_tab A data frame containing segment information with the following columns:
#'   - chromosome: integer, chromosome number
#'   - start.pos: integer, start position of the segment
#'   - end.pos: integer, end position of the segment
#'   - Bf: numeric, B allele frequency
#'   - N.BAF: numeric, number of B allele frequency data points
#'   - sd.BAF: numeric, standard deviation of B allele frequency
#'   - depth.ratio: numeric, depth ratio
#'   - N.ratio: numeric, number of depth ratio data points
#'   - sd.ratio: numeric, standard deviation of depth ratio
#'   - CNt: integer, total copy number
#'   - A: integer, copy number of allele A
#'   - B: integer, copy number of allele B
#'   - LPP: numeric, log posterior probability
#' @param CCFs A numeric vector of possible CCF values.
#' @param cellularity A numeric value representing the cellularity of the sample.
#' @param ploidy A numeric value representing the ploidy of the sample.
#' @param avg_depth A numeric value representing the average sequencing depth.
#' @param CNn A numeric value representing the copy number in normal cells.
#' @param avg.depth.ratio A numeric value representing the average depth ratio.
#' @param CNst_diff A numeric value representing the difference in copy number states (default is 1).
#' @param normal.ploidy A numeric value representing the normal ploidy (default is 2).
#' @param ... Additional arguments (currently not used).
#'
#' @return A data frame with the original segment information and an additional column 'CCF' containing the calculated CCF values for each segment.
#'
#' @examples
#' # Example usage:
#' segs_tab <- data.frame(
#'   chromosome = c(1, 1),
#'   start.pos = c(10000, 20000),
#'   end.pos = c(15000, 25000),
#'   Bf = c(0.5, 0.6),
#'   N.BAF = c(0.4, 0.5),
#'   sd.BAF = c(0.01, 0.02),
#'   depth.ratio = c(1.2, 1.3),
#'   N.ratio = c(1.1, 1.2),
#'   sd.ratio = c(0.05, 0.06),
#'   CNt = c(2, 3),
#'   A = c(1, 1),
#'   B = c(1, 2),
#'   LPP = c(0.9, 0.8)
#' )
#' CCFs <- seq(0, 1, by = 0.1)
#' cellularity <- 0.8
#' ploidy <- 2
#' avg_depth <- 30
#' CNn <- 2
#' avg.depth.ratio <- 1.1
#' result <- ccf_per_segment(segs_tab, CCFs, cellularity, ploidy, avg_depth, CNn, avg.depth.ratio)
#' print(result)
#'
#' @export
ccf_per_segment <- function(segs_tab, CCFs, cellularity, ploidy,
    avg_depth, CNn, avg.depth.ratio, CNst_diff = 1, normal.ploidy = 2,
    ...) {
    # Generate BAF/ratio model points
    types_matrix <- baf.types.matrix(0, max(segs_tab$CNt) + CNst_diff,
        CNn)
    types_values <- baf.model.points(cellularity, ploidy, types_matrix,
        avg.depth.ratio)
    types_data <- cbind(types_matrix, types_values)

    # Vectorized operations
    CNt_values <- segs_tab$CNt
    depth_ratios <- segs_tab$depth.ratio
    depths <- avg_depth * depth_ratios

    # Get type ratios for each CNt value
    type_ratios <- sapply(CNt_values, function(cnt) {
        ratios <- types_data$depth.ratio[types_data$CNt == cnt]
        if (length(ratios) > 0) {
            ratios[1]
        } else {
            NA
        }
    })

    # Calculate subclonal CNt values
    i_sCNt <- ifelse(type_ratios >= depth_ratios, CNt_values -
        CNst_diff, CNt_values + CNst_diff)

    # Calculate CCF likelihoods for all segments
    ccf_results <- t(mapply(function(ratio, cnt, scnt, depth) {
        likelihood <- depth_ratio_ccf_likelihood(x = ratio, CCFs = CCFs,
            CNn = CNn, CNt = cnt, CNst = scnt, cellularity = cellularity,
            ploidy = ploidy, depth = depth, avg.depth.ratio = avg.depth.ratio,
            normal.ploidy = normal.ploidy)
        c(CCF = CCFs[which.max(likelihood)], CNst = scnt, llp = -log(max(likelihood)))
    }, depth_ratios, CNt_values, i_sCNt, depths))

    return(ccf_results)
}

# Example usage for CCF likelihood plot: plot( seq(0, 1,
# 0.001), depth_ratio_ccf_likelihood( x = 2.2, CCFs =
# seq(0, 1, 0.001), CNn = 2, CNt = 4, CNst = 5, cellularity
# = 1, ploidy = 2, depth = 100 ), type = 'l' )

# Example usage for subclonal segment simulation: par(mfrow
# = c(4, 1)) for (ccf in seq(1, 0.1, -0.3)) { m <-
# sim_subclonal_segs( CNn = 2, CNt = 3, CNst = 4, N = 1000,
# cellularity = 1, ploidy = 2, CCF = ccf, depth = 60, xlim
# = c(0, 350), ylim = c(0, 0.04), main = paste('simulated
# ccf', ccf), xlab = 'coverage' ) message(m) }
