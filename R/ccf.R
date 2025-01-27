

# mean_0 = 2 mean_1 = 2 mean_2 = 3

# cellularity = 0.6

# n = 10000 w = 0.5

# normal_dist = rpois( round(n * (1 - cellularity), 0), mean_0 ) dist_1 =
# rpois( round(n * callularity * (1 - w), 0), mean_1 ) dist_2 = rpois( round(n
# * callularity * (w), 0), mean_2 )




# sim_subclonal_segs <- function(CNn, CNt, CNst, N, cellularity, ploidy, CCF,
# depth, ...) {

# dist_1 <- rpois( round(N * cellularity * (CCF), 0), CNt * depth )

# if (cellularity < 1) { normal_dist <- rpois( round(N * (1 - cellularity), 0),
# CNn * depth ) } else { normal_dist <- NA } if (CCF < 1) { dist_2 <- rpois(
# round(N * cellularity * (1 - CCF), 0), CNst * depth ) } else { dist_2 <- NA }
# plot( density(dist_1, adjust = CNt), ...  ) if (sum( is.na(normal_dist) == 0
# )) { lines(density(normal_dist, adjust = CNn)) } if (sum(is.na(dist_2)) == 0)
# { lines(density(dist_2, adjust = CNst)) } mean( c(dist_1, dist_2,
# normal_dist)/(ploidy * depth), na.rm = T )

# }

# par(mfrow = c(4, 1)) for (ccf in seq(1, 0.1, -0.3)) {

# m <- sim_subclonal_segs( 2, 3, 4, 1000, 1, 2, ccf, 60, xlim = c(0, 350), ylim
# = c(0, 0.04), main = paste('simulated ccf', ccf), xlab = 'coverage' )
# message(m) }


depth_ratio_dpois_ccf <- function(x, CNn, CNt, CNst, cellularity, ploidy, CCF, depth,
    avg.depth.ratio = 1, normal.ploidy = 2, ...) {
    dr_t <- theoretical.depth.ratio(CNt = CNt, cellularity = cellularity, ploidy = ploidy,
        CNn = CNn, normal.ploidy = normal.ploidy, avg.depth.ratio = avg.depth.ratio)
    dr_s <- theoretical.depth.ratio(CNt = CNst, cellularity = cellularity, ploidy = ploidy,
        CNn = CNn, normal.ploidy = normal.ploidy, avg.depth.ratio = avg.depth.ratio)
    dr_x <- weighted.mean(c(dr_t, dr_s), c(CCF, 1 - CCF))
    depth_x <- round(x * depth, 0)
    dpois(depth_x, round(dr_x * depth, 0), ...)
}

depth_ratio_ccf_likelihood <- function(x, CCFs, CNn, CNt, CNst, cellularity, ploidy,
    depth, avg.depth.ratio = 1, normal.ploidy = 2) {
    sapply(CCFs, FUN = function(x, r, CNn, CNt, CNst, cellularity, ploidy, depth,
        avg.depth.ratio, normal.ploidy) {
        depth_ratio_dpois_ccf(r, CNn, CNt, CNst, cellularity, ploidy, x, depth, avg.depth.ratio,
            normal.ploidy)

    }, r = x, CNn = CNn, CNt = CNt, CNst = CNst, cellularity, ploidy = ploidy, depth = depth,
        avg.depth.ratio = avg.depth.ratio, normal.ploidy = normal.ploidy)

}

# plot( seq(0, 1, 0.001), depth_ratio_ccf_likelihood( x = 2.2, CCFs = seq(0, 1,
# 0.001), CNn = 2, CNt = 4, CNst = 5, cellularity = 1, ploidy = 2, depth = 100
# ), type = 'l' )
