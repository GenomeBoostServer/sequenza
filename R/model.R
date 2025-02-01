#' Calculate theoretical depth ratio
#'
#' @param CNt Copy number in tumor
#' @param cellularity Tumor purity (0-1)
#' @param ploidy Tumor ploidy
#' @param CNn Normal copy number (default: 2)
#' @param normal.ploidy Normal tissue ploidy (default: 2)
#' @param avg.depth.ratio Average depth ratio (default: 1)
#' @return Theoretical depth ratio
theoretical.depth.ratio <- function(CNt, cellularity, ploidy,
    CNn = 2, normal.ploidy = 2, avg.depth.ratio = 1) {
    # Validate input parameters
    if (cellularity < 0 || cellularity > 1) {
        stop("Cellularity must be between 0 and 1")
    }

    # Calculate copy number term accounting for cellularity
    cellu_copy_term <- (1 - cellularity) + (CNt/CNn * cellularity)

    # Calculate ploidy term accounting for cellularity
    ploidy_cellu_term <- (ploidy/normal.ploidy * cellularity) +
        (1 - cellularity)

    # Return normalized depth ratio
    avg.depth.ratio * cellu_copy_term/ploidy_cellu_term
}

#' Calculate theoretical B-allele frequency
#'
#' @param CNt Copy number in tumor
#' @param B Number of B alleles
#' @param cellularity Tumor purity (0-1)
#' @param CNn Normal copy number (default: 2)
#' @return Vector of theoretical B-allele frequencies
theoretical.baf <- function(CNt, B, cellularity, CNn = 2) {
    # Validate input parameters
    if (cellularity < 0 || cellularity > 1) {
        stop("Cellularity must be between 0 and 1")
    }

    # Calculate B-allele frequency
    baf <- ((B * cellularity) + (1 - cellularity))/((CNt * cellularity) +
        CNn * (1 - cellularity))

    # Set BAF to NA for cases with CNn <= 1
    baf[CNn <= 1] <- NA

    return(baf)
}

#' Calculate theoretical mutation frequency
#'
#' @param CNt Copy number in tumor
#' @param Mt Number of mutant alleles
#' @param cellularity Tumor purity (0-1)
#' @param CNn Normal copy number (default: 2)
#' @return Vector of theoretical mutation frequencies
theoretical.mufreq <- function(CNt, Mt, cellularity, CNn = 2) {
    # Validate input parameters
    if (cellularity < 0 || cellularity > 1) {
        stop("Cellularity must be between 0 and 1")
    }

    # Calculate number of normal alleles
    normal_alleles <- (CNt - Mt) * cellularity + CNn * (1 - cellularity)

    # Calculate total number of alleles
    all_alleles <- (CNt * cellularity) + CNn * (1 - cellularity)

    # Return mutation frequency
    1 - (normal_alleles/all_alleles)
}

#' Generate matrix of possible B-allele frequency types
#'
#' @param CNt.min Minimum copy number in tumor
#' @param CNt.max Maximum copy number in tumor
#' @param CNn Normal copy number (default: 2)
#' @return Data frame with columns CNn, CNt, and B
baf.types.matrix <- function(CNt.min, CNt.max, CNn = 2) {
    # Validate input parameters
    if (CNt.min > CNt.max) {
        stop("CNt.min must be less than or equal to CNt.max")
    }

    # Calculate copy number ratios
    cn_ratio_vect <- seq(from = CNt.min/CNn, to = CNt.max/CNn,
        by = 1/CNn)
    CNt <- cn_ratio_vect * CNn

    # Generate B allele combinations
    if (CNn < 2) {
        # For CNn < 2, only consider 0 B alleles
        b_comb <- lapply(CNt, function(x) 0)
    } else {
        # Generate possible B allele counts for each CNt
        b_comb <- lapply(CNt, function(x) {
            seq(from = 0, to = trunc(x/2))
        })
    }

    # Create final data frame
    times_b <- sapply(b_comb, length)
    data.frame(CNn = CNn, CNt = rep(CNt, times = times_b), B = unlist(b_comb))
}

#' Generate matrix of possible mutation frequency types
#'
#' @param CNt.min Minimum copy number in tumor
#' @param CNt.max Maximum copy number in tumor
#' @param CNn Normal copy number (default: 2)
#' @return Data frame with columns CNn, CNt, and Mt
mufreq.types.matrix <- function(CNt.min, CNt.max, CNn = 2) {
    # Validate input parameters
    if (CNt.min > CNt.max) {
        stop("CNt.min must be less than or equal to CNt.max")
    }

    # Calculate copy number ratios
    cn_ratio_vect <- seq(from = CNt.min/CNn, to = CNt.max/CNn,
        by = 1/CNn)
    CNt <- cn_ratio_vect * CNn

    # Generate possible mutation counts for each CNt
    mut_comb <- lapply(CNt, function(x) seq(from = 0, to = x))

    # Create final data frame
    times_muts <- sapply(mut_comb, length)
    data.frame(CNn = CNn, CNt = rep(CNt, times = times_muts),
        Mt = unlist(mut_comb))
}

#' Calculate BAF model points
#'
#' @param cellularity Tumor purity (0-1)
#' @param ploidy Tumor ploidy
#' @param baf_types Data frame from baf.types.matrix()
#' @param avg.depth.ratio Average depth ratio
#' @return Data frame with BAF and depth ratio values
baf.model.points <- function(cellularity, ploidy, baf_types,
    avg.depth.ratio) {
    # Calculate theoretical depth ratios
    depth_ratio <- theoretical.depth.ratio(cellularity = cellularity,
        ploidy = ploidy, CNn = baf_types[, "CNn"], CNt = baf_types[,
            "CNt"], avg.depth.ratio = avg.depth.ratio)

    # Calculate theoretical BAF values
    baf <- theoretical.baf(cellularity = cellularity, CNn = baf_types[,
        "CNn"], CNt = baf_types[, "CNt"], B = baf_types[, "B"])

    # Return results as data frame
    data.frame(BAF = baf, depth.ratio = depth_ratio)
}

#' Calculate mutation frequency model points
#'
#' @param cellularity Tumor purity (0-1)
#' @param ploidy Tumor ploidy
#' @param mufreq_types Data frame from mufreq.types.matrix()
#' @param avg.depth.ratio Average depth ratio
#' @return Data frame with mutation frequencies and depth ratio values
mufreq.model.points <- function(cellularity, ploidy, mufreq_types,
    avg.depth.ratio) {
    # Calculate theoretical mutation frequencies
    mufreqs <- theoretical.mufreq(cellularity = cellularity,
        CNn = mufreq_types[, "CNn"], CNt = mufreq_types[, "CNt"],
        Mt = mufreq_types[, "Mt"])

    # Calculate theoretical depth ratios
    depth_ratio <- theoretical.depth.ratio(cellularity = cellularity,
        ploidy = ploidy, CNn = mufreq_types[, "CNn"], CNt = mufreq_types[,
            "CNt"], avg.depth.ratio = avg.depth.ratio)

    # Return results as data frame
    data.frame(mufreqs = mufreqs, depth.ratio = depth_ratio)
}


#' Function to calculate B-allele frequency with confidence intervals
#'
#' @param Af A-allele frequencies
#' @param Bf B-allele frequencies
#' @param good.readsNumber of good quality reads
#' @param conf Confidence level (default: 0.95)
#' @return Vector of length 3: (lower CI, B-allele frequency, upper CI)
b_allele_freq <- function(Af, Bf, good.reads, conf = 0.95) {
    if (length(Bf) > 1) {
        # Calculate weighted density of allele frequencies
        weights <- good.reads/(2 * sum(good.reads))
        dd <- density(x = c(Bf, Af), weights = c(weights, weights))

        # Find local maxima in density
        density_diff <- diff(dd$y)
        sign_changes <- diff(sign(density_diff))

        # Try to find strong local maxima first
        points.max <- which(sign_changes == -2) + 1

        # If no strong maxima found or global maximum not
        # among them, look for weaker local maxima
        if (length(points.max) < 1 || !which.max(dd$y) %in% points.max) {
            points.max <- which(sign_changes == -1) + 1
        }

        # Extract x and y coordinates of maxima
        l.max <- dd$x[points.max]
        d.max <- dd$y[points.max]

        # Find the B-allele value at the highest density
        # peak
        peak_indices <- which(dd$x %in% l.max)
        highest_peak <- which.max(dd$y[peak_indices])
        b.val <- l.max[highest_peak]

        # Handle case where no peak is found
        if (length(b.val) < 1) {
            message("WARNING: No clear peak found in density estimation")
            message("l.max:", paste(l.max, collapse = ", "))
            message("d.max:", paste(d.max, collapse = ", "))
            b.val <- min(l.max)
        }

        # Get density value at the chosen peak
        d.val <- d.max[which(l.max == b.val)]

        # Calculate confidence interval bounds
        density_threshold <- d.val - (d.val * (1 - conf))
        b.range <- range(dd$x[dd$y >= density_threshold])

        # Adjust B-allele frequency if it's greater than
        # 0.5
        if (b.val > 0.5) {
            b.val <- 1 - b.val
        }

        # Calculate symmetric confidence interval
        max_diff <- max(b.range) - b.val
        min_diff <- b.val - min(b.range)
        interval_width <- min(max_diff, min_diff)

        # Return lower bound, B-allele frequency, and upper
        # bound Lower CI
        return(c(b.val - interval_width, b.val, b.val + interval_width))
    } else if (length(Bf) == 1) {
        # If only one B-allele frequency, return it with no
        # uncertainty
        return(c(Bf, Bf, Bf))
    } else {
        # If no B-allele frequencies, return NA
        return(c(NA, NA, NA))
    }
}
