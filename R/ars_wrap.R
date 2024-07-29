f_y <- function(y, K = 3, U = 10) {
    y * (K - 3/2) - 1/(2 * exp(y)) +
        lgamma(exp(y)) -
        lgamma(
            exp(y) +
                U
        )
}

f_y_prima <- function(y, K = 3, U = 10) {
    (K - 3/2) + 1/(2 * exp(y)) +
        exp(y) *
            digamma(exp(y)) -
        exp(y) *
            digamma(
                exp(y) +
                  U
            )
}

sample_alpha <- function(K = 3, U = 10) {
    gen <- ars.new(f_y, f_y_prima, -Inf, Inf, K = K, U = U)
    exp(ur(gen, 1))
}

beta.f_y <- function(y, S, W) {

    F <- dim(W)[1]
    R <- dim(S)[3]

    sum_SW <- 0
    for (r in 1:R) {
        Sr <- S[, , r]
        Sr <- as.matrix(Sr)  # compatibility with one-dimensional case
        sum_SW <- sum_SW + log(
            det(Sr) *
                det(W)
        ) -
            sum(diag(Sr %*% W))
    }

    sum_gamma <- 0
    for (d in 1:F) {
        sum_gamma <- sum_gamma + lgamma(
            (exp(y) +
                d - F)/2
        )
    }
    r <- y
    r <- r - R * sum_gamma
    r <- r - F/(2 * (exp(y) -
        F + 1))
    r <- r - (3/2) * log(
        exp(y) -
            F + 1
    )
    r <- r + ((R * exp(y) *
        F)/2) * (y - log(2))
    r <- r + (exp(y)/2) *
        sum_SW

    return(r)
}

beta.f_y_prima <- function(y, S, W) {

    F <- dim(W)[1]
    R <- dim(S)[3]

    sum_SW <- 0
    for (r in 1:R) {
        Sr <- S[, , r]
        Sr <- as.matrix(Sr)  # compatibility with one-dimensional case
        sum_SW <- sum_SW + log(
            det(Sr) *
                det(W)
        ) -
            sum(diag(Sr %*% W))
    }

    sum_digamma <- 0
    for (d in 1:F) {
        sum_digamma <- sum_digamma + digamma(
            (exp(y) +
                d - F)/2
        )
    }
    r <- 1
    r <- r - R * exp(y) *
        0.5 * sum_digamma
    r <- r + exp(y) *
        F/(2 * ((exp(y) -
        F + 1))^2)
    r <- r - (3/2) * exp(y)/(exp(y) -
        F + 1) + (R * F * exp(y))/2 +
        R * F * exp(y) *
            (y - log(2))/2
    r <- r + (exp(y)/2) *
        sum_SW

    return(r)
}

sample_beta_a0 <- function(S, W) {
    # Sample the degrees of freedom of the Wishart given a
    # scale matrix W and some observed precision matrices S
    # make it compatible with the one-dimensional case in
    # case we want to use this version for both the multi
    # and the uni-dimensional
    if (is.matrix(W)) {
        xlb <- log(
            dim(W)[1] -
                0.9
        )
    }
    gen <- ars.new(
        beta.f_y, beta.f_y_prima, lb = xlb, ub = Inf, S = S,
        W = W
    )
    exp(ur(gen, 1))
}
