## #################################################################
## Simulation used in the manuscript
## Summarized from Simulation\ data\ generator.R and ndata2.R
## #################################################################

n <- 10

#' External variables; t0 = 1, 3, 5, 7
#' @param n is sample size
#' @param cp is censoring rate
#' @param sce scenario 1 assumes no treatment effect
#' @param q is the quantile
#'
#' Internal variables
#' @param b is true parameter matrix
#' for the simulation in the manuscript it is 4 by 2; 4 t0's 2 X's
#' @param k shape parameter of weibull
#'
#' beta_0 is calculated by
#'  t0 <- 0:3
#'  p0 <- (-log(q))^.5 / 5
#'  log((-log(q) + (p0 * t0)^2)^.5 / p0 - t0)
#'
#'
#' Return variables:
#' X.1 X.2 are covariates
#' Z.1 -- Z.4 are the residual times at t0 = 1, 2, 3, 4
#' d.1 -- d.4 are the censoring indicators for the corresponding Z's

n <- 100
cen <- .2
sce <- 2
q <- .25

datGen <- function(n, cp, sce, q) {
    k <- 2
    X <- model.matrix(~ sample(0:1, n, TRUE))
    ## Calculate b
    b0 <- log(5)
    t0 <- 0:3
    if (sce == 1) b <- matrix(c(1.609438, 1.410748, 1.219403, 1.040613 , 0, 0, 0, 0), 4)
    else b <- matrix(c(1.609438, 1.410748, 1.219403, 1.040613,
                       0.6931472, 0.7974189, 0.9070615, 1.0174711), 4)
    u <- runif(n)
    if (cp == Inf) cen <- Inf
    if (cp == .1) cen <- t(matrix(runif(4 * n, 0, c(123.69, 70.39, 64.86, 61.17)), 4))
    if (cp == .3) cen <- t(matrix(runif(4 * n, 0, c(41.27, 24.35, 23.34, 22.53)), 4))
    if (cp == .5) cen <- t(matrix(runif(4 * n, 0, c(23.55, 14.07, 13.62, 13.43)), 4))
    if (cp == .7) cen <- t(matrix(runif(4 * n, 0, c(14.2, 8.49, 8.36, 8.5)), 4))
    Y <- (-log(u))^(1 / k) / ((-log(q))^(1 / k) / exp(X %*% t(b)))
    Z <- pmin(Y, cen)
    delta <- 1 * (Y <= cen)
    colnames(X) <- 1:2
    data.frame(Z = Z, X = X, d = delta)
}

