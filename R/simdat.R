## #################################################################
## Simulation used in the manuscript
## Summarized from Simulation\ data\ generator.R and ndata2.R
## #################################################################

library(survival)

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
#' W.1 -- W.4 are the corresponding censoring weights

datGen <- function(n, cp, sce, qt) {
    k <- 2
    X <- model.matrix(~ sample(0:1, n, TRUE))
    if (sce == 1) b <- matrix(c(1.609438, 1.410748, 1.219403, 1.040613 , 0, 0, 0, 0), 4)
    else b <- matrix(c(1.609438, 1.410748, 1.219403, 1.040613,
                       0.6931472, 0.7974189, 0.9070615, 1.0174711), 4)
    u <- runif(n)
    if (cp == 0) cen <- Inf
    if (sce == 2 & qt == .25) {
        if (cp == .1) cen <- t(matrix(runif(4 * n, 0, c(123.69, 115.24, 106.73, 101.62)), 4))
        if (cp == .3) cen <- t(matrix(runif(4 * n, 0, c(41.27, 41.27, 39.17, 37.34, 35.91)), 4))
        if (cp == .5) cen <- t(matrix(runif(4 * n, 0, c(23.55, 22.5, 21.61, 21.19)), 4))
        if (cp == .7) cen <- t(matrix(runif(4 * n, 0, c(14.2, 13.55, 13.07, 13.04)), 4))
    }
    if (sce == 2 & qt == .5) {
        if (cp == .1) cen <- t(matrix(runif(4 * n, 0, c(78.11, 70.39, 64.86, 61.17)), 4))
        if (cp == .3) cen <- t(matrix(runif(4 * n, 0, c(26.36, 24.35, 23.34, 22.53)), 4))
        if (cp == .5) cen <- t(matrix(runif(4 * n, 0, c(15.08, 14.07, 13.62, 13.43)), 4))
        if (cp == .7) cen <- t(matrix(runif(4 * n, 0, c(9.09, 8.49, 8.36, 8.5)), 4))
    }
    if (sce == 2 & qt == .75) {
        if (cp == .1) cen <- t(matrix(runif(4 * n, 0, c(55.68, 48.34, 43.86, 41.23)), 4))
        if (cp == .3) cen <- t(matrix(runif(4 * n, 0, c(18.79, 16.92, 16.08, 15.80)), 4))
        if (cp == .5) cen <- t(matrix(runif(4 * n, 0, c(10.69, 9.86, 9.57, 9.80)), 4))
        if (cp == .7) cen <- t(matrix(runif(4 * n, 0, c(6.49, 5.96, 6.03, 6.44)), 4))
    }
    if (sce == 1 & qt == .25) {
        if (cp == .1) cen <- t(matrix(runif(4 * n, 0, c(81.63, 73.44, 66.69, 61.86)), 4))
        if (cp == .3) cen <- t(matrix(runif(4 * n, 0, c(27.42, 25.32, 23.86, 22.66)), 4))
        if (cp == .5) cen <- t(matrix(runif(4 * n, 0, c(16.30, 15.34, 14.72, 14.33)), 4))
        if (cp == .7) cen <- t(matrix(runif(4 * n, 0, c(10.48, 9.92, 9.64, 9.58)), 4))
    }
    if (sce == 1 & qt == .50) {
        if (cp == .1) cen <- t(matrix(runif(4 * n, 0, c(52.15, 44.78, 39.03, 35)), 4))
        if (cp == .3) cen <- t(matrix(runif(4 * n, 0, c(17.76, 15.85, 14.55, 13.78)), 4))
        if (cp == .5) cen <- t(matrix(runif(4 * n, 0, c(10.50, 9.63, 9.22, 9.06)), 4))
        if (cp == .7) cen <- t(matrix(runif(4 * n, 0, c(6.8, 6.26, 6.2, 6.36)), 4))
    }
    if (sce == 1 & qt == .75) {
        if (cp == .1) cen <- t(matrix(runif(4 * n, 0, c(36.98, 30.1, 25.65, 21.92)), 4))
        if (cp == .3) cen <- t(matrix(runif(4 * n, 0, c(12.53, 10.82, 9.89, 9.46)), 4))
        if (cp == .5) cen <- t(matrix(runif(4 * n, 0, c(7.47, 6.68, 6.44, 6.59)), 4))
        if (cp == .7) cen <- t(matrix(runif(4 * n, 0, c(4.81, 4.61, 4.50, 4.41)), 4))
    }    
    ## Y <- (-log(u))^(1 / k) / ((-log(qt))^(1 / k) / exp(colSums(b[1,] * apply(X, 1, cumsum))))
    Y <- (-log(u))^(1 / k) / ((-log(qt))^(1 / k) / exp(colSums(b[1,] * t(X))))
    Z <- matrix(pmin(Y, cen), n)
    delta <- 1 * (Y <= cen)
    ## Calculate censoring weights
    W <- sapply(1:4, function(x) {
        sv <- survfit(Surv(Z[,x], 1 - delta[,x]) ~ 1)
        w <- delta[,x] / sv$surv[findInterval(Z[,x], sv$time)]
        w[is.na(w)] <- max(w, na.rm = TRUE)
        return(w)
    })
    delta[t(t(Z) <= 0:3)] <- NA
    W[t(t(Z) <= 0:3)] <- NA
    Z[t(t(Z) <= 0:3)] <- NA
    colnames(X) <- 1:2
    data.frame(Z = Z, X = X, d = delta, W = W)
}

summary(datGen(1e5, .1, 1, .5)[,7:10])
summary(datGen(1e5, .3, 1, .5)[,7:10])
summary(datGen(1e5, .5, 1, .5)[,7:10])
summary(datGen(1e5, .7, 1, .5)[,7:10])

summary(datGen(1e5, .1, 2, .5)[,7:10])
summary(datGen(1e5, .3, 2, .5)[,7:10])
summary(datGen(1e5, .5, 2, .5)[,7:10])
summary(datGen(1e5, .7, 2, .5)[,7:10])
