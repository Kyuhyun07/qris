## #################################################################
## Real data analysis
## #################################################################

library(survival)
library(nleqslv)
library(BB)
library(Rcpp)
library(RcppArmadillo)
library(pracma)

#' External variables; t0 = 1, 3, 5, 7
#' @param n is sample size
#' @param cp is censoring rate
#' @param sce scenario 1 assumes no treatment effect/ 2 = with treatment effect
#' @param qt is the quantile
#' @param t is the followup time
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
#'
#'
datGen <- function(n, cp, sce, qt, t) {
    k <- 2
    X <- model.matrix(~ sample(0:1, n, TRUE))
    if (sce == 1) {
      if (t == 0) {b <- c(1.609438, 0)}
      else if (t == 1) {b <- c(1.410748, 0)}
      else if (t == 2) {b <- c(1.219403, 0)}
      else {b <- c(1.040613, 0)}
    } else {
      if(t == 0) {b <- c(1.609438, 0.6931472)}
      else if (t == 1) {b <- c(1.410748, 0.7974189)}
      else if (t == 2) {b <- c(1.219403, 0.9070615)}
      else {b <- c(1.040613, 1.0174711)}
    }

    u <- runif(n)
    if (cp == 0) cen <- Inf
    if (sce == 2 & qt == .25) {
        if (cp == .1 & t == 0) cen <- t(runif(n, 0, 123.69))
        if (cp == .1 & t == 1) cen <- t(runif(n, 0, 115.24))
        if (cp == .1 & t == 2) cen <- t(runif(n, 0, 106.73))
        if (cp == .1 & t == 3) cen <- t(runif(n, 0, 101.62))
        if (cp == .3 & t == 0) cen <- t(runif(n, 0, 41.27))
        if (cp == .3 & t == 1) cen <- t(runif(n, 0, 39.17))
        if (cp == .3 & t == 2) cen <- t(runif(n, 0, 37.34))
        if (cp == .3 & t == 3) cen <- t(runif(n, 0, 35.91))
        if (cp == .5 & t == 0) cen <- t(runif(n, 0, 23.55))
        if (cp == .5 & t == 1) cen <- t(runif(n, 0, 22.5))
        if (cp == .5 & t == 2) cen <- t(runif(n, 0, 21.61))
        if (cp == .5 & t == 3) cen <- t(runif(n, 0, 21.19))
        if (cp == .7 & t == 0) cen <- t(runif(n, 0, 14.2))
        if (cp == .7 & t == 1) cen <- t(runif(n, 0, 13.55))
        if (cp == .7 & t == 2) cen <- t(runif(n, 0, 13.07))
        if (cp == .7 & t == 3) cen <- t(runif(n, 0, 13.04))
    }
    if (sce == 2 & qt == .5) {
        if (cp == .1 & t == 0) cen <- t(runif(n, 0, 78.11))
        if (cp == .1 & t == 1) cen <- t(runif(n, 0, 70.39))
        if (cp == .1 & t == 2) cen <- t(runif(n, 0, 64.86))
        if (cp == .1 & t == 3) cen <- t(runif(n, 0, 61.17))
        if (cp == .3 & t == 0) cen <- t(runif(n, 0, 26.36))
        if (cp == .3 & t == 1) cen <- t(runif(n, 0, 24.35))
        if (cp == .3 & t == 2) cen <- t(runif(n, 0, 23.34))
        if (cp == .3 & t == 3) cen <- t(runif(n, 0, 22.53))
        if (cp == .5 & t == 0) cen <- t(runif(n, 0, 15.08))
        if (cp == .5 & t == 1) cen <- t(runif(n, 0, 14.07))
        if (cp == .5 & t == 2) cen <- t(runif(n, 0, 13.62))
        if (cp == .5 & t == 3) cen <- t(runif(n, 0, 13.43))
        if (cp == .7 & t == 0) cen <- t(runif(n, 0, 9.09))
        if (cp == .7 & t == 1) cen <- t(runif(n, 0, 8.49))
        if (cp == .7 & t == 2) cen <- t(runif(n, 0, 8.36))
        if (cp == .7 & t == 3) cen <- t(runif(n, 0, 8.5))
    }
    if (sce == 2 & qt == .75) {
        if (cp == .1 & t == 0) cen <- t(runif(n, 0, 55.68))
        if (cp == .1 & t == 1) cen <- t(runif(n, 0, 48.34))
        if (cp == .1 & t == 2) cen <- t(runif(n, 0, 43.86))
        if (cp == .1 & t == 3) cen <- t(runif(n, 0, 41.23))
        if (cp == .3 & t == 0) cen <- t(runif(n, 0, 18.79))
        if (cp == .3 & t == 1) cen <- t(runif(n, 0, 16.92))
        if (cp == .3 & t == 2) cen <- t(runif(n, 0, 16.08))
        if (cp == .3 & t == 3) cen <- t(runif(n, 0, 15.80))
        if (cp == .5 & t == 0) cen <- t(runif(n, 0, 10.69))
        if (cp == .5 & t == 1) cen <- t(runif(n, 0, 9.86))
        if (cp == .5 & t == 2) cen <- t(runif(n, 0, 9.57))
        if (cp == .5 & t == 3) cen <- t(runif(n, 0, 9.80))
        if (cp == .7 & t == 0) cen <- t(runif(n, 0, 6.49))
        if (cp == .7 & t == 1) cen <- t(runif(n, 0, 5.96))
        if (cp == .7 & t == 2) cen <- t(runif(n, 0, 6.03))
        if (cp == .7 & t == 3) cen <- t(runif(n, 0, 6.44))
    }
    if (sce == 1 & qt == .25) {
        if (cp == .1 & t == 0) cen <- t(runif(n, 0, 81.636))
        if (cp == .1 & t == 1) cen <- t(runif(n, 0, 73.44))
        if (cp == .1 & t == 2) cen <- t(runif(n, 0, 66.69))
        if (cp == .1 & t == 3) cen <- t(runif(n, 0, 61.86))
        if (cp == .3 & t == 0) cen <- t(runif(n, 0, 27.42))
        if (cp == .3 & t == 1) cen <- t(runif(n, 0, 25.32))
        if (cp == .3 & t == 2) cen <- t(runif(n, 0, 23.86))
        if (cp == .3 & t == 3) cen <- t(runif(n, 0, 22.66))
        if (cp == .5 & t == 0) cen <- t(runif(n, 0, 16.30))
        if (cp == .5 & t == 1) cen <- t(runif(n, 0, 15.34))
        if (cp == .5 & t == 2) cen <- t(runif(n, 0, 14.72))
        if (cp == .5 & t == 3) cen <- t(runif(n, 0, 14.33))
        if (cp == .7 & t == 0) cen <- t(runif(n, 0, 10.48))
        if (cp == .7 & t == 1) cen <- t(runif(n, 0, 9.92))
        if (cp == .7 & t == 2) cen <- t(runif(n, 0, 9.64))
        if (cp == .7 & t == 3) cen <- t(runif(n, 0, 9.58))
    }
    if (sce == 1 & qt == .50) {
        if (cp == .1 & t == 0) cen <- t(runif(n, 0, 52.15))
        if (cp == .1 & t == 1) cen <- t(runif(n, 0, 44.78))
        if (cp == .1 & t == 2) cen <- t(runif(n, 0, 39.03))
        if (cp == .1 & t == 3) cen <- t(runif(n, 0, 35))
        if (cp == .3 & t == 0) cen <- t(runif(n, 0, 17.76))
        if (cp == .3 & t == 1) cen <- t(runif(n, 0, 15.85))
        if (cp == .3 & t == 2) cen <- t(runif(n, 0, 14.55))
        if (cp == .3 & t == 3) cen <- t(runif(n, 0, 13.78))
        if (cp == .5 & t == 0) cen <- t(runif(n, 0, 10.50))
        if (cp == .5 & t == 1) cen <- t(runif(n, 0, 9.63))
        if (cp == .5 & t == 2) cen <- t(runif(n, 0, 9.22))
        if (cp == .5 & t == 3) cen <- t(runif(n, 0, 9.06))
        if (cp == .7 & t == 0) cen <- t(runif(n, 0, 6.8))
        if (cp == .7 & t == 1) cen <- t(runif(n, 0, 6.26))
        if (cp == .7 & t == 2) cen <- t(runif(n, 0, 6.2))
        if (cp == .7 & t == 3) cen <- t(runif(n, 0, 6.36))
    }
    if (sce == 1 & qt == .75) {
        if (cp == .1 & t == 0) cen <- t(runif(n, 0, 36.98))
        if (cp == .1 & t == 1) cen <- t(runif(n, 0, 30.1))
        if (cp == .1 & t == 2) cen <- t(runif(n, 0, 25.65))
        if (cp == .1 & t == 3) cen <- t(runif(n, 0, 21.92))
        if (cp == .3 & t == 0) cen <- t(runif(n, 0, 12.53))
        if (cp == .3 & t == 1) cen <- t(runif(n, 0, 10.82))
        if (cp == .3 & t == 2) cen <- t(runif(n, 0, 9.89))
        if (cp == .3 & t == 3) cen <- t(runif(n, 0, 9.46))
        if (cp == .5 & t == 0) cen <- t(runif(n, 0, 7.47))
        if (cp == .5 & t == 1) cen <- t(runif(n, 0, 6.68))
        if (cp == .5 & t == 2) cen <- t(runif(n, 0, 6.44))
        if (cp == .5 & t == 3) cen <- t(runif(n, 0, 6.59))
        if (cp == .7 & t == 0) cen <- t(runif(n, 0, 4.81))
        if (cp == .7 & t == 1) cen <- t(runif(n, 0, 4.61))
        if (cp == .7 & t == 2) cen <- t(runif(n, 0, 4.50))
        if (cp == .7 & t == 3) cen <- t(runif(n, 0, 4.41))
    }
    ## Y <- (-log(u))^(1 / k) / ((-log(qt))^(1 / k) / exp(colSums(b[1,] * apply(X, 1, cumsum))))
    Y <- (-log(u))^(1 / k) / ((-log(qt))^(1 / k) / exp(colSums(b * t(X))))
    Z <- matrix(pmin(Y, cen), n)
    delta <- 1 * (Y <= cen)
    # ## Calculate censoring weights
    # sv <- survfit(Surv(Z, t(1 - delta)) ~ 1)
    # ## Changed to Li's weight (210415) by Q
    # W <- delta / sv$surv[findInterval(Z, sv$time)]*sv$surv[min(which(floor(sv$time)==(t)))]
    # W[is.na(W)] <- max(W, na.rm = TRUE)
    # delta = t(delta)
    # W = t(W)
    # delta[t(t(Z) <= t)] <- NA
    # W[t(t(Z) <= t)] <- NA
    # Z[t(t(Z) <= t)] <- NA
    colnames(X) <- 1:2
    data = data.frame(Z = Z, X = X, d = t(delta))
}

sourceCpp(code = '
    #include <RcppArmadillo.h>
    // [[Rcpp::depends(RcppArmadillo)]]
    using namespace arma;
    // [[Rcpp::export]]
    arma::mat isObj(arma::vec b, arma::mat X, arma::vec W, arma::mat H, arma::vec I,
                    arma::vec logT, double Q) {
    arma::mat m1 = X % repmat(I, 1, X.n_cols);
    arma::mat m2 = normcdf((X * b - logT) / sqrt(diagvec(X * H * X.t()))) % W - Q;
    return m1.t() * m2;
  }')

#' #' Induce smoothing objective equation, adopted from is_optim_objectF()
#' sourceCpp(code = '
#'   #include <RcppArmadillo.h>
#'   // [[Rcpp::depends(RcppArmadillo)]]
#'   using namespace arma;
#'   // [[Rcpp::export]]
#'   arma::mat isOpm(arma::vec b, arma::mat X, arma::vec W, arma::mat H,
#'                   arma::vec Z, double t0, double Q) {
#'   arma::mat se = sqrt(diagvec(X * H * X.t()));
#'   arma::mat xdif = X * b - log(Z - t0);
#'   arma::mat m1 = W;
#'   arma::mat m2 = xdif % (normcdf(xdif / se) - Q);
#'   arma::mat m3 = normpdf(xdif / se) % se;
#'   return m1.t() * (m2 + m3);
#' }')
#'
#' #' Non smoothed estimating equation, adopted from rq_objectF
#' sourceCpp(code = '
#'   #include <RcppArmadillo.h>
#'   // [[Rcpp::depends(RcppArmadillo)]]
#'   using namespace arma;
#'   // [[Rcpp::export]]
#'   arma::mat rqObj(arma::vec b, arma::mat X, arma::vec W,
#'                   arma::vec Z, double t0, double Q) {
#'   arma::mat m1 = X % repmat(1, X.n_cols);
#'   arma::mat m2 = Q - arma::conv_to<arma::mat>::from((log(Z - t0) <= X * b)*W);
#'   return m1.t() * m2;
#' }')
#'
#' #' Non smoothed estimating equation, adopted from rq_optim_objectF
#' sourceCpp(code = '
#'   #include <RcppArmadillo.h>
#'   // [[Rcpp::depends(RcppArmadillo)]]
#'   using namespace arma;
#'   // [[Rcpp::export]]
#'   arma::mat rqOpm(arma::vec b, arma::mat X, arma::vec W,
#'                   arma::vec Z, double t0, double Q) {
#'   arma::mat xdif = X * b - log(Z - t0);
#'   arma::mat m1 = W % xdif;
#'   arma::mat m2 = arma::conv_to<arma::mat>::from((log(Z - t0) <= X * b)) - Q;
#'   return m1.t() * m2;
#' }')


# Assumes data is created by `datGen`
# Example data generation
n = 200
cp = .1
sce = 2
qt = .5
t = 1
dat <- datGen(n, cp, sce, qt, t)

#' @param z a vector of z = min(T,C)
#' @param x a matrix of covariate (n x p) (n = number of subject, p = number of covariate)
#' @param w weight
#' @param b0 initial guess of regression parameter (length should be p)
#' @param t followup time

# Test parameter
z <- dat[,1]
x <- as.matrix(dat[,3])
d <- dat[,4]
b0 <- rep(0.5, p+1)

ismb <- function() {z, x, d, b0, t}
  # Li's weight
  sv <- survfit(Surv(z, t(1 - t(d))) ~ 1)
  w <- t(d) / sv$surv[findInterval(z, sv$time)]*sv$surv[min(which(floor(sv$time)==(t)))]
  w[is.na(w)] <- max(w, na.rm = TRUE)
  w = t(w)
  d[t(t(z) <= t)] <- NA
  w[t(t(z) <= t)] <- NA
  z[t(t(z) <= t)] <- NA
  # reorganize data
  dat0 = as.matrix(cbind(z,rep(1,nrow(x)),x,d,w))

  # number of subject
  n <- length(z)
  # number of covariate (exclue intercept)
  p <- ncol(x)
  # H matrix
  H <- diag(p+1) / n
  logT <- log(dat0[,1]-t)
  W <- dat0[,5]
  X <- dat0[,2:3]
  I <- as.numeric(dat0[,1]>=t)

  f1 <-nleqslv(b0, function(b) isObj(b, X, W, H, I, logT, qt))
  # f2 <- optim(b0, function(b) isOpm(b, X, W, H, Z, 0, qt))
  # f3 <- nleqslv(b0, function(b) rqObj(b, X, W, Z, 0, qt))
  # f4 <- optim(b0, function(b) rqOpm(b, X, W, Z, 0, qt))
  c(f1$x, f1$termcd) # f2$par, f3$x, f4$par, f1$termcd, f2$convergence, f3$termcd, f4$convergence)
}

do(200, .1, 1, .5, 1)
do(1000, .1, 1, .5)
