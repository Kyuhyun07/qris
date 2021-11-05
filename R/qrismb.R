#' Estimate a quantile regression estimator of residual lifetime from survival data
#'
#' Using two estimation methods
#' 1. L1-minimization(non-smooth estimating equation)
#' 2. Induced smoothing approach (smooth estimating equation)
#'
#' @param Z is a vector of observed time, which is minimum of failure time and censored time
#' @param nc is a number of covariates used in analysis
#' @param covariate is a matrix of covariate (# row = # of subject, # of column = # of covariate(nc))
#' @param D is a vector of censoring indicator (1 = not censored, 0 = censored)
#' @param t_0 is the followup time(or basetime of analysis)
#' @param Q is the quantile
#' @param ne is number of multiplier bootstrapping for V matrix estimation
#' @param init is option for initial guess of regression parameter ("random" assumes all coefficients as random numbers, "one" assumes all coefficients as 1s, otherwise a solution from rq function)
#' @param method is option how to estimate coefficient and standard error of it
#' ("nonsmooth" uses non-smooth estimating equation : L1-minimization method in coefficient estimation, and full multiplier bootstrap in standard error estimation.
#' "smooth" uses induced smoothed estimating equation : nonlinear equation solver in coefficient estimation and partial multiplier bootstrap in standard error estimation
#' "iterative" uses induced smoothed estmating eqution and iterative calculation updating coefficient and SE).
#' @return An object of class "\code{qrismb}" representing the fit.
#' The \code{qrismb} object is a list containing at least the following components:
#' \describe{
#'   \item{coefficient}{a vector of point estimates}
#'   \item{stderr}{a vector of standard error of point estiamtes}
#'   \item{iterno}{a number of itertation until convergence (only for iterative procedure)}
#'   }
#'
#' @export
#' @importFrom survival Surv
#' @importFrom survival survfit
#' @importFrom quantreg rq.wfit
#' @importFrom nleqslv nleqslv
#' @example inst/examples/ex_qrismb.R
qrismb <- function(Z, covariate, D, t_0 = 0, Q = 0.5, ne = 100, init="rq", method="smooth"){
    nc <- ncol(covariate)
    if(nc < 1){
        stop("Use at least one covariate")
    }
    if(sum(D!=0&D!=1)>=1) {
        stop("Delta must consist of only 0 and 1")
    }
    if(t_0<0) {
        stop("basetime must be 0 and positive number")
    }
    if(length(Q) > 1) {
        stop("Multiple taus not allowed in qrismb")
    }
    if(Q<=0|Q>=1) {
        stop("Tau must be scalar number between 0 and 1")
    }
    if(ne<=1) {
        stop("number of multiplier bootstrapping must greater than 1")
    }

    n = length(Z)
    data <- matrix(NA, n, nc+4)
    data[,1] <- Z
    ## Suppress warning message
    data[,2] <-suppressWarnings(log(Z-t_0))
    data[,3] <- as.numeric(Z>=t_0)
    data[,4:(nc+3)] <- as.matrix(covariate)
    data[is.na(data[,2]),2] <- -10
    data[data[,2]==-Inf,2] <- -10
    data[,(nc+4)] <- D
    data[n,(nc+4)] <- 1
    data <- as.data.frame(data) #head(data,20)
    colnames(data)[1:3] <- c("Z", "log(Z-t_0)", "I[Z>t_0]")
    covar <- paste("covariate",1:nc,sep = "")
    colnames(data)[4:(nc+3)] <- covar
    colnames(data)[(nc+4)] <- c("delta")
    # colnames(data)[(nc+4):(nc+5)] <- c("delta","Weight")
    data = na.omit(data)
    # data <- data[ complete.cases(data[ , c(4:3+nc)]), ]
    n = nrow(data)

    ## Rcpp IPCW with jump weight
    sv <- survfit(Surv(data[,1], 1 - data[,(nc+4)]) ~ 1)
    if (t_0<=sv$time[1]){
        ghatt_0 = 1
    } else {
        ghatt_0 = sv$surv[max(which(floor(sv$time)<(t_0)))]
    }
    W <- data[,(nc+4)] / sv$surv[findInterval(data[,1], sv$time)]*ghatt_0
    W[is.na(W)] <- max(W, na.rm = TRUE)
    data[,(nc+5)] <- W
    colnames(data)[(nc+5)] <- c("weight")

    # Covariate setting
    X <- as.matrix(cbind(c(rep(1,n)),data[,4:(nc+3)]))
    logT <- data[,2]
    I <- data[,3]
    H <- diag(1/n, nc+1, nc+1)

    ## Guess beta option (rq : solution from rq, 0 : no effect of covariate, others : covariate)
    if (init == "random"){
        betastart <- rnorm(nc+1)
    } else if (init == "one"){
        betastart <- c(1,rep(0,nc))
    } else {
        betastart <- as.vector(rq.wfit(X, data[,2], tau=Q, weight=W)$coef)
    }

    if (method == "nonsmooth"){
        ## Estimating equation for estimaing beta (using rq)
        pseudo1 <- -apply(X*I*W,2,sum)
        pseudo2 <- 2*apply(X*I*Q,2,sum)
        M <- 10^6
        Y.reg <- c(data[,2],M,M)
        X.reg <- rbind(X,rbind(pseudo1,pseudo2))
        W.reg <- c(I*W,1,1)
        Li.fit <- rq.wfit(X.reg,Y.reg,weights=W.reg)
        coefficient <- as.vector(Li.fit$coefficients)
        ## Full multiplier bootstrap
        fb_result <- c()
        ## if (Li.fit$code == 1 | Li.fit$code == 2){
        if (all(Li.fit$coefficients<=10)){
            for (j in 1:ne){
                # generating perturbation variable
                eta <- rexp(n,1)
                if (all(data$delta==rep(1,n))){
                    W_star <- rep(1,n)
                } else {
                    Gest <- ghat(data$Z,1-data$delta,eta)
                    W_star <- data$delta / Gest$survp[findInterval(data$Z, Gest$deathtime)]*Gest$survp[max(which(floor(Gest$deathtime)<(t_0)))]
                    W_star[is.na(W_star)] <- max(W_star, na.rm = TRUE)
                }
                rev_pseudo1 <- -apply(X*I*W_star*eta,2,sum)
                rev_pseudo2 <- 2*apply(X*I*eta*Q,2,sum)
                M <- 10^6
                rev_Y.reg <- c(data[,2],M,M)
                rev_X.reg <- rbind(X,rbind(rev_pseudo1,rev_pseudo2))
                rev_W.reg <- c(W_star*I*eta,1,1)
                full_boot <- rq.wfit(rev_X.reg,rev_Y.reg,weights=rev_W.reg)
                if (all(full_boot$coef<=2)){
                    fb_result <- cbind(fb_result,as.vector(full_boot$coef))
                } else {
                    fb_result <- cbind(fb_result,rep(NA,nc+1))
                }
            }
            sigma <- try(cov(t(fb_result), use="complete.obs"), silent = T)
            if(class(sigma)[1] == "try-error"|sum(!is.na(apply(fb_result,2,sum)))<=1){
                stop("More resampling iterations are necessary")
            } else {
                se <- sqrt(diag(sigma))
                list(coefficient = coefficient, stderr = se)
            }
        } else {
            coefficient <- c(NA,rep(NA,nc))
            se <- c(NA,rep(NA,nc))
            list(coefficient=coefficient, stderr = se)
        }
    } else if (method == "iterative") {
          iter_beta_result <- old_beta <- new_beta <- betastart
          new_sigma <- old_sigma <- H
          iter_SE_result <- sqrt(diag(old_sigma))

          for (k in 1:10){
              number_iter <- k
              old_beta <- new_beta
              old_sigma <- new_sigma
              slope_a <- Amat(old_beta, X, W, old_sigma, I, logT, Q)
              # Step 1 : Update beta()
              new_beta <- old_beta - qr.solve(slope_a) %*% isObj(old_beta, X, W, H, I, logT, Q)
              iter_beta_result <- rbind(iter_beta_result, t(new_beta))
              # Step 2 : Update Sigma()
              result.ismb <- c()
              for (j in 1:ne){
                  # generating perturbation variable
                  eta = rexp(n,1)
                  if (D==rep(1,n)){
                      W_star <- rep(1,n)
                      } else {
                          Gest <- ghat(data$Z,1-data$delta, eta)
                          W_star <- data$delta / Gest$survp[findInterval(data$Z, Gest$deathtime)]*Gest$survp[max(which(floor(Gest$deathtime)<(t_0)))]
                          W_star[is.na(W_star)] <- max(W_star, na.rm = TRUE)
                      }
                  result <- rev_isObj(old_beta, X, W_star, H, eta, I, logT, Q)
                  result.ismb <- cbind(result.ismb,result)
              }
              new_V <- cov(t(result.ismb))
              new_sigma <- t(qr.solve(slope_a)) %*% new_V %*% qr.solve(slope_a)
              iter_SE_result <- rbind(iter_SE_result , sqrt(diag(new_sigma)))
              if(norm(new_beta-old_beta, "F")<1e-6) break
          }

          # Last iteration
          old_beta <- new_beta
          old_sigma <- new_sigma
          slope_a <- Amat(old_beta, X, W, old_sigma, I, logT, Q)
          new_beta <- old_beta - qr.solve(slope_a) %*% isObj(old_beta, X, W, H, I, logT, Q)
          iter_beta_result <- rbind(iter_beta_result, t(new_beta))
          result.ismb <- c()
          for (j in 1:ne){
              # generating perturbation variable
              eta <- rexp(n,1)
              if (D==rep(1,n)){
              W_star <- rep(1,n)
              } else {
                  Gest <- ghat(data$Z,1-data$delta,eta)
                  W_star <- data$delta / Gest$survp[findInterval(data$Z, Gest$deathtime)]*Gest$survp[max(which(floor(Gest$deathtime)<(t_0)))]
                  W_star[is.na(W_star)] <- max(W_star, na.rm = TRUE)
              }
              result <- rev_isObj(old_beta, X, W_star, H, eta, I, logT, Q)
              result.ismb = cbind(result.ismb,result)
          }
          new_V = cov(t(result.ismb))
          new_sigma = t(qr.solve(slope_a)) %*% new_V %*% qr.solve(slope_a)
          iter_SE_result = rbind(iter_SE_result , sqrt(diag(new_sigma)))
          list(coefficient=iter_beta_result, stderr = iter_SE_result, iterno = k+1)
        } else {
        rcpp.fit <- nleqslv(betastart, function(b) isObj(b, X, W, H, I, logT, Q))
        if (rcpp.fit$termcd == 1 | rcpp.fit$termcd == 2){
            coefficient <- rcpp.fit$x
            ## Variance estimation : ISMB
            rcpp.result.ismb <- c()
            for (j in 1:ne){
                ## generating perturbation variable
                E <- rexp(n,1)
                if (all(data$delta==rep(1,n))){
                    W_star <- rep(1,n)
                } else {
                    Gest <- ghat(data$Z,1-data$delta,E)
                    W_star <- data$delta / Gest$survp[findInterval(data$Z, Gest$deathtime)]*Gest$survp[max(which(floor(Gest$deathtime)<(t_0)))]
                    W_star[is.na(W_star)] <- max(W_star, na.rm = TRUE)
                }
                result <- rev_isObj(coefficient, X, W_star, H, E, I, logT, Q)
                rcpp.result.ismb <- cbind(rcpp.result.ismb,result)
            }
            v <- cov(t(rcpp.result.ismb))
            rcpp.a <- Amat(coefficient, X, W, H, I, logT, Q)
            inva <- try(solve(rcpp.a))
            if(class(inva)[1] == "try-error"){
                se <- rep(NA,nc+1)
                stop("Slope matrix is singular matrix")
            } else {
                sigma <- t(inva) %*% v %*% inva
                se <- sqrt(diag(sigma))
            }
            list(coefficient=coefficient, stderr = se)
        } else {
            coefficient <- c(NA,rep(NA,nc))
            se <- c(NA,rep(NA,nc))
            list(coefficient=coefficient, stderr = se)
        }
    }
}

#' Estimate Kaplan Meier estimate of the survival function of the censoring time C
#' ghat
#'
#' @param T is a vector of observed time, which is minimum of failure time and censored time
#' @param event is a vector of censoring indicator (not censored = 1, censored = 0)
#' @param wgt is a vector of weight
#'
#'@return An object of class "\code{ghat}" estimate KM estimator of survival function of
#' The \code{ghat} object is a list containing at least the following components:
#' \describe{
#'   \item{deathtime}{the observed time}
#'   \item{ndeath}{a vector of number of subject who experienced event at deathtime}
#'   \item{nrisk}{a vector of number of subject who are possible to experience event at deathtime}
#'   \item{survp}{a vector of survival probability at deathtime}
#'   }
#' @export
#'
ghat <- function(T,censor,wgt=1){
    deathtime <- c(0,unique(sort(T[censor[]==1])))
    nrisk <- ndeath <- survp <- rep(0,length(deathtime))
    nrisk[1] <- sum((0<=T)*wgt)
    ndeath[1] <- sum((T==0)*censor*wgt)
    for(i in 2:length(deathtime)){
        nrisk[i] <- sum((deathtime[i-1]<=T)*wgt)
        ndeath[i] <- sum((T==deathtime[i-1])*censor*wgt)
    }
    prodobj <- 1-ndeath/nrisk
    for(i in 1:length(deathtime)){survp[i] <- prod(prodobj[1:i])}
    return(data.frame(cbind(deathtime,ndeath,nrisk,survp)))}
