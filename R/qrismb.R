#' Estimate a quantile regression estimator of residual lifetime from survival data
#'
#' Using two estimation methods
#' 1. L1-minimization(non-smooth estimating equation)
#' 2. Induced smoothing approach (smooth estimating equation)
#' @useDynLib qrismb
#' @importFrom Rcpp sourceCpp
NULL
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
#' (otherwise uses induced smoothed estimating equation : nonlinear equation solver in coefficient estimation and partial multiplier bootstrap in standard error estimation).
#' @return An object of class "\code{qrismb}" representing the fit.
#' The \code{qrismb} object is a list containing at least the following components:
#' \describe{
#'   \item{coefficient}{a vector of point estimates}
#'   \item{stderr}{a vector of standard error of point estiamtes}
#'   }
#' @examples
#' data("retinopathy")
#' reti = retinopathy
#' reti_rev = reti
#' reti_rev$risk = rescale(reti$risk, to=c(0,1), from=c(0,12))
#' reti_rev$type = as.numeric(reti$type)-1
#' reti_rev$typextrt = reti_rev$type*reti$trt
#' Z=reti_rev$futime
#' nc=5
#' covariate=as.matrix(reti_rev[,c(9,4,5,6,10)])
#' D=reti_rev$status
#' t_0=1
#' Q=0.25
#' ne=200
#' qrismb(Z, nc, covariate, D, t_0, Q, ne, "rq", "smooth")
#' qrismb(Z, nc, covariate, D, t_0, Q, ne, "random", "nonsmooth")
qrismb= function(Z, nc, covariate, D, t_0 = 0, Q = 0.5, ne = 100, init, method){
  if(nc < 1)
    stop("Use at least one covariate")
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
  data = matrix(NA, n, nc+5)
  data[,1] = Z
  # Suppress warning message
  options(warn=-1)
  data[,2] = log(Z-t_0)
  data[,3] = as.numeric(Z>=t_0)
  data[,4:(nc+3)] = as.matrix(covariate)
  data[is.na(data[,2]),2]=-10
  data[,(nc+4)] = D
  data[n,(nc+4)] = 1
  data = as.data.frame(data) #head(data,20)
  colnames(data)[1:3]=c("Z", "log(Z-t_0)", "I[Z>t_0]")
  covar = paste("covariate",1:nc,sep = "")
  colnames(data)[4:(nc+3)] = covar
  colnames(data)[(nc+4):(nc+5)] = c("delta","Weight")

  # Rcpp Li's weight with jump weight
  sv <- survfit(Surv(data[,1], 1 - data[,(nc+4)]) ~ 1)
  W <- data[,(nc+4)] / sv$surv[findInterval(data[,1], sv$time)]*sv$surv[min(which(floor(sv$time)==(t_0)))]
  W[is.na(W)] <- max(W, na.rm = TRUE)
  data[,(nc+5)] = W

  # Covariate setting (1 covariate)

  X = as.matrix(cbind(c(rep(1,n)),data[,4:(nc+3)]))
  W = data[,(nc+5)]
  logT = data[,2]
  I = data[,3]
  H = diag(1/n, nc+1, nc+1)

  # Guess beta option (rq : solution from rq, 0 : no effect of covariate, others : covariate)
  if (init == "random"){
    betastart = rnorm(nc+1)
  } else if (init == "one"){
    betastart = c(1,rep(0,nc))
  } else {
    betastart = as.vector(rq.wfit(X,data[,2], tau=Q, weight=W)$coef)
  }

  if (method == "nonsmooth"){
    # Estimating equation for estimaing beta (using rq)
    pseudo1 = -apply(X*I*W,2,sum)
    pseudo2 = 2*apply(X*I*Q,2,sum)
    M = 10^6
    Y.reg = c(data[,2],M,M)
    X.reg = rbind(X,rbind(pseudo1,pseudo2))
    W.reg = c(I*W,1,1)
    Li.fit = rq.wfit(X.reg,Y.reg,weights=W.reg)
    coefficient = as.vector(Li.fit$coefficients)

    # Full multiplier bootstrap
    fb_result = c()
    # if (Li.fit$code == 1 | Li.fit$code == 2){
    if (all(Li.fit$coefficients<=10)){
      # Variance estimation : Li's method (Full multiplier bootstrap)
      SC.func=function(T,censor,wgt=1){
        deathtime=unique(sort(T[censor[]==1]))
        nrisk=ndeath=rep(0,length(deathtime))
        for(i in 1:length(deathtime)){
          nrisk[i]=sum((deathtime[i]<=T)*wgt)
          ndeath[i]=sum((T==deathtime[i])*censor*wgt)
        }
        prodobj=1-ndeath/nrisk
        survp=rep(0,length(deathtime))
        for(i in 1:length(deathtime)){survp[i]=prod(prodobj[1:i])}
        return(data.frame(cbind(deathtime,ndeath,nrisk,survp)))}
      for (j in 1:ne){
        # generating perturbation variable
        eta = rexp(n,1)
        if (all(D==rep(1,n))){
          W_star = rep(1,n)
        } else {
          Gest=SC.func(Z,1-D,eta)
          SC=stepfun(Gest$deathtime,c(1,Gest$survp))
          W_star <- D / Gest$survp[findInterval(Z, Gest$deathtime)]*Gest$survp[min(which(floor(Gest$deathtime)==(t_0)))]
          W_star[is.na(W_star)] <- max(W_star, na.rm = TRUE)
        }
        rev_pseudo1 = -apply(X*I*W_star*eta,2,sum)
        rev_pseudo2 = 2*apply(X*I*eta*Q,2,sum)
        M = 10^6
        rev_Y.reg = c(data[,2],M,M)
        rev_X.reg = rbind(X,rbind(rev_pseudo1,rev_pseudo2))
        rev_W.reg = c(W_star*I*eta,1,1)
        full_boot = rq.wfit(rev_X.reg,rev_Y.reg,weights=rev_W.reg)
        if (all(full_boot$coef<=2)){
          fb_result = cbind(fb_result,as.vector(full_boot$coef))
        } else {
          fb_result = cbind(fb_result,rep(NA,nc+1))
        }
      }
      sigma = cov(t(fb_result), use="complete.obs")
      se = sqrt(diag(sigma))
      list(coefficient=coefficient, stderr = se)
    } else {
      coefficient = c(NA,rep(NA,nc))
      se = c(NA,rep(NA,nc))
      list(coefficient=coefficient, stderr = se)
    }
  } else {
    #' sourceCpp(code = '
    #'       #include <RcppArmadillo.h>
    #'       // [[Rcpp::depends(RcppArmadillo)]]
    #'       using namespace arma;
    #'       // [[Rcpp::export]]
    #'       arma::mat isObj(arma::vec b, arma::mat X, arma::vec W, arma::mat H, arma::vec I,
    #'       arma::vec logT, double Q) {
    #'       arma::mat m1 = X % repmat(I, 1, X.n_cols);
    #'       arma::mat m2 = normcdf((X * b - logT) / sqrt(diagvec(X * H * X.t()))) % W - Q;
    #'       return m1.t() * m2;
    #'       }')
    #'
    #'
    #' #' Induce smoothing estimating equation, adopted from rev_is_objectF()
    #' sourceCpp(code = '
    #'       #include <RcppArmadillo.h>
    #'       // [[Rcpp::depends(RcppArmadillo)]]
    #'       using namespace arma;
    #'       // [[Rcpp::export]]
    #'       arma::mat rev_isObj(arma::vec b, arma::mat X, arma::vec W, arma::mat H, arma::vec E,arma::vec I,
    #'       arma::vec logT, double Q) {
    #'       arma::mat m1 = X % repmat(I, 1, X.n_cols) % repmat(E, 1, X.n_cols);
    #'       arma::mat m2 = normcdf((X * b - logT) / sqrt(diagvec(X * H * X.t()))) % W - Q;
    #'       return m1.t() * m2;
    #'       }')
    #'
    #' #' Induce smoothing estimating equation, adopted from Amat
    #' sourceCpp(code = '
    #'       #include <RcppArmadillo.h>
    #'       // [[Rcpp::depends(RcppArmadillo)]]
    #'       using namespace arma;
    #'       // [[Rcpp::export]]
    #'       arma::mat Amat(arma::vec b, arma::mat X, arma::vec W_star, arma::mat H, arma::vec E, arma::vec I,
    #'       arma::vec logT, double Q) {
    #'       arma::mat m1 = X % repmat(I, 1, X.n_cols) % repmat(W_star, 1, X.n_cols);
    #'       arma::mat m2 =  (-X / repmat(sqrt(diagvec(X * H * X.t())), 1, X.n_cols)) % repmat(normpdf((X * b - logT) / sqrt(diagvec(X * H * X.t()))), 1, X.n_cols);
    #'       return m1.t() * m2;
    #'       }')
    rcpp.fit <- nleqslv(betastart, function(b) isObj(b, X, W, H, I, logT, Q))

    if (rcpp.fit$termcd == 1 | rcpp.fit$termcd == 2){
      coefficient = rcpp.fit$x
      # Variance estimation : ISMB
      rcpp.result.ismb=c()

      SC.func=function(T,censor,wgt=1){
        deathtime=unique(sort(T[censor[]==1]))
        nrisk=ndeath=rep(0,length(deathtime))
        for(i in 1:length(deathtime)){
          nrisk[i]=sum((deathtime[i]<=T)*wgt)
          ndeath[i]=sum((T==deathtime[i])*censor*wgt)
        }
        prodobj=1-ndeath/nrisk
        survp=rep(0,length(deathtime))
        for(i in 1:length(deathtime)){survp[i]=prod(prodobj[1:i])}
        return(data.frame(cbind(deathtime,ndeath,nrisk,survp)))}

      for (j in 1:ne){
        # generating perturbation variable
        E = rexp(n,1)
        if (D==rep(1,n)){
          W_star = rep(1,n)
        } else {
          Gest=SC.func(Z,1-D,E)
          SC=stepfun(Gest$deathtime,c(1,Gest$survp))
          W_star <- D / Gest$survp[findInterval(Z, Gest$deathtime)]*Gest$survp[min(which(floor(Gest$deathtime)==(t_0)))]
          W_star[is.na(W_star)] <- max(W_star, na.rm = TRUE)
        }
        result = rev_isObj(coefficient, X, W_star, H, E, I, logT, Q)
        rcpp.result.ismb = cbind(rcpp.result.ismb,result)
      }
      v = cov(t(rcpp.result.ismb))
      rcpp.a = Amat(coefficient, X, W, H, E, I, logT, Q)
      inva <- try(solve(rcpp.a))
      if(class(inva)[1] == "try-error"){
        se = rep(NA,nc+1)
        stop("Slope matrix is singular matrix")
      } else {
        sigma = t(inva) %*% v %*% inva
        se = sqrt(diag(sigma))
      }
      list(coefficient=coefficient, stderr = se)
    } else {
      coefficient = c(NA,rep(NA,nc))
      se = c(NA,rep(NA,nc))
      list(coefficient=coefficient, stderr = se)
    }
  }
}
