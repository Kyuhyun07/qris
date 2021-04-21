#### R package ####
library(quantreg)
library(survival)
library(emplik)
library(xtable)
library(nleqslv)
library(BB)
library(Rcpp)
library(RcppArmadillo)
library(scales)

#### method 1 : Original code (ISMB + rcpp's weight + weight-in, only R use) ####
ismbw1.est= function(Z, nc, covariate, D, t_0, Q, ne, init){
  n = length(Z)
  data = matrix(NA, n, nc+5)
  data[,1] = Z
  data[,2] = log(Z-t_0)
  data[,3] = as.numeric(Z>=t_0)
  data[,4:(nc+3)] = covariate
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
  
  # Objective equation
  objectF = function(beta){
    result = t(X*I) %*% {W*(pnorm((X%*%beta-logT)/sqrt(diag(X %*% H %*% t(X)))))-Q}
    result = as.vector(result)
  }
  
  #### revised object equation ####
  rev.objectF = function(beta){
    result = t(eta*X*I) %*% {W*(pnorm((X%*%beta-logT)/sqrt(diag(X %*% H %*% t(X)))))-Q}
    result = as.vector(result)
  }
  
  # Covariate setting (1 covariate)
  X = as.matrix(cbind(c(rep(1,n)),data[,4:(nc+3)]))
  W = data[,(nc+5)]
  logT = data[,2]
  I = data[,3]
  H = diag(1/n, nc+1, nc+1)
  
  # Guess beta option (rq : solution from rq, 0 : no effect of covariate, others : covariate)
  if (init == "rq"){
    betastart = as.vector(rq(data[,2] ~ as.matrix(data[,4:(nc+3)]), tau=Q, weight=W)$coefficient)
  } else if (init == 0){
    betastart = c(1,rep(0,nc))
  } else {
    betastart = rnorm(nc+1)}
  is.fit = nleqslv(betastart, objectF)
  
  if (is.fit$termcd == 1 | is.fit$termcd == 2){
    coefficient = is.fit$x
    # Variance estimation : ISMB
    result.ismb=c()
    for (j in 1:ne){
      eta = rexp(n,1)
      result = t(eta*X*I) %*% {W*(pnorm((X%*%coefficient-logT)/sqrt(diag(X %*% H %*% t(X)))))-Q}
      result.ismb = cbind(result.ismb,result)
    }
    v = cov(t(result.ismb))
    a.beta = t(X*I*W*as.vector(dnorm((logT-X%*%coefficient)/sqrt(diag(X %*% H %*% t(X))))))%*%(-X/sqrt(diag(X %*% H %*% t(X))))
    inva.beta = qr.solve(a.beta)
    sigma = t(inva.beta) %*% v %*% inva.beta
    SE = sqrt(diag(sigma))
    beta.se = cbind(coefficient, SE)
    print(beta.se)
  } else {
    coefficient = c(NA,rep(NA,nc))
    SE = c(NA,rep(NA,nc))
    beta.se = cbind(coefficient, SE)
    print(beta.se)
  }
}

#### method 2 : Rcpp ISMB (R package) ####
#' @param Z is a vector of observed time, which is minimum of failure time and censored time
#' @param nc is a number of covariates used in analysis
#' @param covariate is a matrix of covariate (# row = # of subject, # of column = # of covariate(nc))
#' @param D is a vector of censoring indicator (1 = not censored, 0 = censored)
#' @param t_0 is the followup time(or basetime of analysis)
#' @param Q is the quantile
#' @param ne is number of multiplier bootstrapping for V matrix estimation
#' @param init is option for initial guess of regression parameter ("rq" = a solution of rq function, 0 assumes no treatment effect, otherwise assume treatment effect as 1)
#' Induced smoothing estimating equation, adopted from is_objectF()
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

#' Induce smoothing estimating equation, adopted from rev_is_objectF()
sourceCpp(code = '
          #include <RcppArmadillo.h>
          // [[Rcpp::depends(RcppArmadillo)]]
          using namespace arma;
          // [[Rcpp::export]]
          arma::mat rev_isObj(arma::vec b, arma::mat X, arma::vec W, arma::mat H, arma::vec E,arma::vec I,
          arma::vec logT, double Q) {
          arma::mat m1 = X % repmat(I, 1, X.n_cols) % repmat(E, 1, X.n_cols);
          arma::mat m2 = normcdf((X * b - logT) / sqrt(diagvec(X * H * X.t()))) % W - Q;
          return m1.t() * m2;
          }')

#' Induce smoothing estimating equation, adopted from Amat
sourceCpp(code = '
          #include <RcppArmadillo.h>
          // [[Rcpp::depends(RcppArmadillo)]]
          using namespace arma;
          // [[Rcpp::export]]
          arma::mat Amat(arma::vec b, arma::mat X, arma::vec W, arma::mat H, arma::vec E, arma::vec I,
          arma::vec logT, double Q) {
          arma::mat m1 = X % repmat(I, 1, X.n_cols) % repmat(W, 1, X.n_cols);
          arma::mat m2 =  (-X / repmat(sqrt(diagvec(X * H * X.t())), 1, X.n_cols)) % repmat(normpdf((X * b - logT) / sqrt(diagvec(X * H * X.t()))), 1, X.n_cols);
          return m1.t() * m2;
          }')

rcpp.est= function(Z, nc, covariate, D, t_0, Q, ne, init){
  n = length(Z)
  data = matrix(NA, n, nc+5)
  data[,1] = Z
  data[,2] = log(Z-t_0)
  data[,3] = as.numeric(Z>t_0)
  data[,4:(nc+3)] = covariate
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
  if (init == "rq"){
    betastart = as.vector(rq(data[,2] ~ as.matrix(data[,4:(nc+3)]), tau=Q, weight=W)$coefficient)
  } else if (init == 0){
    betastart = c(1,rep(0,nc))
  } else {
    betastart = rnorm(nc+1)}
  rcpp.fit <- nleqslv(betastart, function(b) isObj(b, X, W, H, I, logT, Q))
  
  if (rcpp.fit$termcd == 1 | rcpp.fit$termcd == 2){
    coefficient = rcpp.fit$x
    # Variance estimation : ISMB
    rcpp.result.ismb=c()
    for (j in 1:ne){
      E = rexp(n,1)
      result = rev_isObj(coefficient, X, W, H, E, I, logT, Q)
      rcpp.result.ismb = cbind(rcpp.result.ismb,result)
    }
    v = cov(t(rcpp.result.ismb))
    rcpp.a = Amat(coefficient, X, W, H, E, I, logT, Q)
    inva = qr.solve(rcpp.a)
    sigma = t(inva) %*% v %*% inva
    SE = sqrt(diag(sigma))
    beta.se = cbind(coefficient, SE)
    rownames(beta.se) = c("Intercept", covar)
    print(beta.se)
  } else {
    coefficient = c(NA,rep(NA,nc))
    SE = c(NA,rep(NA,nc))
    beta.se = cbind(coefficient, SE)
    rownames(beta.se) = c("Intercept", covar)
    print(beta.se)
  }
}

#### Simulation data (1 covariate) ####
#### Data Generation function ####
data.gen<-function(samplesize, censor){
  sim=matrix(NA,samplesize,5)
  colnames(sim) = c("T","C","Z","X","censored")
  # Generate C_i
  sim[,2] = runif(samplesize,0,censor)
  # Covariates (Control=0, Treatment=1)
  sim[,4] = rbinom(samplesize,size=1,p=0.5)
  # Generate T_i (Given Condition r=rho_0, k=2, exp(beta_0)=5, exp(beta_1)=2))
  unif = runif(n=samplesize ,min = 0,max = 1)
  for (q in 1:samplesize){
    if (sim[q,4]==0){
      sim[q,1]={{-log(1-unif[q])}^(1/k)}/r.initial.0
    } else {
      sim[q,1]={{-log(1-unif[q])}^(1/k)}/r.initial.1
    }
  }
  # Generate Y_i (min(T,C))
  sim[,3] = apply(sim[,1:2], 1, FUN=min)
  # Censoring indicator (Censored=0, Not censored=1)
  sim[,5]=I(sim[,1]<sim[,2])
  # Ordering
  sim = sim[order(sim[,3]),]
  n = nrow(sim)
  sim = as.data.frame(sim)
  return(sim)
}

#### Given information ####
exp.beta.initial.0=5
exp.beta.initial.1=10
k=2
r.initial.0=(log(10/5))^(1/k)/exp.beta.initial.0
r.initial.1=(log(10/5))^(1/k)/exp.beta.initial.1

#### Example parameter ####
c.0=5000000
c.1=70.39
c.3=24.35
c.5=14.07
c.7=8.49
a<-data.gen(200,c.3)
Z=a[,3]
nc=1
covariate=a[,4]
D=a[,5]
t_0=1
Q=0.5
ne=200

#### simulation result ###
# Original code
# guess w/o covariate effect : 0.11 sec
tic()
ismbw1.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200, 0)
toc()
# guess w/ covariate effect : 0.15 sec
tic()
ismbw1.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200, 1)
toc()
# Use rq solution : 0.27 sec
tic()
ismbw1.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200, "rq")
toc()

# Rcpp code
# guess w/o covariate effect : 0.03 sec
tic()
rcpp.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200, 0)
toc()
# guess w covariate effect : 0.03 sec
tic()
rcpp.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200, 1)
toc()
# Use rq solution : 0.04 sec
tic()
rcpp.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200, "rq")
toc()


#### Real data : Diabetic retinopathy study data ####
data("retinopathy")
reti = retinopathy
reti_rev = reti
# covariate 1 : Risk group (rescale 6 to 12 -> 0.5 to 1)
reti_rev$risk = rescale(reti$risk, to=c(0,1), from=c(0,12))
# covariate 2 : Age (no change, continuous variable)
# covariate 3 : Type of Diabetes (adult = 1, juvenile = 0)
reti_rev$type = as.numeric(reti$type)-1
# covariate 4 : treatment (no change)
# covariate 5 : interaction (type x trt)
reti_rev$typextrt = reti_rev$type*reti$trt

Z=reti_rev$futime
nc=5
covariate=as.matrix(reti_rev[,c(9,4,5,6,10)])
D=reti_rev$status
t_0=1
Q=0.25
ne=200

# Original code
# guess w/o covariate effect : 0.39 sec
tic()
ismbw1.est(Z, nc, covariate, D, t_0, Q, ne, 0)
toc()
# guess w/ covariate effect : 0.36 sec
tic()
ismbw1.est(Z, nc, covariate, D, t_0, Q, ne, 1)
toc()
# Use rq solution : 0.32 sec
tic()
ismbw1.est(Z, nc, covariate, D, t_0, Q, ne, "rq")
toc()

# Rcpp code
# guess w/o covariate effect : 0.33 sec
tic()
rcpp.est(Z, nc, covariate, D, t_0, Q, ne, 0)
toc()
# guess w covariate effect : 0.29 sec
tic()
rcpp.est(Z, nc, covariate, D, t_0, Q, ne, 1)
toc()
# Use rq solution : 0.24 sec
tic()
rcpp.est(Z, nc, covariate, D, t_0, Q, ne, "rq")
toc()