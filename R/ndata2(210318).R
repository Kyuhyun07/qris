#### ndata2(210318) : Condition ####
# data size = 200
# beta0, beta1 effective
# Quantile 50%
# simulation 2000
# eta = 200
# method 1 : Original simulation code (ISMB + Li's weight + weight-in + dfsane)
# method 2 : rcpp's method (Li's weight + weight-in + dfsane)

library(quantreg)
library(survival)
library(emplik)
library(xtable)
library(nleqslv)
library(BB)
library(Rcpp)
library(RcppArmadillo)
library(tictoc)


#### True Beta ####
#beta_0    beta_1
#t_0=0 1.609438 0.6931472
#t_0=1 1.410748 0.7974189
#t_0=2 1.219403 0.9070615
#t_0=3 1.040613 1.0174711

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

#### Given information(My assumption) ####
exp.beta.initial.0=5
exp.beta.initial.1=10
k=2
r.initial.0=(log(10/5))^(1/k)/exp.beta.initial.0
r.initial.1=(log(10/5))^(1/k)/exp.beta.initial.1

#### Example parameter ####
# c.0=5000000
# c.1=70.39
# c.3=24.35
# c.5=14.07
# c.7=8.49
# a<-data.gen(200,c.3)
# Z=a[,3]
# nc=1
# covariate=a[,4]
# D=a[,5]
# t_0=1
# Q=0.5
# ne=200

#### Estimation functions ####
# method 1 : Our method (ISMB + rcpp's weight + weight-in)
ismbw1.est= function(Z, nc, covariate, D, t_0, Q, ne){
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
  for (j in 4:(nc+3)){
    colnames(data)[j]="covariate"
  }
  colnames(data)[(nc+4):(nc+5)] = c("delta","Li's weight")

  # Original Li's weight with jump weight
  fit = WKM(data[,1], (data[,(nc+4)]), zc = rep(1,n), w = rep(1,n))
  fit2 = survfit(Surv(data[,1],1-(data[,(nc+4)])) ~ 1)
  G_t <- fit2$surv[min(which(floor(fit2$time)==t_0))]
  data[,(nc+5)] = fit$jump*n*G_t

  # # Rcpp Li's weight with jump weight
  # sv <- survfit(Surv(data[,1], 1 - data[,(nc+4)]) ~ 1)
  # W <- data[,(nc+4)] / sv$surv[findInterval(data[,1], sv$time)]*sv$surv[min(which(floor(sv$time)==(t_0)))]
  # W[is.na(W)] <- max(W, na.rm = TRUE)
  # data[,(nc+5)] = W

  # # Li's weight with survfit
  fit = WKM(data[,1], 1-(data[,(nc+4)]), zc = rep(1,n), w = rep(1,n))
  fit2 = survfit(Surv(data[,1],1-(data[,(nc+4)])) ~ 1)
  G_t <- fit2$surv[min(which(floor(fit2$time)==t_0))]
  # data[,(nc+5)] = data[,(nc+4)]/fit$surv*G_t
  dummy2 = data[,(nc+4)]/fit$surv*G_t
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

  # Change betastart when real data analysis c(1,rep(1,nc))
  # In simulation, betastart = true beta
  if (t_0 == 0){
    betastart = c(1.609438,0.6931472)
  } else if (t_0 == 1) {
    betastart = c(1.410748, 0.7974189)
  } else if (t_0 == 2) {
    betastart = c(1.219403, 0.9070615)
  } else {
    betastart = c(1.040613, 1.0174711)}
  is.fit = nleqslv(betastart, objectF)

  if (is.fit$termcd == 1){
    solbeta = is.fit$x
    # Variance estimation : ISMB
    result.ismb=c()
    set.seed(11)
    for (j in 1:ne){
      eta = rexp(n,1)
      result = t(eta*X*I) %*% {W*(pnorm((X%*%solbeta-logT)/sqrt(diag(X %*% H %*% t(X)))))-Q}
      result.ismb = cbind(result.ismb,result)
    }
    v = cov(t(result.ismb))
    a.beta = t(X*I*W*as.vector(dnorm((logT-X%*%solbeta)/sqrt(diag(X %*% H %*% t(X))))))%*%(-X/sqrt(diag(X %*% H %*% t(X))))
    inva.beta = qr.solve(a.beta)
    sigma = t(inva.beta) %*% v %*% inva.beta
    sd = sqrt(diag(sigma))
    beta.sd = cbind(solbeta, sd)
    print(beta.sd)
  } else {
    solbeta = c(NA,NA)
    sd = c(NA,NA)
    beta.sd = cbind(solbeta, sd)
    print(beta.sd)
  }
}


# method 2. Rcpp ISMB
#' Induce smoothing estimating equation, adopted from is_objectF()
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
result = t(eta*X*I) %*% {W*(pnorm((X%*%solbeta-logT)/sqrt(diag(X %*% H %*% t(X)))))-Q}

a.beta = t(X*I*W*as.vector(dnorm((X%*%solbeta-logT)/sqrt(diag(X %*% H %*% t(X))))))%*%(-X/sqrt(diag(X %*% H %*% t(X))))
#' Induce smoothing estimating equation, adopted from Amat
sourceCpp(code = '
    #include <RcppArmadillo.h>
    // [[Rcpp::depends(RcppArmadillo)]]
    using namespace arma;
    // [[Rcpp::export]]
    arma::mat Amat(arma::vec b, arma::mat X, arma::vec W, arma::mat H, arma::vec E, arma::vec I,
                    arma::vec logT, double Q) {
    arma::mat m1 = X % repmat(I, 1, X.n_cols) % repmat(W, 1, X.n_cols);
    arma::mat m2 = (X / sqrt(diagvec(X * H * X.t())))) % repmat(normpdf((X * b - logT) / sqrt(diagvec(X * H * X.t()))), 1, X.n_cols);
    return m1.t() * m2;
  }')

rcpp.est= function(Z, nc, covariate, D, t_0, Q, ne){
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
  for (j in 4:(nc+3)){
    colnames(data)[j]="covariate"
  }
  colnames(data)[(nc+4):(nc+5)] = c("delta","rcpp's weight")

  # Original Li's weight with jump weight
  fit = WKM(data[,1], (data[,(nc+4)]), zc = rep(1,n), w = rep(1,n))
  fit2 = survfit(Surv(data[,1],1-(data[,(nc+4)])) ~ 1)
  G_t <- fit2$surv[min(which(floor(fit2$time)==t_0))]
  data[,(nc+5)] = fit$jump*n*G_t

  # # Rcpp Li's weight with jump weight
  # sv <- survfit(Surv(data[,1], 1 - data[,(nc+4)]) ~ 1)
  # W <- data[,(nc+4)] / sv$surv[findInterval(data[,1], sv$time)]*sv$surv[min(which(floor(sv$time)==(t_0)))]
  # W[is.na(W)] <- max(W, na.rm = TRUE)
  # data[,(nc+5)] = W

  # # Li's weight with survfit
  # fit = WKM(data[,1], 1-(data[,(nc+4)]), zc = rep(1,n), w = rep(1,n))
  # fit2 = survfit(Surv(data[,1],1-(data[,(nc+4)])) ~ 1)
  # G_t <- fit2$surv[min(which(floor(fit2$time)==t_0))]
  # data[,(nc+5)] = data[,(nc+4)]/fit$surv*G_t

  # Covariate setting (1 covariate)
  X = as.matrix(cbind(c(rep(1,n)),data[,4:(nc+3)]))
  W = data[,(nc+5)]
  logT = data[,2]
  I = data[,3]
  H = diag(1/n, nc+1, nc+1)

  # Change betastart when real data analysis c(1,rep(1,nc))
  if (t_0 == 0){
    betastart = c(1.609438,0.6931472)
  } else if (t_0 == 1) {
    betastart = c(1.410748, 0.7974189)
  } else if (t_0 == 2) {
    betastart = c(1.219403, 0.9070615)
  } else {
    betastart = c(1.040613, 1.0174711)}
  rcpp.fit <- nleqslv(betastart, function(b) isObj(b, X, W, H, I, logT, Q))

  if (rcpp.fit$termcd == 1){
    solbeta = rcpp.fit$x
    # Variance estimation : ISMB
    result.ismb=c()
    set.seed(11)
    for (j in 1:ne){
      E = rexp(n,1)
      result = rev_isObj(solbeta, X, W, H, E, I, logT, Q)
      result.ismb = cbind(result.ismb,result)
    }
    v = cov(t(result.ismb))
    a.beta = t(X*I*W*as.vector(dnorm((X%*%solbeta-logT)/sqrt(diag(X %*% H %*% t(X))))))%*%(-X/sqrt(diag(X %*% H %*% t(X))))
    inva.beta = qr.solve(a.beta)
    sigma = t(inva.beta) %*% v %*% inva.beta
    sd = sqrt(diag(sigma))
    beta.sd = cbind(solbeta, sd)
    print(beta.sd)
  } else {
    solbeta = c(NA,NA)
    sd = c(NA,NA)
    beta.sd = cbind(solbeta, sd)
    print(beta.sd)
  }
}


#### Indicator function ####
ind=function(a,b,c){
  if (a>=b&a<=c) {
    result=1
  } else {
    result=0
  }
  print(result)
}

table0.isw1<-matrix(NA,5,8)
rownames(table0.isw1)<-c(0,10,30,50,70)
colnames(table0.isw1)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table1.isw1<-matrix(NA,5,8)
rownames(table1.isw1)<-c(0,10,30,50,70)
colnames(table1.isw1)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table2.isw1<-matrix(NA,5,8)
rownames(table2.isw1)<-c(0,10,30,50,70)
colnames(table2.isw1)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table3.isw1<-matrix(NA,5,8)
rownames(table3.isw1)<-c(0,10,30,50,70)
colnames(table3.isw1)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")

table0.rcpp<-matrix(NA,5,8)
rownames(table0.rcpp)<-c(0,10,30,50,70)
colnames(table0.rcpp)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table1.rcpp<-matrix(NA,5,8)
rownames(table1.rcpp)<-c(0,10,30,50,70)
colnames(table1.rcpp)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table2.rcpp<-matrix(NA,5,8)
rownames(table2.rcpp)<-c(0,10,30,50,70)
colnames(table2.rcpp)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table3.rcpp<-matrix(NA,5,8)
rownames(table3.rcpp)<-c(0,10,30,50,70)
colnames(table3.rcpp)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")


#### censoring point at t_0=0 ####
c.0=5000000
c.1=78.11
c.3=26.36
c.5=15.08
c.7=9.09
#### t_0=0 & c=0% ####
b0.isw1.00 = c()
b0.isw1.sd.00 = c()
b1.isw1.00 = c()
b1.isw1.sd.00 = c()
cover.isw1.00=matrix(NA,2000,8)
colnames(cover.isw1.00)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.rcpp.00 = c()
b0.rcpp.sd.00 = c()
b1.rcpp.00 = c()
b1.rcpp.sd.00 = c()
cover.rcpp.00=matrix(NA,2000,8)
colnames(cover.rcpp.00)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.0)
  # method 1 : Our method (ISMB + rcpp's weight + weight-in + dfsane)
  tryCatch({
    tic()
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    toc()
    b0.isw1.00[i] = ismb.fit[1,1]
    b0.isw1.sd.00[i] = ismb.fit[1,2]
    b1.isw1.00[i] = ismb.fit[2,1]
    b1.isw1.sd.00[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.00[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.00[i,2] = ismb.fit[1,1]
    cover.isw1.00[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.00[i,4] = ind(1.609438, cover.isw1.00[i,1], cover.isw1.00[i,3])
    cover.isw1.00[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.00[i,6] = ismb.fit[2,1]
    cover.isw1.00[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.00[i,8] = ind(0.6931472, cover.isw1.00[i,5], cover.isw1.00[i,7])}
    , error=function(e){
      b0.isw1.00[i] = NA
      b0.isw1.sd.00[i] = NA
      b1.isw1.00[i] = NA
      b1.isw1.sd.00[i] = NA
      # Coverage
      cover.isw1.00[i,1] = NA
      cover.isw1.00[i,2] = NA
      cover.isw1.00[i,3] = NA
      cover.isw1.00[i,4] = NA
      cover.isw1.00[i,5] = NA
      cover.isw1.00[i,6] = NA
      cover.isw1.00[i,7] = NA
      cover.isw1.00[i,8] = NA
    })

  # method 2 : rcpp's method (rcpp's weight + weight-in + dfsane)
  tryCatch({
    rcpp.fit = rcpp.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.rcpp.00[i] = rcpp.fit[1,1]
    b0.rcpp.sd.00[i] = rcpp.fit[1,2]
    b1.rcpp.00[i] = rcpp.fit[2,1]
    b1.rcpp.sd.00[i] = rcpp.fit[2,2]
    # Coverage
    cover.rcpp.00[i,1] = rcpp.fit[1,1]-1.96*rcpp.fit[1,2]
    cover.rcpp.00[i,2] = rcpp.fit[1,1]
    cover.rcpp.00[i,3] = rcpp.fit[1,1]+1.96*rcpp.fit[1,2]
    cover.rcpp.00[i,4] = ind(1.609438, cover.rcpp.00[i,1], cover.rcpp.00[i,3])
    cover.rcpp.00[i,5] = rcpp.fit[2,1]-1.96*rcpp.fit[2,2]
    cover.rcpp.00[i,6] = rcpp.fit[2,1]
    cover.rcpp.00[i,7] = rcpp.fit[2,1]+1.96*rcpp.fit[2,2]
    cover.rcpp.00[i,8] = ind(0.6931472, cover.rcpp.00[i,5], cover.rcpp.00[i,7])}
    , error=function(e){
      b0.rcpp.00[i] = NA
      b0.rcpp.sd.00[i] = NA
      b1.rcpp.00[i] = NA
      b1.rcpp.sd.00[i] = NA
      # Coverage
      cover.rcpp.00[i,1] = NA
      cover.rcpp.00[i,2] = NA
      cover.rcpp.00[i,3] = NA
      cover.rcpp.00[i,4] = NA
      cover.rcpp.00[i,5] = NA
      cover.rcpp.00[i,6] = NA
      cover.rcpp.00[i,7] = NA
      cover.rcpp.00[i,8] = NA
    })
}
# IS beta table
table0.isw1[1,1]<-mean(b0.isw1.00,na.rm=TRUE)
table0.isw1[1,2]<-mean(b0.isw1.sd.00,na.rm=TRUE)
table0.isw1[1,3]<-sd(b0.isw1.00,na.rm=TRUE)
table0.isw1[1,4]<-mean(cover.isw1.00[,4],na.rm=TRUE)
table0.isw1[1,5]<-mean(b1.isw1.00,na.rm=TRUE)
table0.isw1[1,6]<-mean(b1.isw1.sd.00,na.rm=TRUE)
table0.isw1[1,7]<-sd(b1.isw1.00,na.rm=TRUE)
table0.isw1[1,8]<-mean(cover.isw1.00[,8],na.rm=TRUE)

table0.rcpp[1,1]<-mean(b0.rcpp.00,na.rm=TRUE)
table0.rcpp[1,2]<-mean(b0.rcpp.sd.00,na.rm=TRUE)
table0.rcpp[1,3]<-sd(b0.rcpp.00,na.rm=TRUE)
table0.rcpp[1,4]<-mean(cover.rcpp.00[,4],na.rm=TRUE)
table0.rcpp[1,5]<-mean(b1.rcpp.00,na.rm=TRUE)
table0.rcpp[1,6]<-mean(b1.rcpp.sd.00,na.rm=TRUE)
table0.rcpp[1,7]<-sd(b1.rcpp.00,na.rm=TRUE)
table0.rcpp[1,8]<-mean(cover.rcpp.00[,8],na.rm=TRUE)

#### t_0=0 & c=10% ####
b0.isw1.01 = c()
b0.isw1.sd.01 = c()
b1.isw1.01 = c()
b1.isw1.sd.01 = c()
cover.isw1.01=matrix(NA,2000,8)
colnames(cover.isw1.01)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.rcpp.01 = c()
b0.rcpp.sd.01 = c()
b1.rcpp.01 = c()
b1.rcpp.sd.01 = c()
cover.rcpp.01=matrix(NA,2000,8)
colnames(cover.rcpp.01)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.1)
  # method 1 : Our method (ISMB + rcpp's weight + weight-in + dfsane)
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.isw1.01[i] = ismb.fit[1,1]
    b0.isw1.sd.01[i] = ismb.fit[1,2]
    b1.isw1.01[i] = ismb.fit[2,1]
    b1.isw1.sd.01[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.01[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.01[i,2] = ismb.fit[1,1]
    cover.isw1.01[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.01[i,4] = ind(1.609438, cover.isw1.01[i,1], cover.isw1.01[i,3])
    cover.isw1.01[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.01[i,6] = ismb.fit[2,1]
    cover.isw1.01[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.01[i,8] = ind(0.6931472, cover.isw1.01[i,5], cover.isw1.01[i,7])}
    , error=function(e){
      b0.isw1.01[i] = NA
      b0.isw1.sd.01[i] = NA
      b1.isw1.01[i] = NA
      b1.isw1.sd.01[i] = NA
      # Coverage
      cover.isw1.01[i,1] = NA
      cover.isw1.01[i,2] = NA
      cover.isw1.01[i,3] = NA
      cover.isw1.01[i,4] = NA
      cover.isw1.01[i,5] = NA
      cover.isw1.01[i,6] = NA
      cover.isw1.01[i,7] = NA
      cover.isw1.01[i,8] = NA
    })

  # method 2 : rcpp's method (rcpp's weight + weight-in + dfsane)
  tryCatch({
    rcpp.fit = rcpp.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.rcpp.01[i] = rcpp.fit[1,1]
    b0.rcpp.sd.01[i] = rcpp.fit[1,2]
    b1.rcpp.01[i] = rcpp.fit[2,1]
    b1.rcpp.sd.01[i] = rcpp.fit[2,2]
    # Coverage
    cover.rcpp.01[i,1] = rcpp.fit[1,1]-1.96*rcpp.fit[1,2]
    cover.rcpp.01[i,2] = rcpp.fit[1,1]
    cover.rcpp.01[i,3] = rcpp.fit[1,1]+1.96*rcpp.fit[1,2]
    cover.rcpp.01[i,4] = ind(1.609438, cover.rcpp.01[i,1], cover.rcpp.01[i,3])
    cover.rcpp.01[i,5] = rcpp.fit[2,1]-1.96*rcpp.fit[2,2]
    cover.rcpp.01[i,6] = rcpp.fit[2,1]
    cover.rcpp.01[i,7] = rcpp.fit[2,1]+1.96*rcpp.fit[2,2]
    cover.rcpp.01[i,8] = ind(0.6931472, cover.rcpp.01[i,5], cover.rcpp.01[i,7])}
    , error=function(e){
      b0.rcpp.01[i] = NA
      b0.rcpp.sd.01[i] = NA
      b1.rcpp.01[i] = NA
      b1.rcpp.sd.01[i] = NA
      # Coverage
      cover.rcpp.01[i,1] = NA
      cover.rcpp.01[i,2] = NA
      cover.rcpp.01[i,3] = NA
      cover.rcpp.01[i,4] = NA
      cover.rcpp.01[i,5] = NA
      cover.rcpp.01[i,6] = NA
      cover.rcpp.01[i,7] = NA
      cover.rcpp.01[i,8] = NA
    })
}
# IS beta table
table0.isw1[2,1]<-mean(b0.isw1.01,na.rm=TRUE)
table0.isw1[2,2]<-mean(b0.isw1.sd.01,na.rm=TRUE)
table0.isw1[2,3]<-sd(b0.isw1.01,na.rm=TRUE)
table0.isw1[2,4]<-mean(cover.isw1.01[,4],na.rm=TRUE)
table0.isw1[2,5]<-mean(b1.isw1.01,na.rm=TRUE)
table0.isw1[2,6]<-mean(b1.isw1.sd.01,na.rm=TRUE)
table0.isw1[2,7]<-sd(b1.isw1.01,na.rm=TRUE)
table0.isw1[2,8]<-mean(cover.isw1.01[,8],na.rm=TRUE)

table0.rcpp[2,1]<-mean(b0.rcpp.01,na.rm=TRUE)
table0.rcpp[2,2]<-mean(b0.rcpp.sd.01,na.rm=TRUE)
table0.rcpp[2,3]<-sd(b0.rcpp.01,na.rm=TRUE)
table0.rcpp[2,4]<-mean(cover.rcpp.01[,4],na.rm=TRUE)
table0.rcpp[2,5]<-mean(b1.rcpp.01,na.rm=TRUE)
table0.rcpp[2,6]<-mean(b1.rcpp.sd.01,na.rm=TRUE)
table0.rcpp[2,7]<-sd(b1.rcpp.01,na.rm=TRUE)
table0.rcpp[2,8]<-mean(cover.rcpp.01[,8],na.rm=TRUE)

#### t_0=0 & c=30% ####
b0.isw1.03 = c()
b0.isw1.sd.03 = c()
b1.isw1.03 = c()
b1.isw1.sd.03 = c()
cover.isw1.03=matrix(NA,2000,8)
colnames(cover.isw1.03)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.rcpp.03 = c()
b0.rcpp.sd.03 = c()
b1.rcpp.03 = c()
b1.rcpp.sd.03 = c()
cover.rcpp.03=matrix(NA,2000,8)
colnames(cover.rcpp.03)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.3)
  # method 1 : Our method (ISMB + rcpp's weight + weight-in + dfsane)
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.isw1.03[i] = ismb.fit[1,1]
    b0.isw1.sd.03[i] = ismb.fit[1,2]
    b1.isw1.03[i] = ismb.fit[2,1]
    b1.isw1.sd.03[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.03[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.03[i,2] = ismb.fit[1,1]
    cover.isw1.03[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.03[i,4] = ind(1.609438, cover.isw1.03[i,1], cover.isw1.03[i,3])
    cover.isw1.03[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.03[i,6] = ismb.fit[2,1]
    cover.isw1.03[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.03[i,8] = ind(0.6931472, cover.isw1.03[i,5], cover.isw1.03[i,7])}
    , error=function(e){
      b0.isw1.03[i] = NA
      b0.isw1.sd.03[i] = NA
      b1.isw1.03[i] = NA
      b1.isw1.sd.03[i] = NA
      # Coverage
      cover.isw1.03[i,1] = NA
      cover.isw1.03[i,2] = NA
      cover.isw1.03[i,3] = NA
      cover.isw1.03[i,4] = NA
      cover.isw1.03[i,5] = NA
      cover.isw1.03[i,6] = NA
      cover.isw1.03[i,7] = NA
      cover.isw1.03[i,8] = NA
    })

  # method 2 : rcpp's method (rcpp's weight + weight-in + dfsane)
  tryCatch({
    rcpp.fit = rcpp.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.rcpp.03[i] = rcpp.fit[1,1]
    b0.rcpp.sd.03[i] = rcpp.fit[1,2]
    b1.rcpp.03[i] = rcpp.fit[2,1]
    b1.rcpp.sd.03[i] = rcpp.fit[2,2]
    # Coverage
    cover.rcpp.03[i,1] = rcpp.fit[1,1]-1.96*rcpp.fit[1,2]
    cover.rcpp.03[i,2] = rcpp.fit[1,1]
    cover.rcpp.03[i,3] = rcpp.fit[1,1]+1.96*rcpp.fit[1,2]
    cover.rcpp.03[i,4] = ind(1.609438, cover.rcpp.03[i,1], cover.rcpp.03[i,3])
    cover.rcpp.03[i,5] = rcpp.fit[2,1]-1.96*rcpp.fit[2,2]
    cover.rcpp.03[i,6] = rcpp.fit[2,1]
    cover.rcpp.03[i,7] = rcpp.fit[2,1]+1.96*rcpp.fit[2,2]
    cover.rcpp.03[i,8] = ind(0.6931472, cover.rcpp.03[i,5], cover.rcpp.03[i,7])}
    , error=function(e){
      b0.rcpp.03[i] = NA
      b0.rcpp.sd.03[i] = NA
      b1.rcpp.03[i] = NA
      b1.rcpp.sd.03[i] = NA
      # Coverage
      cover.rcpp.03[i,1] = NA
      cover.rcpp.03[i,2] = NA
      cover.rcpp.03[i,3] = NA
      cover.rcpp.03[i,4] = NA
      cover.rcpp.03[i,5] = NA
      cover.rcpp.03[i,6] = NA
      cover.rcpp.03[i,7] = NA
      cover.rcpp.03[i,8] = NA
    })
}

# IS beta table
table0.isw1[3,1]<-mean(b0.isw1.03,na.rm=TRUE)
table0.isw1[3,2]<-mean(b0.isw1.sd.03,na.rm=TRUE)
table0.isw1[3,3]<-sd(b0.isw1.03,na.rm=TRUE)
table0.isw1[3,4]<-mean(cover.isw1.03[,4],na.rm=TRUE)
table0.isw1[3,5]<-mean(b1.isw1.03,na.rm=TRUE)
table0.isw1[3,6]<-mean(b1.isw1.sd.03,na.rm=TRUE)
table0.isw1[3,7]<-sd(b1.isw1.03,na.rm=TRUE)
table0.isw1[3,8]<-mean(cover.isw1.03[,8],na.rm=TRUE)

table0.rcpp[3,1]<-mean(b0.rcpp.03,na.rm=TRUE)
table0.rcpp[3,2]<-mean(b0.rcpp.sd.03,na.rm=TRUE)
table0.rcpp[3,3]<-sd(b0.rcpp.03,na.rm=TRUE)
table0.rcpp[3,4]<-mean(cover.rcpp.03[,4],na.rm=TRUE)
table0.rcpp[3,5]<-mean(b1.rcpp.03,na.rm=TRUE)
table0.rcpp[3,6]<-mean(b1.rcpp.sd.03,na.rm=TRUE)
table0.rcpp[3,7]<-sd(b1.rcpp.03,na.rm=TRUE)
table0.rcpp[3,8]<-mean(cover.rcpp.03[,8],na.rm=TRUE)

#### t_0=0 & c=50% ####
b0.isw1.05 = c()
b0.isw1.sd.05 = c()
b1.isw1.05 = c()
b1.isw1.sd.05 = c()
cover.isw1.05=matrix(NA,2000,8)
colnames(cover.isw1.05)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.rcpp.05 = c()
b0.rcpp.sd.05 = c()
b1.rcpp.05 = c()
b1.rcpp.sd.05 = c()
cover.rcpp.05=matrix(NA,2000,8)
colnames(cover.rcpp.05)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.5)
  # method 1 : Our method (ISMB + rcpp's weight + weight-in + dfsane)
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.isw1.05[i] = ismb.fit[1,1]
    b0.isw1.sd.05[i] = ismb.fit[1,2]
    b1.isw1.05[i] = ismb.fit[2,1]
    b1.isw1.sd.05[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.05[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.05[i,2] = ismb.fit[1,1]
    cover.isw1.05[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.05[i,4] = ind(1.609438, cover.isw1.05[i,1], cover.isw1.05[i,3])
    cover.isw1.05[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.05[i,6] = ismb.fit[2,1]
    cover.isw1.05[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.05[i,8] = ind(0.6931472, cover.isw1.05[i,5], cover.isw1.05[i,7])}
    , error=function(e){
      b0.isw1.05[i] = NA
      b0.isw1.sd.05[i] = NA
      b1.isw1.05[i] = NA
      b1.isw1.sd.05[i] = NA
      # Coverage
      cover.isw1.05[i,1] = NA
      cover.isw1.05[i,2] = NA
      cover.isw1.05[i,3] = NA
      cover.isw1.05[i,4] = NA
      cover.isw1.05[i,5] = NA
      cover.isw1.05[i,6] = NA
      cover.isw1.05[i,7] = NA
      cover.isw1.05[i,8] = NA
    })

  # method 2 : rcpp's method (rcpp's weight + weight-in + dfsane)
  tryCatch({
    rcpp.fit = rcpp.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.rcpp.05[i] = rcpp.fit[1,1]
    b0.rcpp.sd.05[i] = rcpp.fit[1,2]
    b1.rcpp.05[i] = rcpp.fit[2,1]
    b1.rcpp.sd.05[i] = rcpp.fit[2,2]
    # Coverage
    cover.rcpp.05[i,1] = rcpp.fit[1,1]-1.96*rcpp.fit[1,2]
    cover.rcpp.05[i,2] = rcpp.fit[1,1]
    cover.rcpp.05[i,3] = rcpp.fit[1,1]+1.96*rcpp.fit[1,2]
    cover.rcpp.05[i,4] = ind(1.609438, cover.rcpp.05[i,1], cover.rcpp.05[i,3])
    cover.rcpp.05[i,5] = rcpp.fit[2,1]-1.96*rcpp.fit[2,2]
    cover.rcpp.05[i,6] = rcpp.fit[2,1]
    cover.rcpp.05[i,7] = rcpp.fit[2,1]+1.96*rcpp.fit[2,2]
    cover.rcpp.05[i,8] = ind(0.6931472, cover.rcpp.05[i,5], cover.rcpp.05[i,7])}
    , error=function(e){
      b0.rcpp.05[i] = NA
      b0.rcpp.sd.05[i] = NA
      b1.rcpp.05[i] = NA
      b1.rcpp.sd.05[i] = NA
      # Coverage
      cover.rcpp.05[i,1] = NA
      cover.rcpp.05[i,2] = NA
      cover.rcpp.05[i,3] = NA
      cover.rcpp.05[i,4] = NA
      cover.rcpp.05[i,5] = NA
      cover.rcpp.05[i,6] = NA
      cover.rcpp.05[i,7] = NA
      cover.rcpp.05[i,8] = NA
    })
}

# IS beta table
table0.isw1[4,1]<-mean(b0.isw1.05,na.rm=TRUE)
table0.isw1[4,2]<-mean(b0.isw1.sd.05,na.rm=TRUE)
table0.isw1[4,3]<-sd(b0.isw1.05,na.rm=TRUE)
table0.isw1[4,4]<-mean(cover.isw1.05[,4],na.rm=TRUE)
table0.isw1[4,5]<-mean(b1.isw1.05,na.rm=TRUE)
table0.isw1[4,6]<-mean(b1.isw1.sd.05,na.rm=TRUE)
table0.isw1[4,7]<-sd(b1.isw1.05,na.rm=TRUE)
table0.isw1[4,8]<-mean(cover.isw1.05[,8],na.rm=TRUE)

table0.rcpp[4,1]<-mean(b0.rcpp.05,na.rm=TRUE)
table0.rcpp[4,2]<-mean(b0.rcpp.sd.05,na.rm=TRUE)
table0.rcpp[4,3]<-sd(b0.rcpp.05,na.rm=TRUE)
table0.rcpp[4,4]<-mean(cover.rcpp.05[,4],na.rm=TRUE)
table0.rcpp[4,5]<-mean(b1.rcpp.05,na.rm=TRUE)
table0.rcpp[4,6]<-mean(b1.rcpp.sd.05,na.rm=TRUE)
table0.rcpp[4,7]<-sd(b1.rcpp.05,na.rm=TRUE)
table0.rcpp[4,8]<-mean(cover.rcpp.05[,8],na.rm=TRUE)

#### t_0=0 & c=70% ####
b0.isw1.07 = c()
b0.isw1.sd.07 = c()
b1.isw1.07 = c()
b1.isw1.sd.07 = c()
cover.isw1.07=matrix(NA,2000,8)
colnames(cover.isw1.07)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.rcpp.07 = c()
b0.rcpp.sd.07 = c()
b1.rcpp.07 = c()
b1.rcpp.sd.07 = c()
cover.rcpp.07=matrix(NA,2000,8)
colnames(cover.rcpp.07)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.7)
  # method 1 : Our method (ISMB + rcpp's weight + weight-in + dfsane)
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.isw1.07[i] = ismb.fit[1,1]
    b0.isw1.sd.07[i] = ismb.fit[1,2]
    b1.isw1.07[i] = ismb.fit[2,1]
    b1.isw1.sd.07[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.07[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.07[i,2] = ismb.fit[1,1]
    cover.isw1.07[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.07[i,4] = ind(1.609438, cover.isw1.07[i,1], cover.isw1.07[i,3])
    cover.isw1.07[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.07[i,6] = ismb.fit[2,1]
    cover.isw1.07[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.07[i,8] = ind(0.6931472, cover.isw1.07[i,5], cover.isw1.07[i,7])}
    , error=function(e){
      b0.isw1.07[i] = NA
      b0.isw1.sd.07[i] = NA
      b1.isw1.07[i] = NA
      b1.isw1.sd.07[i] = NA
      # Coverage
      cover.isw1.07[i,1] = NA
      cover.isw1.07[i,2] = NA
      cover.isw1.07[i,3] = NA
      cover.isw1.07[i,4] = NA
      cover.isw1.07[i,5] = NA
      cover.isw1.07[i,6] = NA
      cover.isw1.07[i,7] = NA
      cover.isw1.07[i,8] = NA
    })

  # method 2 : rcpp's method (rcpp's weight + weight-in + dfsane)
  tryCatch({
    rcpp.fit = rcpp.est(a[,3], 1, a[,4], a[,5], 0, 0.5, 200)
    b0.rcpp.07[i] = rcpp.fit[1,1]
    b0.rcpp.sd.07[i] = rcpp.fit[1,2]
    b1.rcpp.07[i] = rcpp.fit[2,1]
    b1.rcpp.sd.07[i] = rcpp.fit[2,2]
    # Coverage
    cover.rcpp.07[i,1] = rcpp.fit[1,1]-1.96*rcpp.fit[1,2]
    cover.rcpp.07[i,2] = rcpp.fit[1,1]
    cover.rcpp.07[i,3] = rcpp.fit[1,1]+1.96*rcpp.fit[1,2]
    cover.rcpp.07[i,4] = ind(1.609438, cover.rcpp.07[i,1], cover.rcpp.07[i,3])
    cover.rcpp.07[i,5] = rcpp.fit[2,1]-1.96*rcpp.fit[2,2]
    cover.rcpp.07[i,6] = rcpp.fit[2,1]
    cover.rcpp.07[i,7] = rcpp.fit[2,1]+1.96*rcpp.fit[2,2]
    cover.rcpp.07[i,8] = ind(0.6931472, cover.rcpp.07[i,5], cover.rcpp.07[i,7])}
    , error=function(e){
      b0.rcpp.07[i] = NA
      b0.rcpp.sd.07[i] = NA
      b1.rcpp.07[i] = NA
      b1.rcpp.sd.07[i] = NA
      # Coverage
      cover.rcpp.07[i,1] = NA
      cover.rcpp.07[i,2] = NA
      cover.rcpp.07[i,3] = NA
      cover.rcpp.07[i,4] = NA
      cover.rcpp.07[i,5] = NA
      cover.rcpp.07[i,6] = NA
      cover.rcpp.07[i,7] = NA
      cover.rcpp.07[i,8] = NA
    })
}
# IS beta table
table0.isw1[5,1]<-mean(b0.isw1.07,na.rm=TRUE)
table0.isw1[5,2]<-mean(b0.isw1.sd.07,na.rm=TRUE)
table0.isw1[5,3]<-sd(b0.isw1.07,na.rm=TRUE)
table0.isw1[5,4]<-mean(cover.isw1.07[,4],na.rm=TRUE)
table0.isw1[5,5]<-mean(b1.isw1.07,na.rm=TRUE)
table0.isw1[5,6]<-mean(b1.isw1.sd.07,na.rm=TRUE)
table0.isw1[5,7]<-sd(b1.isw1.07,na.rm=TRUE)
table0.isw1[5,8]<-mean(cover.isw1.07[,8],na.rm=TRUE)

table0.rcpp[5,1]<-mean(b0.rcpp.07,na.rm=TRUE)
table0.rcpp[5,2]<-mean(b0.rcpp.sd.07,na.rm=TRUE)
table0.rcpp[5,3]<-sd(b0.rcpp.07,na.rm=TRUE)
table0.rcpp[5,4]<-mean(cover.rcpp.07[,4],na.rm=TRUE)
table0.rcpp[5,5]<-mean(b1.rcpp.07,na.rm=TRUE)
table0.rcpp[5,6]<-mean(b1.rcpp.sd.07,na.rm=TRUE)
table0.rcpp[5,7]<-sd(b1.rcpp.07,na.rm=TRUE)
table0.rcpp[5,8]<-mean(cover.rcpp.07[,8],na.rm=TRUE)

#### censoring point at t_0=1 ####
c.0=5000000
c.1=70.39
c.3=24.35
c.5=14.07
c.7=8.49
#### t_0=1 & c=0% ####
b0.isw1.10 = c()
b0.isw1.sd.10 = c()
b1.isw1.10 = c()
b1.isw1.sd.10 = c()
cover.isw1.10=matrix(NA,2000,8)
colnames(cover.isw1.10)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.rcpp.10 = c()
b0.rcpp.sd.10 = c()
b1.rcpp.10 = c()
b1.rcpp.sd.10 = c()
cover.rcpp.10=matrix(NA,2000,8)
colnames(cover.rcpp.10)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.0)
  # method 1 : Our method (ISMB + rcpp's weight + weight-in + dfsane)
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw1.10[i] = ismb.fit[1,1]
    b0.isw1.sd.10[i] = ismb.fit[1,2]
    b1.isw1.10[i] = ismb.fit[2,1]
    b1.isw1.sd.10[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.10[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.10[i,2] = ismb.fit[1,1]
    cover.isw1.10[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.10[i,4] = ind(1.410748, cover.isw1.10[i,1], cover.isw1.10[i,3])
    cover.isw1.10[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.10[i,6] = ismb.fit[2,1]
    cover.isw1.10[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.10[i,8] = ind(0.7974189, cover.isw1.10[i,5], cover.isw1.10[i,7])}
    , error=function(e){
      b0.isw1.10[i] = NA
      b0.isw1.sd.10[i] = NA
      b1.isw1.10[i] = NA
      b1.isw1.sd.10[i] = NA
      # Coverage
      cover.isw1.10[i,1] = NA
      cover.isw1.10[i,2] = NA
      cover.isw1.10[i,3] = NA
      cover.isw1.10[i,4] = NA
      cover.isw1.10[i,5] = NA
      cover.isw1.10[i,6] = NA
      cover.isw1.10[i,7] = NA
      cover.isw1.10[i,8] = NA
    })

  # method 2 : rcpp's method (rcpp's weight + weight-in + dfsane)
  tryCatch({
    rcpp.fit = rcpp.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.rcpp.10[i] = rcpp.fit[1,1]
    b0.rcpp.sd.10[i] = rcpp.fit[1,2]
    b1.rcpp.10[i] = rcpp.fit[2,1]
    b1.rcpp.sd.10[i] = rcpp.fit[2,2]
    # Coverage
    cover.rcpp.10[i,1] = rcpp.fit[1,1]-1.96*rcpp.fit[1,2]
    cover.rcpp.10[i,2] = rcpp.fit[1,1]
    cover.rcpp.10[i,3] = rcpp.fit[1,1]+1.96*rcpp.fit[1,2]
    cover.rcpp.10[i,4] = ind(1.410748, cover.rcpp.10[i,1], cover.rcpp.10[i,3])
    cover.rcpp.10[i,5] = rcpp.fit[2,1]-1.96*rcpp.fit[2,2]
    cover.rcpp.10[i,6] = rcpp.fit[2,1]
    cover.rcpp.10[i,7] = rcpp.fit[2,1]+1.96*rcpp.fit[2,2]
    cover.rcpp.10[i,8] = ind(0.7974189, cover.rcpp.10[i,5], cover.rcpp.10[i,7])}
    , error=function(e){
      b0.rcpp.10[i] = NA
      b0.rcpp.sd.10[i] = NA
      b1.rcpp.10[i] = NA
      b1.rcpp.sd.10[i] = NA
      # Coverage
      cover.rcpp.10[i,1] = NA
      cover.rcpp.10[i,2] = NA
      cover.rcpp.10[i,3] = NA
      cover.rcpp.10[i,4] = NA
      cover.rcpp.10[i,5] = NA
      cover.rcpp.10[i,6] = NA
      cover.rcpp.10[i,7] = NA
      cover.rcpp.10[i,8] = NA
    })
}

# IS beta table
table1.isw1[1,1]<-mean(b0.isw1.10,na.rm=TRUE)
table1.isw1[1,2]<-mean(b0.isw1.sd.10,na.rm=TRUE)
table1.isw1[1,3]<-sd(b0.isw1.10,na.rm=TRUE)
table1.isw1[1,4]<-mean(cover.isw1.10[,4],na.rm=TRUE)
table1.isw1[1,5]<-mean(b1.isw1.10,na.rm=TRUE)
table1.isw1[1,6]<-mean(b1.isw1.sd.10,na.rm=TRUE)
table1.isw1[1,7]<-sd(b1.isw1.10,na.rm=TRUE)
table1.isw1[1,8]<-mean(cover.isw1.10[,8],na.rm=TRUE)

table1.rcpp[1,1]<-mean(b0.rcpp.10,na.rm=TRUE)
table1.rcpp[1,2]<-mean(b0.rcpp.sd.10,na.rm=TRUE)
table1.rcpp[1,3]<-sd(b0.rcpp.10,na.rm=TRUE)
table1.rcpp[1,4]<-mean(cover.rcpp.10[,4],na.rm=TRUE)
table1.rcpp[1,5]<-mean(b1.rcpp.10,na.rm=TRUE)
table1.rcpp[1,6]<-mean(b1.rcpp.sd.10,na.rm=TRUE)
table1.rcpp[1,7]<-sd(b1.rcpp.10,na.rm=TRUE)
table1.rcpp[1,8]<-mean(cover.rcpp.10[,8],na.rm=TRUE)

#### t_0=1 & c=10% ####
b0.isw1.11 = c()
b0.isw1.sd.11 = c()
b1.isw1.11 = c()
b1.isw1.sd.11 = c()
cover.isw1.11=matrix(NA,2000,8)
colnames(cover.isw1.11)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.rcpp.11 = c()
b0.rcpp.sd.11 = c()
b1.rcpp.11 = c()
b1.rcpp.sd.11 = c()
cover.rcpp.11=matrix(NA,2000,8)
colnames(cover.rcpp.11)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.1)
  # method 1 : Our method (ISMB + rcpp's weight + weight-in + dfsane)
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw1.11[i] = ismb.fit[1,1]
    b0.isw1.sd.11[i] = ismb.fit[1,2]
    b1.isw1.11[i] = ismb.fit[2,1]
    b1.isw1.sd.11[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.11[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.11[i,2] = ismb.fit[1,1]
    cover.isw1.11[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.11[i,4] = ind(1.410748, cover.isw1.11[i,1], cover.isw1.11[i,3])
    cover.isw1.11[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.11[i,6] = ismb.fit[2,1]
    cover.isw1.11[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.11[i,8] = ind(0.7974189, cover.isw1.11[i,5], cover.isw1.11[i,7])}
    , error=function(e){
      b0.isw1.11[i] = NA
      b0.isw1.sd.11[i] = NA
      b1.isw1.11[i] = NA
      b1.isw1.sd.11[i] = NA
      # Coverage
      cover.isw1.11[i,1] = NA
      cover.isw1.11[i,2] = NA
      cover.isw1.11[i,3] = NA
      cover.isw1.11[i,4] = NA
      cover.isw1.11[i,5] = NA
      cover.isw1.11[i,6] = NA
      cover.isw1.11[i,7] = NA
      cover.isw1.11[i,8] = NA
    })

  # method 2 : rcpp's method (rcpp's weight + weight-in + dfsane)
  tryCatch({
    rcpp.fit = rcpp.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.rcpp.11[i] = rcpp.fit[1,1]
    b0.rcpp.sd.11[i] = rcpp.fit[1,2]
    b1.rcpp.11[i] = rcpp.fit[2,1]
    b1.rcpp.sd.11[i] = rcpp.fit[2,2]
    # Coverage
    cover.rcpp.11[i,1] = rcpp.fit[1,1]-1.96*rcpp.fit[1,2]
    cover.rcpp.11[i,2] = rcpp.fit[1,1]
    cover.rcpp.11[i,3] = rcpp.fit[1,1]+1.96*rcpp.fit[1,2]
    cover.rcpp.11[i,4] = ind(1.410748, cover.rcpp.11[i,1], cover.rcpp.11[i,3])
    cover.rcpp.11[i,5] = rcpp.fit[2,1]-1.96*rcpp.fit[2,2]
    cover.rcpp.11[i,6] = rcpp.fit[2,1]
    cover.rcpp.11[i,7] = rcpp.fit[2,1]+1.96*rcpp.fit[2,2]
    cover.rcpp.11[i,8] = ind(0.7974189, cover.rcpp.11[i,5], cover.rcpp.11[i,7])}
    , error=function(e){
      b0.rcpp.11[i] = NA
      b0.rcpp.sd.11[i] = NA
      b1.rcpp.11[i] = NA
      b1.rcpp.sd.11[i] = NA
      # Coverage
      cover.rcpp.11[i,1] = NA
      cover.rcpp.11[i,2] = NA
      cover.rcpp.11[i,3] = NA
      cover.rcpp.11[i,4] = NA
      cover.rcpp.11[i,5] = NA
      cover.rcpp.11[i,6] = NA
      cover.rcpp.11[i,7] = NA
      cover.rcpp.11[i,8] = NA
    })
}

# IS beta table
table1.isw1[2,1]<-mean(b0.isw1.11,na.rm=TRUE)
table1.isw1[2,2]<-mean(b0.isw1.sd.11,na.rm=TRUE)
table1.isw1[2,3]<-sd(b0.isw1.11,na.rm=TRUE)
table1.isw1[2,4]<-mean(cover.isw1.11[,4],na.rm=TRUE)
table1.isw1[2,5]<-mean(b1.isw1.11,na.rm=TRUE)
table1.isw1[2,6]<-mean(b1.isw1.sd.11,na.rm=TRUE)
table1.isw1[2,7]<-sd(b1.isw1.11,na.rm=TRUE)
table1.isw1[2,8]<-mean(cover.isw1.11[,8],na.rm=TRUE)

table1.rcpp[2,1]<-mean(b0.rcpp.11,na.rm=TRUE)
table1.rcpp[2,2]<-mean(b0.rcpp.sd.11,na.rm=TRUE)
table1.rcpp[2,3]<-sd(b0.rcpp.11,na.rm=TRUE)
table1.rcpp[2,4]<-mean(cover.rcpp.11[,4],na.rm=TRUE)
table1.rcpp[2,5]<-mean(b1.rcpp.11,na.rm=TRUE)
table1.rcpp[2,6]<-mean(b1.rcpp.sd.11,na.rm=TRUE)
table1.rcpp[2,7]<-sd(b1.rcpp.11,na.rm=TRUE)
table1.rcpp[2,8]<-mean(cover.rcpp.11[,8],na.rm=TRUE)

#### t_0=1 & c=30% ####
b0.isw1.13 = c()
b0.isw1.sd.13 = c()
b1.isw1.13 = c()
b1.isw1.sd.13 = c()
cover.isw1.13=matrix(NA,2000,8)
colnames(cover.isw1.13)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.rcpp.13 = c()
b0.rcpp.sd.13 = c()
b1.rcpp.13 = c()
b1.rcpp.sd.13 = c()
cover.rcpp.13=matrix(NA,2000,8)
colnames(cover.rcpp.13)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.3)
  # method 1 : Our method (ISMB + rcpp's weight + weight-in + dfsane)
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw1.13[i] = ismb.fit[1,1]
    b0.isw1.sd.13[i] = ismb.fit[1,2]
    b1.isw1.13[i] = ismb.fit[2,1]
    b1.isw1.sd.13[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.13[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.13[i,2] = ismb.fit[1,1]
    cover.isw1.13[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.13[i,4] = ind(1.410748, cover.isw1.13[i,1], cover.isw1.13[i,3])
    cover.isw1.13[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.13[i,6] = ismb.fit[2,1]
    cover.isw1.13[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.13[i,8] = ind(0.7974189, cover.isw1.13[i,5], cover.isw1.13[i,7])}
    , error=function(e){
      b0.isw1.13[i] = NA
      b0.isw1.sd.13[i] = NA
      b1.isw1.13[i] = NA
      b1.isw1.sd.13[i] = NA
      # Coverage
      cover.isw1.13[i,1] = NA
      cover.isw1.13[i,2] = NA
      cover.isw1.13[i,3] = NA
      cover.isw1.13[i,4] = NA
      cover.isw1.13[i,5] = NA
      cover.isw1.13[i,6] = NA
      cover.isw1.13[i,7] = NA
      cover.isw1.13[i,8] = NA
    })

  # method 2 : rcpp's method (rcpp's weight + weight-in + dfsane)
  tryCatch({
    rcpp.fit = rcpp.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.rcpp.13[i] = rcpp.fit[1,1]
    b0.rcpp.sd.13[i] = rcpp.fit[1,2]
    b1.rcpp.13[i] = rcpp.fit[2,1]
    b1.rcpp.sd.13[i] = rcpp.fit[2,2]
    # Coverage
    cover.rcpp.13[i,1] = rcpp.fit[1,1]-1.96*rcpp.fit[1,2]
    cover.rcpp.13[i,2] = rcpp.fit[1,1]
    cover.rcpp.13[i,3] = rcpp.fit[1,1]+1.96*rcpp.fit[1,2]
    cover.rcpp.13[i,4] = ind(1.410748, cover.rcpp.13[i,1], cover.rcpp.13[i,3])
    cover.rcpp.13[i,5] = rcpp.fit[2,1]-1.96*rcpp.fit[2,2]
    cover.rcpp.13[i,6] = rcpp.fit[2,1]
    cover.rcpp.13[i,7] = rcpp.fit[2,1]+1.96*rcpp.fit[2,2]
    cover.rcpp.13[i,8] = ind(0.7974189, cover.rcpp.13[i,5], cover.rcpp.13[i,7])}
    , error=function(e){
      b0.rcpp.13[i] = NA
      b0.rcpp.sd.13[i] = NA
      b1.rcpp.13[i] = NA
      b1.rcpp.sd.13[i] = NA
      # Coverage
      cover.rcpp.13[i,1] = NA
      cover.rcpp.13[i,2] = NA
      cover.rcpp.13[i,3] = NA
      cover.rcpp.13[i,4] = NA
      cover.rcpp.13[i,5] = NA
      cover.rcpp.13[i,6] = NA
      cover.rcpp.13[i,7] = NA
      cover.rcpp.13[i,8] = NA
    })
}
# IS beta table
table1.isw1[3,1]<-mean(b0.isw1.13,na.rm=TRUE)
table1.isw1[3,2]<-mean(b0.isw1.sd.13,na.rm=TRUE)
table1.isw1[3,3]<-sd(b0.isw1.13,na.rm=TRUE)
table1.isw1[3,4]<-mean(cover.isw1.13[,4],na.rm=TRUE)
table1.isw1[3,5]<-mean(b1.isw1.13,na.rm=TRUE)
table1.isw1[3,6]<-mean(b1.isw1.sd.13,na.rm=TRUE)
table1.isw1[3,7]<-sd(b1.isw1.13,na.rm=TRUE)
table1.isw1[3,8]<-mean(cover.isw1.13[,8],na.rm=TRUE)

table1.rcpp[3,1]<-mean(b0.rcpp.13,na.rm=TRUE)
table1.rcpp[3,2]<-mean(b0.rcpp.sd.13,na.rm=TRUE)
table1.rcpp[3,3]<-sd(b0.rcpp.13,na.rm=TRUE)
table1.rcpp[3,4]<-mean(cover.rcpp.13[,4],na.rm=TRUE)
table1.rcpp[3,5]<-mean(b1.rcpp.13,na.rm=TRUE)
table1.rcpp[3,6]<-mean(b1.rcpp.sd.13,na.rm=TRUE)
table1.rcpp[3,7]<-sd(b1.rcpp.13,na.rm=TRUE)
table1.rcpp[3,8]<-mean(cover.rcpp.13[,8],na.rm=TRUE)

#### t_0=1 & c=50% ####
b0.isw1.15 = c()
b0.isw1.sd.15 = c()
b1.isw1.15 = c()
b1.isw1.sd.15 = c()
cover.isw1.15=matrix(NA,2000,8)
colnames(cover.isw1.15)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.rcpp.15 = c()
b0.rcpp.sd.15 = c()
b1.rcpp.15 = c()
b1.rcpp.sd.15 = c()
cover.rcpp.15=matrix(NA,2000,8)
colnames(cover.rcpp.15)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.5)
  # method 1 : Our method (ISMB + rcpp's weight + weight-in + dfsane)
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw1.15[i] = ismb.fit[1,1]
    b0.isw1.sd.15[i] = ismb.fit[1,2]
    b1.isw1.15[i] = ismb.fit[2,1]
    b1.isw1.sd.15[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.15[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.15[i,2] = ismb.fit[1,1]
    cover.isw1.15[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.15[i,4] = ind(1.410748, cover.isw1.15[i,1], cover.isw1.15[i,3])
    cover.isw1.15[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.15[i,6] = ismb.fit[2,1]
    cover.isw1.15[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.15[i,8] = ind(0.7974189, cover.isw1.15[i,5], cover.isw1.15[i,7])}
    , error=function(e){
      b0.isw1.15[i] = NA
      b0.isw1.sd.15[i] = NA
      b1.isw1.15[i] = NA
      b1.isw1.sd.15[i] = NA
      # Coverage
      cover.isw1.15[i,1] = NA
      cover.isw1.15[i,2] = NA
      cover.isw1.15[i,3] = NA
      cover.isw1.15[i,4] = NA
      cover.isw1.15[i,5] = NA
      cover.isw1.15[i,6] = NA
      cover.isw1.15[i,7] = NA
      cover.isw1.15[i,8] = NA
    })

  # method 2 : rcpp's method (rcpp's weight + weight-in + dfsane)
  tryCatch({
    rcpp.fit = rcpp.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.rcpp.15[i] = rcpp.fit[1,1]
    b0.rcpp.sd.15[i] = rcpp.fit[1,2]
    b1.rcpp.15[i] = rcpp.fit[2,1]
    b1.rcpp.sd.15[i] = rcpp.fit[2,2]
    # Coverage
    cover.rcpp.15[i,1] = rcpp.fit[1,1]-1.96*rcpp.fit[1,2]
    cover.rcpp.15[i,2] = rcpp.fit[1,1]
    cover.rcpp.15[i,3] = rcpp.fit[1,1]+1.96*rcpp.fit[1,2]
    cover.rcpp.15[i,4] = ind(1.410748, cover.rcpp.15[i,1], cover.rcpp.15[i,3])
    cover.rcpp.15[i,5] = rcpp.fit[2,1]-1.96*rcpp.fit[2,2]
    cover.rcpp.15[i,6] = rcpp.fit[2,1]
    cover.rcpp.15[i,7] = rcpp.fit[2,1]+1.96*rcpp.fit[2,2]
    cover.rcpp.15[i,8] = ind(0.7974189, cover.rcpp.15[i,5], cover.rcpp.15[i,7])}
    , error=function(e){
      b0.rcpp.15[i] = NA
      b0.rcpp.sd.15[i] = NA
      b1.rcpp.15[i] = NA
      b1.rcpp.sd.15[i] = NA
      # Coverage
      cover.rcpp.15[i,1] = NA
      cover.rcpp.15[i,2] = NA
      cover.rcpp.15[i,3] = NA
      cover.rcpp.15[i,4] = NA
      cover.rcpp.15[i,5] = NA
      cover.rcpp.15[i,6] = NA
      cover.rcpp.15[i,7] = NA
      cover.rcpp.15[i,8] = NA
    })
}

# IS beta table
table1.isw1[4,1]<-mean(b0.isw1.15,na.rm=TRUE)
table1.isw1[4,2]<-mean(b0.isw1.sd.15,na.rm=TRUE)
table1.isw1[4,3]<-sd(b0.isw1.15,na.rm=TRUE)
table1.isw1[4,4]<-mean(cover.isw1.15[,4],na.rm=TRUE)
table1.isw1[4,5]<-mean(b1.isw1.15,na.rm=TRUE)
table1.isw1[4,6]<-mean(b1.isw1.sd.15,na.rm=TRUE)
table1.isw1[4,7]<-sd(b1.isw1.15,na.rm=TRUE)
table1.isw1[4,8]<-mean(cover.isw1.15[,8],na.rm=TRUE)

table1.rcpp[4,1]<-mean(b0.rcpp.15,na.rm=TRUE)
table1.rcpp[4,2]<-mean(b0.rcpp.sd.15,na.rm=TRUE)
table1.rcpp[4,3]<-sd(b0.rcpp.15,na.rm=TRUE)
table1.rcpp[4,4]<-mean(cover.rcpp.15[,4],na.rm=TRUE)
table1.rcpp[4,5]<-mean(b1.rcpp.15,na.rm=TRUE)
table1.rcpp[4,6]<-mean(b1.rcpp.sd.15,na.rm=TRUE)
table1.rcpp[4,7]<-sd(b1.rcpp.15,na.rm=TRUE)
table1.rcpp[4,8]<-mean(cover.rcpp.15[,8],na.rm=TRUE)

#### t_0=1 & c=70% ####
b0.isw1.17 = c()
b0.isw1.sd.17 = c()
b1.isw1.17 = c()
b1.isw1.sd.17 = c()
cover.isw1.17=matrix(NA,2000,8)
colnames(cover.isw1.17)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.rcpp.17 = c()
b0.rcpp.sd.17 = c()
b1.rcpp.17 = c()
b1.rcpp.sd.17 = c()
cover.rcpp.17=matrix(NA,2000,8)
colnames(cover.rcpp.17)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.7)
  # method 1 : Our method (ISMB + rcpp's weight + weight-in + dfsane)
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.isw1.17[i] = ismb.fit[1,1]
    b0.isw1.sd.17[i] = ismb.fit[1,2]
    b1.isw1.17[i] = ismb.fit[2,1]
    b1.isw1.sd.17[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.17[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.17[i,2] = ismb.fit[1,1]
    cover.isw1.17[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.17[i,4] = ind(1.410748, cover.isw1.17[i,1], cover.isw1.17[i,3])
    cover.isw1.17[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.17[i,6] = ismb.fit[2,1]
    cover.isw1.17[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.17[i,8] = ind(0.7974189, cover.isw1.17[i,5], cover.isw1.17[i,7])}
    , error=function(e){
      b0.isw1.17[i] = NA
      b0.isw1.sd.17[i] = NA
      b1.isw1.17[i] = NA
      b1.isw1.sd.17[i] = NA
      # Coverage
      cover.isw1.17[i,1] = NA
      cover.isw1.17[i,2] = NA
      cover.isw1.17[i,3] = NA
      cover.isw1.17[i,4] = NA
      cover.isw1.17[i,5] = NA
      cover.isw1.17[i,6] = NA
      cover.isw1.17[i,7] = NA
      cover.isw1.17[i,8] = NA
    })

  # method 2 : rcpp's method (rcpp's weight + weight-in + dfsane)
  tryCatch({
    rcpp.fit = rcpp.est(a[,3], 1, a[,4], a[,5], 1, 0.5, 200)
    b0.rcpp.17[i] = rcpp.fit[1,1]
    b0.rcpp.sd.17[i] = rcpp.fit[1,2]
    b1.rcpp.17[i] = rcpp.fit[2,1]
    b1.rcpp.sd.17[i] = rcpp.fit[2,2]
    # Coverage
    cover.rcpp.17[i,1] = rcpp.fit[1,1]-1.96*rcpp.fit[1,2]
    cover.rcpp.17[i,2] = rcpp.fit[1,1]
    cover.rcpp.17[i,3] = rcpp.fit[1,1]+1.96*rcpp.fit[1,2]
    cover.rcpp.17[i,4] = ind(1.410748, cover.rcpp.17[i,1], cover.rcpp.17[i,3])
    cover.rcpp.17[i,5] = rcpp.fit[2,1]-1.96*rcpp.fit[2,2]
    cover.rcpp.17[i,6] = rcpp.fit[2,1]
    cover.rcpp.17[i,7] = rcpp.fit[2,1]+1.96*rcpp.fit[2,2]
    cover.rcpp.17[i,8] = ind(0.7974189, cover.rcpp.17[i,5], cover.rcpp.17[i,7])}
    , error=function(e){
      b0.rcpp.17[i] = NA
      b0.rcpp.sd.17[i] = NA
      b1.rcpp.17[i] = NA
      b1.rcpp.sd.17[i] = NA
      # Coverage
      cover.rcpp.17[i,1] = NA
      cover.rcpp.17[i,2] = NA
      cover.rcpp.17[i,3] = NA
      cover.rcpp.17[i,4] = NA
      cover.rcpp.17[i,5] = NA
      cover.rcpp.17[i,6] = NA
      cover.rcpp.17[i,7] = NA
      cover.rcpp.17[i,8] = NA
    })
}
# IS beta table
table1.isw1[5,1]<-mean(b0.isw1.17,na.rm=TRUE)
table1.isw1[5,2]<-mean(b0.isw1.sd.17,na.rm=TRUE)
table1.isw1[5,3]<-sd(b0.isw1.17,na.rm=TRUE)
table1.isw1[5,4]<-mean(cover.isw1.17[,4],na.rm=TRUE)
table1.isw1[5,5]<-mean(b1.isw1.17,na.rm=TRUE)
table1.isw1[5,6]<-mean(b1.isw1.sd.17,na.rm=TRUE)
table1.isw1[5,7]<-sd(b1.isw1.17,na.rm=TRUE)
table1.isw1[5,8]<-mean(cover.isw1.17[,8],na.rm=TRUE)

table1.rcpp[5,1]<-mean(b0.rcpp.17,na.rm=TRUE)
table1.rcpp[5,2]<-mean(b0.rcpp.sd.17,na.rm=TRUE)
table1.rcpp[5,3]<-sd(b0.rcpp.17,na.rm=TRUE)
table1.rcpp[5,4]<-mean(cover.rcpp.17[,4],na.rm=TRUE)
table1.rcpp[5,5]<-mean(b1.rcpp.17,na.rm=TRUE)
table1.rcpp[5,6]<-mean(b1.rcpp.sd.17,na.rm=TRUE)
table1.rcpp[5,7]<-sd(b1.rcpp.17,na.rm=TRUE)
table1.rcpp[5,8]<-mean(cover.rcpp.17[,8],na.rm=TRUE)

#### censoring point at t_0=2 ####
c.0=5000000
c.1=64.86
c.3=23.34
c.5=13.62
c.7=8.36
#### t_0=2 & c=0% ####
b0.isw1.20 = c()
b0.isw1.sd.20 = c()
b1.isw1.20 = c()
b1.isw1.sd.20 = c()
cover.isw1.20=matrix(NA,2000,8)
colnames(cover.isw1.20)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.rcpp.20 = c()
b0.rcpp.sd.20 = c()
b1.rcpp.20 = c()
b1.rcpp.sd.20 = c()
cover.rcpp.20=matrix(NA,2000,8)
colnames(cover.rcpp.20)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.0)
  # method 1 : Our method (ISMB + rcpp's weight + weight-in + dfsane)
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw1.20[i] = ismb.fit[1,1]
    b0.isw1.sd.20[i] = ismb.fit[1,2]
    b1.isw1.20[i] = ismb.fit[2,1]
    b1.isw1.sd.20[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.20[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.20[i,2] = ismb.fit[1,1]
    cover.isw1.20[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.20[i,4] = ind(1.219403, cover.isw1.20[i,1], cover.isw1.20[i,3])
    cover.isw1.20[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.20[i,6] = ismb.fit[2,1]
    cover.isw1.20[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.20[i,8] = ind(0.9070615, cover.isw1.20[i,5], cover.isw1.20[i,7])}
    , error=function(e){
      b0.isw1.20[i] = NA
      b0.isw1.sd.20[i] = NA
      b1.isw1.20[i] = NA
      b1.isw1.sd.20[i] = NA
      # Coverage
      cover.isw1.20[i,1] = NA
      cover.isw1.20[i,2] = NA
      cover.isw1.20[i,3] = NA
      cover.isw1.20[i,4] = NA
      cover.isw1.20[i,5] = NA
      cover.isw1.20[i,6] = NA
      cover.isw1.20[i,7] = NA
      cover.isw1.20[i,8] = NA
    })

    # method 2 : rcpp's method (rcpp's weight + weight-in + dfsane)
  tryCatch({
    rcpp.fit = rcpp.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.rcpp.20[i] = rcpp.fit[1,1]
    b0.rcpp.sd.20[i] = rcpp.fit[1,2]
    b1.rcpp.20[i] = rcpp.fit[2,1]
    b1.rcpp.sd.20[i] = rcpp.fit[2,2]
    # Coverage
    cover.rcpp.20[i,1] = rcpp.fit[1,1]-1.96*rcpp.fit[1,2]
    cover.rcpp.20[i,2] = rcpp.fit[1,1]
    cover.rcpp.20[i,3] = rcpp.fit[1,1]+1.96*rcpp.fit[1,2]
    cover.rcpp.20[i,4] = ind(1.219403, cover.rcpp.20[i,1], cover.rcpp.20[i,3])
    cover.rcpp.20[i,5] = rcpp.fit[2,1]-1.96*rcpp.fit[2,2]
    cover.rcpp.20[i,6] = rcpp.fit[2,1]
    cover.rcpp.20[i,7] = rcpp.fit[2,1]+1.96*rcpp.fit[2,2]
    cover.rcpp.20[i,8] = ind(0.9070615, cover.rcpp.20[i,5], cover.rcpp.20[i,7])}
    , error=function(e){
      b0.rcpp.20[i] = NA
      b0.rcpp.sd.20[i] = NA
      b1.rcpp.20[i] = NA
      b1.rcpp.sd.20[i] = NA
      # Coverage
      cover.rcpp.20[i,1] = NA
      cover.rcpp.20[i,2] = NA
      cover.rcpp.20[i,3] = NA
      cover.rcpp.20[i,4] = NA
      cover.rcpp.20[i,5] = NA
      cover.rcpp.20[i,6] = NA
      cover.rcpp.20[i,7] = NA
      cover.rcpp.20[i,8] = NA
    })
}

# IS beta table
table2.isw1[1,1]<-mean(b0.isw1.20,na.rm=TRUE)
table2.isw1[1,2]<-mean(b0.isw1.sd.20,na.rm=TRUE)
table2.isw1[1,3]<-sd(b0.isw1.20,na.rm=TRUE)
table2.isw1[1,4]<-mean(cover.isw1.20[,4],na.rm=TRUE)
table2.isw1[1,5]<-mean(b1.isw1.20,na.rm=TRUE)
table2.isw1[1,6]<-mean(b1.isw1.sd.20,na.rm=TRUE)
table2.isw1[1,7]<-sd(b1.isw1.20,na.rm=TRUE)
table2.isw1[1,8]<-mean(cover.isw1.20[,8],na.rm=TRUE)

table2.rcpp[1,1]<-mean(b0.rcpp.20,na.rm=TRUE)
table2.rcpp[1,2]<-mean(b0.rcpp.sd.20,na.rm=TRUE)
table2.rcpp[1,3]<-sd(b0.rcpp.20,na.rm=TRUE)
table2.rcpp[1,4]<-mean(cover.rcpp.20[,4],na.rm=TRUE)
table2.rcpp[1,5]<-mean(b1.rcpp.20,na.rm=TRUE)
table2.rcpp[1,6]<-mean(b1.rcpp.sd.20,na.rm=TRUE)
table2.rcpp[1,7]<-sd(b1.rcpp.20,na.rm=TRUE)
table2.rcpp[1,8]<-mean(cover.rcpp.20[,8],na.rm=TRUE)

#### t_0=2 & c=10% ####
b0.isw1.21 = c()
b0.isw1.sd.21 = c()
b1.isw1.21 = c()
b1.isw1.sd.21 = c()
cover.isw1.21=matrix(NA,2000,8)
colnames(cover.isw1.21)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.rcpp.21 = c()
b0.rcpp.sd.21 = c()
b1.rcpp.21 = c()
b1.rcpp.sd.21 = c()
cover.rcpp.21=matrix(NA,2000,8)
colnames(cover.rcpp.21)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.1)
  # method 1 : Our method (ISMB + rcpp's weight + weight-in + dfsane)
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw1.21[i] = ismb.fit[1,1]
    b0.isw1.sd.21[i] = ismb.fit[1,2]
    b1.isw1.21[i] = ismb.fit[2,1]
    b1.isw1.sd.21[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.21[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.21[i,2] = ismb.fit[1,1]
    cover.isw1.21[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.21[i,4] = ind(1.219403, cover.isw1.21[i,1], cover.isw1.21[i,3])
    cover.isw1.21[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.21[i,6] = ismb.fit[2,1]
    cover.isw1.21[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.21[i,8] = ind(0.9070615, cover.isw1.21[i,5], cover.isw1.21[i,7])}
    , error=function(e){
      b0.isw1.21[i] = NA
      b0.isw1.sd.21[i] = NA
      b1.isw1.21[i] = NA
      b1.isw1.sd.21[i] = NA
      # Coverage
      cover.isw1.21[i,1] = NA
      cover.isw1.21[i,2] = NA
      cover.isw1.21[i,3] = NA
      cover.isw1.21[i,4] = NA
      cover.isw1.21[i,5] = NA
      cover.isw1.21[i,6] = NA
      cover.isw1.21[i,7] = NA
      cover.isw1.21[i,8] = NA
    })

  # method 2 : rcpp's method (rcpp's weight + weight-in + dfsane)
  tryCatch({
    rcpp.fit = rcpp.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.rcpp.21[i] = rcpp.fit[1,1]
    b0.rcpp.sd.21[i] = rcpp.fit[1,2]
    b1.rcpp.21[i] = rcpp.fit[2,1]
    b1.rcpp.sd.21[i] = rcpp.fit[2,2]
    # Coverage
    cover.rcpp.21[i,1] = rcpp.fit[1,1]-1.96*rcpp.fit[1,2]
    cover.rcpp.21[i,2] = rcpp.fit[1,1]
    cover.rcpp.21[i,3] = rcpp.fit[1,1]+1.96*rcpp.fit[1,2]
    cover.rcpp.21[i,4] = ind(1.219403, cover.rcpp.21[i,1], cover.rcpp.21[i,3])
    cover.rcpp.21[i,5] = rcpp.fit[2,1]-1.96*rcpp.fit[2,2]
    cover.rcpp.21[i,6] = rcpp.fit[2,1]
    cover.rcpp.21[i,7] = rcpp.fit[2,1]+1.96*rcpp.fit[2,2]
    cover.rcpp.21[i,8] = ind(0.9070615, cover.rcpp.21[i,5], cover.rcpp.21[i,7])}
    , error=function(e){
      b0.rcpp.21[i] = NA
      b0.rcpp.sd.21[i] = NA
      b1.rcpp.21[i] = NA
      b1.rcpp.sd.21[i] = NA
      # Coverage
      cover.rcpp.21[i,1] = NA
      cover.rcpp.21[i,2] = NA
      cover.rcpp.21[i,3] = NA
      cover.rcpp.21[i,4] = NA
      cover.rcpp.21[i,5] = NA
      cover.rcpp.21[i,6] = NA
      cover.rcpp.21[i,7] = NA
      cover.rcpp.21[i,8] = NA
    })
}

# IS beta table
table2.isw1[2,1]<-mean(b0.isw1.21,na.rm=TRUE)
table2.isw1[2,2]<-mean(b0.isw1.sd.21,na.rm=TRUE)
table2.isw1[2,3]<-sd(b0.isw1.21,na.rm=TRUE)
table2.isw1[2,4]<-mean(cover.isw1.21[,4],na.rm=TRUE)
table2.isw1[2,5]<-mean(b1.isw1.21,na.rm=TRUE)
table2.isw1[2,6]<-mean(b1.isw1.sd.21,na.rm=TRUE)
table2.isw1[2,7]<-sd(b1.isw1.21,na.rm=TRUE)
table2.isw1[2,8]<-mean(cover.isw1.21[,8],na.rm=TRUE)

table2.rcpp[2,1]<-mean(b0.rcpp.21,na.rm=TRUE)
table2.rcpp[2,2]<-mean(b0.rcpp.sd.21,na.rm=TRUE)
table2.rcpp[2,3]<-sd(b0.rcpp.21,na.rm=TRUE)
table2.rcpp[2,4]<-mean(cover.rcpp.21[,4],na.rm=TRUE)
table2.rcpp[2,5]<-mean(b1.rcpp.21,na.rm=TRUE)
table2.rcpp[2,6]<-mean(b1.rcpp.sd.21,na.rm=TRUE)
table2.rcpp[2,7]<-sd(b1.rcpp.21,na.rm=TRUE)
table2.rcpp[2,8]<-mean(cover.rcpp.21[,8],na.rm=TRUE)

#### t_0=2 & c=30% ####
b0.isw1.23 = c()
b0.isw1.sd.23 = c()
b1.isw1.23 = c()
b1.isw1.sd.23 = c()
cover.isw1.23=matrix(NA,2000,8)
colnames(cover.isw1.23)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.rcpp.23 = c()
b0.rcpp.sd.23 = c()
b1.rcpp.23 = c()
b1.rcpp.sd.23 = c()
cover.rcpp.23=matrix(NA,2000,8)
colnames(cover.rcpp.23)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.3)
  # method 1 : Our method (ISMB + rcpp's weight + weight-in + dfsane)
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw1.23[i] = ismb.fit[1,1]
    b0.isw1.sd.23[i] = ismb.fit[1,2]
    b1.isw1.23[i] = ismb.fit[2,1]
    b1.isw1.sd.23[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.23[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.23[i,2] = ismb.fit[1,1]
    cover.isw1.23[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.23[i,4] = ind(1.219403, cover.isw1.23[i,1], cover.isw1.23[i,3])
    cover.isw1.23[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.23[i,6] = ismb.fit[2,1]
    cover.isw1.23[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.23[i,8] = ind(0.9070615, cover.isw1.23[i,5], cover.isw1.23[i,7])}
    , error=function(e){
      b0.isw1.23[i] = NA
      b0.isw1.sd.23[i] = NA
      b1.isw1.23[i] = NA
      b1.isw1.sd.23[i] = NA
      # Coverage
      cover.isw1.23[i,1] = NA
      cover.isw1.23[i,2] = NA
      cover.isw1.23[i,3] = NA
      cover.isw1.23[i,4] = NA
      cover.isw1.23[i,5] = NA
      cover.isw1.23[i,6] = NA
      cover.isw1.23[i,7] = NA
      cover.isw1.23[i,8] = NA
    })

  # method 2 : rcpp's method (rcpp's weight + weight-in + dfsane)
  tryCatch({
    rcpp.fit = rcpp.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.rcpp.23[i] = rcpp.fit[1,1]
    b0.rcpp.sd.23[i] = rcpp.fit[1,2]
    b1.rcpp.23[i] = rcpp.fit[2,1]
    b1.rcpp.sd.23[i] = rcpp.fit[2,2]
    # Coverage
    cover.rcpp.23[i,1] = rcpp.fit[1,1]-1.96*rcpp.fit[1,2]
    cover.rcpp.23[i,2] = rcpp.fit[1,1]
    cover.rcpp.23[i,3] = rcpp.fit[1,1]+1.96*rcpp.fit[1,2]
    cover.rcpp.23[i,4] = ind(1.219403, cover.rcpp.23[i,1], cover.rcpp.23[i,3])
    cover.rcpp.23[i,5] = rcpp.fit[2,1]-1.96*rcpp.fit[2,2]
    cover.rcpp.23[i,6] = rcpp.fit[2,1]
    cover.rcpp.23[i,7] = rcpp.fit[2,1]+1.96*rcpp.fit[2,2]
    cover.rcpp.23[i,8] = ind(0.9070615, cover.rcpp.23[i,5], cover.rcpp.23[i,7])}
    , error=function(e){
      b0.rcpp.23[i] = NA
      b0.rcpp.sd.23[i] = NA
      b1.rcpp.23[i] = NA
      b1.rcpp.sd.23[i] = NA
      # Coverage
      cover.rcpp.23[i,1] = NA
      cover.rcpp.23[i,2] = NA
      cover.rcpp.23[i,3] = NA
      cover.rcpp.23[i,4] = NA
      cover.rcpp.23[i,5] = NA
      cover.rcpp.23[i,6] = NA
      cover.rcpp.23[i,7] = NA
      cover.rcpp.23[i,8] = NA
    })
}

# IS beta table
table2.isw1[3,1]<-mean(b0.isw1.23,na.rm=TRUE)
table2.isw1[3,2]<-mean(b0.isw1.sd.23,na.rm=TRUE)
table2.isw1[3,3]<-sd(b0.isw1.23,na.rm=TRUE)
table2.isw1[3,4]<-mean(cover.isw1.23[,4],na.rm=TRUE)
table2.isw1[3,5]<-mean(b1.isw1.23,na.rm=TRUE)
table2.isw1[3,6]<-mean(b1.isw1.sd.23,na.rm=TRUE)
table2.isw1[3,7]<-sd(b1.isw1.23,na.rm=TRUE)
table2.isw1[3,8]<-mean(cover.isw1.23[,8],na.rm=TRUE)

table2.rcpp[3,1]<-mean(b0.rcpp.23,na.rm=TRUE)
table2.rcpp[3,2]<-mean(b0.rcpp.sd.23,na.rm=TRUE)
table2.rcpp[3,3]<-sd(b0.rcpp.23,na.rm=TRUE)
table2.rcpp[3,4]<-mean(cover.rcpp.23[,4],na.rm=TRUE)
table2.rcpp[3,5]<-mean(b1.rcpp.23,na.rm=TRUE)
table2.rcpp[3,6]<-mean(b1.rcpp.sd.23,na.rm=TRUE)
table2.rcpp[3,7]<-sd(b1.rcpp.23,na.rm=TRUE)
table2.rcpp[3,8]<-mean(cover.rcpp.23[,8],na.rm=TRUE)

#### t_0=2 & c=50% ####
b0.isw1.25 = c()
b0.isw1.sd.25 = c()
b1.isw1.25 = c()
b1.isw1.sd.25 = c()
cover.isw1.25=matrix(NA,2000,8)
colnames(cover.isw1.25)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.rcpp.25 = c()
b0.rcpp.sd.25 = c()
b1.rcpp.25 = c()
b1.rcpp.sd.25 = c()
cover.rcpp.25=matrix(NA,2000,8)
colnames(cover.rcpp.25)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.5)
    # method 1 : Our method (ISMB + rcpp's weight + weight-in + dfsane)
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw1.25[i] = ismb.fit[1,1]
    b0.isw1.sd.25[i] = ismb.fit[1,2]
    b1.isw1.25[i] = ismb.fit[2,1]
    b1.isw1.sd.25[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.25[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.25[i,2] = ismb.fit[1,1]
    cover.isw1.25[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.25[i,4] = ind(1.219403, cover.isw1.25[i,1], cover.isw1.25[i,3])
    cover.isw1.25[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.25[i,6] = ismb.fit[2,1]
    cover.isw1.25[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.25[i,8] = ind(0.9070615, cover.isw1.25[i,5], cover.isw1.25[i,7])}
    , error=function(e){
      b0.isw1.25[i] = NA
      b0.isw1.sd.25[i] = NA
      b1.isw1.25[i] = NA
      b1.isw1.sd.25[i] = NA
      # Coverage
      cover.isw1.25[i,1] = NA
      cover.isw1.25[i,2] = NA
      cover.isw1.25[i,3] = NA
      cover.isw1.25[i,4] = NA
      cover.isw1.25[i,5] = NA
      cover.isw1.25[i,6] = NA
      cover.isw1.25[i,7] = NA
      cover.isw1.25[i,8] = NA
    })

  # method 2 : rcpp's method (rcpp's weight + weight-in + dfsane)
  tryCatch({
    rcpp.fit = rcpp.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.rcpp.25[i] = rcpp.fit[1,1]
    b0.rcpp.sd.25[i] = rcpp.fit[1,2]
    b1.rcpp.25[i] = rcpp.fit[2,1]
    b1.rcpp.sd.25[i] = rcpp.fit[2,2]
    # Coverage
    cover.rcpp.25[i,1] = rcpp.fit[1,1]-1.96*rcpp.fit[1,2]
    cover.rcpp.25[i,2] = rcpp.fit[1,1]
    cover.rcpp.25[i,3] = rcpp.fit[1,1]+1.96*rcpp.fit[1,2]
    cover.rcpp.25[i,4] = ind(1.219403, cover.rcpp.25[i,1], cover.rcpp.25[i,3])
    cover.rcpp.25[i,5] = rcpp.fit[2,1]-1.96*rcpp.fit[2,2]
    cover.rcpp.25[i,6] = rcpp.fit[2,1]
    cover.rcpp.25[i,7] = rcpp.fit[2,1]+1.96*rcpp.fit[2,2]
    cover.rcpp.25[i,8] = ind(0.9070615, cover.rcpp.25[i,5], cover.rcpp.25[i,7])}
    , error=function(e){
      b0.rcpp.25[i] = NA
      b0.rcpp.sd.25[i] = NA
      b1.rcpp.25[i] = NA
      b1.rcpp.sd.25[i] = NA
      # Coverage
      cover.rcpp.25[i,1] = NA
      cover.rcpp.25[i,2] = NA
      cover.rcpp.25[i,3] = NA
      cover.rcpp.25[i,4] = NA
      cover.rcpp.25[i,5] = NA
      cover.rcpp.25[i,6] = NA
      cover.rcpp.25[i,7] = NA
      cover.rcpp.25[i,8] = NA
    })
}

# IS beta table
table2.isw1[4,1]<-mean(b0.isw1.25,na.rm=TRUE)
table2.isw1[4,2]<-mean(b0.isw1.sd.25,na.rm=TRUE)
table2.isw1[4,3]<-sd(b0.isw1.25,na.rm=TRUE)
table2.isw1[4,4]<-mean(cover.isw1.25[,4],na.rm=TRUE)
table2.isw1[4,5]<-mean(b1.isw1.25,na.rm=TRUE)
table2.isw1[4,6]<-mean(b1.isw1.sd.25,na.rm=TRUE)
table2.isw1[4,7]<-sd(b1.isw1.25,na.rm=TRUE)
table2.isw1[4,8]<-mean(cover.isw1.25[,8],na.rm=TRUE)

table2.rcpp[4,1]<-mean(b0.rcpp.25,na.rm=TRUE)
table2.rcpp[4,2]<-mean(b0.rcpp.sd.25,na.rm=TRUE)
table2.rcpp[4,3]<-sd(b0.rcpp.25,na.rm=TRUE)
table2.rcpp[4,4]<-mean(cover.rcpp.25[,4],na.rm=TRUE)
table2.rcpp[4,5]<-mean(b1.rcpp.25,na.rm=TRUE)
table2.rcpp[4,6]<-mean(b1.rcpp.sd.25,na.rm=TRUE)
table2.rcpp[4,7]<-sd(b1.rcpp.25,na.rm=TRUE)
table2.rcpp[4,8]<-mean(cover.rcpp.25[,8],na.rm=TRUE)

#### t_0=2 & c=70% ####
b0.isw1.27 = c()
b0.isw1.sd.27 = c()
b1.isw1.27 = c()
b1.isw1.sd.27 = c()
cover.isw1.27=matrix(NA,2000,8)
colnames(cover.isw1.27)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.rcpp.27 = c()
b0.rcpp.sd.27 = c()
b1.rcpp.27 = c()
b1.rcpp.sd.27 = c()
cover.rcpp.27=matrix(NA,2000,8)
colnames(cover.rcpp.27)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.7)
  # method 1 : Our method (ISMB + rcpp's weight + weight-in + dfsane)
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.isw1.27[i] = ismb.fit[1,1]
    b0.isw1.sd.27[i] = ismb.fit[1,2]
    b1.isw1.27[i] = ismb.fit[2,1]
    b1.isw1.sd.27[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.27[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.27[i,2] = ismb.fit[1,1]
    cover.isw1.27[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.27[i,4] = ind(1.219403, cover.isw1.27[i,1], cover.isw1.27[i,3])
    cover.isw1.27[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.27[i,6] = ismb.fit[2,1]
    cover.isw1.27[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.27[i,8] = ind(0.9070615, cover.isw1.27[i,5], cover.isw1.27[i,7])}
    , error=function(e){
      b0.isw1.27[i] = NA
      b0.isw1.sd.27[i] = NA
      b1.isw1.27[i] = NA
      b1.isw1.sd.27[i] = NA
      # Coverage
      cover.isw1.27[i,1] = NA
      cover.isw1.27[i,2] = NA
      cover.isw1.27[i,3] = NA
      cover.isw1.27[i,4] = NA
      cover.isw1.27[i,5] = NA
      cover.isw1.27[i,6] = NA
      cover.isw1.27[i,7] = NA
      cover.isw1.27[i,8] = NA
    })

  # method 2 : rcpp's method (rcpp's weight + weight-in + dfsane)
  tryCatch({
    rcpp.fit = rcpp.est(a[,3], 1, a[,4], a[,5], 2, 0.5, 200)
    b0.rcpp.27[i] = rcpp.fit[1,1]
    b0.rcpp.sd.27[i] = rcpp.fit[1,2]
    b1.rcpp.27[i] = rcpp.fit[2,1]
    b1.rcpp.sd.27[i] = rcpp.fit[2,2]
    # Coverage
    cover.rcpp.27[i,1] = rcpp.fit[1,1]-1.96*rcpp.fit[1,2]
    cover.rcpp.27[i,2] = rcpp.fit[1,1]
    cover.rcpp.27[i,3] = rcpp.fit[1,1]+1.96*rcpp.fit[1,2]
    cover.rcpp.27[i,4] = ind(1.219403, cover.rcpp.27[i,1], cover.rcpp.27[i,3])
    cover.rcpp.27[i,5] = rcpp.fit[2,1]-1.96*rcpp.fit[2,2]
    cover.rcpp.27[i,6] = rcpp.fit[2,1]
    cover.rcpp.27[i,7] = rcpp.fit[2,1]+1.96*rcpp.fit[2,2]
    cover.rcpp.27[i,8] = ind(0.9070615, cover.rcpp.27[i,5], cover.rcpp.27[i,7])}
    , error=function(e){
      b0.rcpp.27[i] = NA
      b0.rcpp.sd.27[i] = NA
      b1.rcpp.27[i] = NA
      b1.rcpp.sd.27[i] = NA
      # Coverage
      cover.rcpp.27[i,1] = NA
      cover.rcpp.27[i,2] = NA
      cover.rcpp.27[i,3] = NA
      cover.rcpp.27[i,4] = NA
      cover.rcpp.27[i,5] = NA
      cover.rcpp.27[i,6] = NA
      cover.rcpp.27[i,7] = NA
      cover.rcpp.27[i,8] = NA
    })
}
# IS beta table
table2.isw1[5,1]<-mean(b0.isw1.27,na.rm=TRUE)
table2.isw1[5,2]<-mean(b0.isw1.sd.27,na.rm=TRUE)
table2.isw1[5,3]<-sd(b0.isw1.27,na.rm=TRUE)
table2.isw1[5,4]<-mean(cover.isw1.27[,4],na.rm=TRUE)
table2.isw1[5,5]<-mean(b1.isw1.27,na.rm=TRUE)
table2.isw1[5,6]<-mean(b1.isw1.sd.27,na.rm=TRUE)
table2.isw1[5,7]<-sd(b1.isw1.27,na.rm=TRUE)
table2.isw1[5,8]<-mean(cover.isw1.27[,8],na.rm=TRUE)

table2.rcpp[5,1]<-mean(b0.rcpp.27,na.rm=TRUE)
table2.rcpp[5,2]<-mean(b0.rcpp.sd.27,na.rm=TRUE)
table2.rcpp[5,3]<-sd(b0.rcpp.27,na.rm=TRUE)
table2.rcpp[5,4]<-mean(cover.rcpp.27[,4],na.rm=TRUE)
table2.rcpp[5,5]<-mean(b1.rcpp.27,na.rm=TRUE)
table2.rcpp[5,6]<-mean(b1.rcpp.sd.27,na.rm=TRUE)
table2.rcpp[5,7]<-sd(b1.rcpp.27,na.rm=TRUE)
table2.rcpp[5,8]<-mean(cover.rcpp.27[,8],na.rm=TRUE)


#### censoring point at t_0=3 ####
c.0=5000000
c.1=61.17
c.3=22.53
c.5=13.43
c.7=8.5
#### t_0=3 & c=0% ####
b0.isw1.30 = c()
b0.isw1.sd.30 = c()
b1.isw1.30 = c()
b1.isw1.sd.30 = c()
cover.isw1.30=matrix(NA,2000,8)
colnames(cover.isw1.30)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.rcpp.30 = c()
b0.rcpp.sd.30 = c()
b1.rcpp.30 = c()
b1.rcpp.sd.30 = c()
cover.rcpp.30=matrix(NA,2000,8)
colnames(cover.rcpp.30)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.0)
  # method 1 : Our method (ISMB + rcpp's weight + weight-in + dfsane)
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw1.30[i] = ismb.fit[1,1]
    b0.isw1.sd.30[i] = ismb.fit[1,2]
    b1.isw1.30[i] = ismb.fit[2,1]
    b1.isw1.sd.30[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.30[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.30[i,2] = ismb.fit[1,1]
    cover.isw1.30[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.30[i,4] = ind(1.040613, cover.isw1.30[i,1], cover.isw1.30[i,3])
    cover.isw1.30[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.30[i,6] = ismb.fit[2,1]
    cover.isw1.30[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.30[i,8] = ind(1.0174711, cover.isw1.30[i,5], cover.isw1.30[i,7])}
    , error=function(e){
      b0.isw1.30[i] = NA
      b0.isw1.sd.30[i] = NA
      b1.isw1.30[i] = NA
      b1.isw1.sd.30[i] = NA
      # Coverage
      cover.isw1.30[i,1] = NA
      cover.isw1.30[i,2] = NA
      cover.isw1.30[i,3] = NA
      cover.isw1.30[i,4] = NA
      cover.isw1.30[i,5] = NA
      cover.isw1.30[i,6] = NA
      cover.isw1.30[i,7] = NA
      cover.isw1.30[i,8] = NA
    })

  # method 2 : rcpp's method (rcpp's weight + weight-in + dfsane)
  tryCatch({
    rcpp.fit = rcpp.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.rcpp.30[i] = rcpp.fit[1,1]
    b0.rcpp.sd.30[i] = rcpp.fit[1,2]
    b1.rcpp.30[i] = rcpp.fit[2,1]
    b1.rcpp.sd.30[i] = rcpp.fit[2,2]
    # Coverage
    cover.rcpp.30[i,1] = rcpp.fit[1,1]-1.96*rcpp.fit[1,2]
    cover.rcpp.30[i,2] = rcpp.fit[1,1]
    cover.rcpp.30[i,3] = rcpp.fit[1,1]+1.96*rcpp.fit[1,2]
    cover.rcpp.30[i,4] = ind(1.040613, cover.rcpp.30[i,1], cover.rcpp.30[i,3])
    cover.rcpp.30[i,5] = rcpp.fit[2,1]-1.96*rcpp.fit[2,2]
    cover.rcpp.30[i,6] = rcpp.fit[2,1]
    cover.rcpp.30[i,7] = rcpp.fit[2,1]+1.96*rcpp.fit[2,2]
    cover.rcpp.30[i,8] = ind(1.0174711, cover.rcpp.30[i,5], cover.rcpp.30[i,7])}
    , error=function(e){
      b0.rcpp.30[i] = NA
      b0.rcpp.sd.30[i] = NA
      b1.rcpp.30[i] = NA
      b1.rcpp.sd.30[i] = NA
      # Coverage
      cover.rcpp.30[i,1] = NA
      cover.rcpp.30[i,2] = NA
      cover.rcpp.30[i,3] = NA
      cover.rcpp.30[i,4] = NA
      cover.rcpp.30[i,5] = NA
      cover.rcpp.30[i,6] = NA
      cover.rcpp.30[i,7] = NA
      cover.rcpp.30[i,8] = NA
    })
}
# IS beta table
table3.isw1[1,1]<-mean(b0.isw1.30,na.rm=TRUE)
table3.isw1[1,2]<-mean(b0.isw1.sd.30,na.rm=TRUE)
table3.isw1[1,3]<-sd(b0.isw1.30,na.rm=TRUE)
table3.isw1[1,4]<-mean(cover.isw1.30[,4],na.rm=TRUE)
table3.isw1[1,5]<-mean(b1.isw1.30,na.rm=TRUE)
table3.isw1[1,6]<-mean(b1.isw1.sd.30,na.rm=TRUE)
table3.isw1[1,7]<-sd(b1.isw1.30,na.rm=TRUE)
table3.isw1[1,8]<-mean(cover.isw1.30[,8],na.rm=TRUE)

table3.rcpp[1,1]<-mean(b0.rcpp.30,na.rm=TRUE)
table3.rcpp[1,2]<-mean(b0.rcpp.sd.30,na.rm=TRUE)
table3.rcpp[1,3]<-sd(b0.rcpp.30,na.rm=TRUE)
table3.rcpp[1,4]<-mean(cover.rcpp.30[,4],na.rm=TRUE)
table3.rcpp[1,5]<-mean(b1.rcpp.30,na.rm=TRUE)
table3.rcpp[1,6]<-mean(b1.rcpp.sd.30,na.rm=TRUE)
table3.rcpp[1,7]<-sd(b1.rcpp.30,na.rm=TRUE)
table3.rcpp[1,8]<-mean(cover.rcpp.30[,8],na.rm=TRUE)

#### t_0=3 & c=10% ####
b0.isw1.31 = c()
b0.isw1.sd.31 = c()
b1.isw1.31 = c()
b1.isw1.sd.31 = c()
cover.isw1.31=matrix(NA,2000,8)
colnames(cover.isw1.31)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.rcpp.31 = c()
b0.rcpp.sd.31 = c()
b1.rcpp.31 = c()
b1.rcpp.sd.31 = c()
cover.rcpp.31=matrix(NA,2000,8)
colnames(cover.rcpp.31)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.1)

  # method 1 : Our method (ISMB + rcpp's weight + weight-in + dfsane)
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw1.31[i] = ismb.fit[1,1]
    b0.isw1.sd.31[i] = ismb.fit[1,2]
    b1.isw1.31[i] = ismb.fit[2,1]
    b1.isw1.sd.31[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.31[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.31[i,2] = ismb.fit[1,1]
    cover.isw1.31[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.31[i,4] = ind(1.040613, cover.isw1.31[i,1], cover.isw1.31[i,3])
    cover.isw1.31[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.31[i,6] = ismb.fit[2,1]
    cover.isw1.31[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.31[i,8] = ind(1.0174711, cover.isw1.31[i,5], cover.isw1.31[i,7])}
    , error=function(e){
      b0.isw1.31[i] = NA
      b0.isw1.sd.31[i] = NA
      b1.isw1.31[i] = NA
      b1.isw1.sd.31[i] = NA
      # Coverage
      cover.isw1.31[i,1] = NA
      cover.isw1.31[i,2] = NA
      cover.isw1.31[i,3] = NA
      cover.isw1.31[i,4] = NA
      cover.isw1.31[i,5] = NA
      cover.isw1.31[i,6] = NA
      cover.isw1.31[i,7] = NA
      cover.isw1.31[i,8] = NA
    })

  # method 2 : rcpp's method (rcpp's weight + weight-in + dfsane)
  tryCatch({
    rcpp.fit = rcpp.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.rcpp.31[i] = rcpp.fit[1,1]
    b0.rcpp.sd.31[i] = rcpp.fit[1,2]
    b1.rcpp.31[i] = rcpp.fit[2,1]
    b1.rcpp.sd.31[i] = rcpp.fit[2,2]
    # Coverage
    cover.rcpp.31[i,1] = rcpp.fit[1,1]-1.96*rcpp.fit[1,2]
    cover.rcpp.31[i,2] = rcpp.fit[1,1]
    cover.rcpp.31[i,3] = rcpp.fit[1,1]+1.96*rcpp.fit[1,2]
    cover.rcpp.31[i,4] = ind(1.040613, cover.rcpp.31[i,1], cover.rcpp.31[i,3])
    cover.rcpp.31[i,5] = rcpp.fit[2,1]-1.96*rcpp.fit[2,2]
    cover.rcpp.31[i,6] = rcpp.fit[2,1]
    cover.rcpp.31[i,7] = rcpp.fit[2,1]+1.96*rcpp.fit[2,2]
    cover.rcpp.31[i,8] = ind(1.0174711, cover.rcpp.31[i,5], cover.rcpp.31[i,7])}
    , error=function(e){
      b0.rcpp.31[i] = NA
      b0.rcpp.sd.31[i] = NA
      b1.rcpp.31[i] = NA
      b1.rcpp.sd.31[i] = NA
      # Coverage
      cover.rcpp.31[i,1] = NA
      cover.rcpp.31[i,2] = NA
      cover.rcpp.31[i,3] = NA
      cover.rcpp.31[i,4] = NA
      cover.rcpp.31[i,5] = NA
      cover.rcpp.31[i,6] = NA
      cover.rcpp.31[i,7] = NA
      cover.rcpp.31[i,8] = NA
    })
}
# IS beta table
table3.isw1[2,1]<-mean(b0.isw1.31,na.rm=TRUE)
table3.isw1[2,2]<-mean(b0.isw1.sd.31,na.rm=TRUE)
table3.isw1[2,3]<-sd(b0.isw1.31,na.rm=TRUE)
table3.isw1[2,4]<-mean(cover.isw1.31[,4],na.rm=TRUE)
table3.isw1[2,5]<-mean(b1.isw1.31,na.rm=TRUE)
table3.isw1[2,6]<-mean(b1.isw1.sd.31,na.rm=TRUE)
table3.isw1[2,7]<-sd(b1.isw1.31,na.rm=TRUE)
table3.isw1[2,8]<-mean(cover.isw1.31[,8],na.rm=TRUE)

table3.rcpp[2,1]<-mean(b0.rcpp.31,na.rm=TRUE)
table3.rcpp[2,2]<-mean(b0.rcpp.sd.31,na.rm=TRUE)
table3.rcpp[2,3]<-sd(b0.rcpp.31,na.rm=TRUE)
table3.rcpp[2,4]<-mean(cover.rcpp.31[,4],na.rm=TRUE)
table3.rcpp[2,5]<-mean(b1.rcpp.31,na.rm=TRUE)
table3.rcpp[2,6]<-mean(b1.rcpp.sd.31,na.rm=TRUE)
table3.rcpp[2,7]<-sd(b1.rcpp.31,na.rm=TRUE)
table3.rcpp[2,8]<-mean(cover.rcpp.31[,8],na.rm=TRUE)

#### t_0=3 & c=30% ####
b0.isw1.33 = c()
b0.isw1.sd.33 = c()
b1.isw1.33 = c()
b1.isw1.sd.33 = c()
cover.isw1.33=matrix(NA,2000,8)
colnames(cover.isw1.33)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.rcpp.33 = c()
b0.rcpp.sd.33 = c()
b1.rcpp.33 = c()
b1.rcpp.sd.33 = c()
cover.rcpp.33=matrix(NA,2000,8)
colnames(cover.rcpp.33)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.3)
  # method 1 : Our method (ISMB + rcpp's weight + weight-in + dfsane)
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw1.33[i] = ismb.fit[1,1]
    b0.isw1.sd.33[i] = ismb.fit[1,2]
    b1.isw1.33[i] = ismb.fit[2,1]
    b1.isw1.sd.33[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.33[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.33[i,2] = ismb.fit[1,1]
    cover.isw1.33[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.33[i,4] = ind(1.040613, cover.isw1.33[i,1], cover.isw1.33[i,3])
    cover.isw1.33[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.33[i,6] = ismb.fit[2,1]
    cover.isw1.33[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.33[i,8] = ind(1.0174711, cover.isw1.33[i,5], cover.isw1.33[i,7])}
    , error=function(e){
      b0.isw1.33[i] = NA
      b0.isw1.sd.33[i] = NA
      b1.isw1.33[i] = NA
      b1.isw1.sd.33[i] = NA
      # Coverage
      cover.isw1.33[i,1] = NA
      cover.isw1.33[i,2] = NA
      cover.isw1.33[i,3] = NA
      cover.isw1.33[i,4] = NA
      cover.isw1.33[i,5] = NA
      cover.isw1.33[i,6] = NA
      cover.isw1.33[i,7] = NA
      cover.isw1.33[i,8] = NA
    })

  # method 2 : rcpp's method (rcpp's weight + weight-in + dfsane)
  tryCatch({
    rcpp.fit = rcpp.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.rcpp.33[i] = rcpp.fit[1,1]
    b0.rcpp.sd.33[i] = rcpp.fit[1,2]
    b1.rcpp.33[i] = rcpp.fit[2,1]
    b1.rcpp.sd.33[i] = rcpp.fit[2,2]
    # Coverage
    cover.rcpp.33[i,1] = rcpp.fit[1,1]-1.96*rcpp.fit[1,2]
    cover.rcpp.33[i,2] = rcpp.fit[1,1]
    cover.rcpp.33[i,3] = rcpp.fit[1,1]+1.96*rcpp.fit[1,2]
    cover.rcpp.33[i,4] = ind(1.040613, cover.rcpp.33[i,1], cover.rcpp.33[i,3])
    cover.rcpp.33[i,5] = rcpp.fit[2,1]-1.96*rcpp.fit[2,2]
    cover.rcpp.33[i,6] = rcpp.fit[2,1]
    cover.rcpp.33[i,7] = rcpp.fit[2,1]+1.96*rcpp.fit[2,2]
    cover.rcpp.33[i,8] = ind(1.0174711, cover.rcpp.33[i,5], cover.rcpp.33[i,7])}
    , error=function(e){
      b0.rcpp.33[i] = NA
      b0.rcpp.sd.33[i] = NA
      b1.rcpp.33[i] = NA
      b1.rcpp.sd.33[i] = NA
      # Coverage
      cover.rcpp.33[i,1] = NA
      cover.rcpp.33[i,2] = NA
      cover.rcpp.33[i,3] = NA
      cover.rcpp.33[i,4] = NA
      cover.rcpp.33[i,5] = NA
      cover.rcpp.33[i,6] = NA
      cover.rcpp.33[i,7] = NA
      cover.rcpp.33[i,8] = NA
    })
}

# IS beta table
table3.isw1[3,1]<-mean(b0.isw1.33,na.rm=TRUE)
table3.isw1[3,2]<-mean(b0.isw1.sd.33,na.rm=TRUE)
table3.isw1[3,3]<-sd(b0.isw1.33,na.rm=TRUE)
table3.isw1[3,4]<-mean(cover.isw1.33[,4],na.rm=TRUE)
table3.isw1[3,5]<-mean(b1.isw1.33,na.rm=TRUE)
table3.isw1[3,6]<-mean(b1.isw1.sd.33,na.rm=TRUE)
table3.isw1[3,7]<-sd(b1.isw1.33,na.rm=TRUE)
table3.isw1[3,8]<-mean(cover.isw1.33[,8],na.rm=TRUE)

table3.rcpp[3,1]<-mean(b0.rcpp.33,na.rm=TRUE)
table3.rcpp[3,2]<-mean(b0.rcpp.sd.33,na.rm=TRUE)
table3.rcpp[3,3]<-sd(b0.rcpp.33,na.rm=TRUE)
table3.rcpp[3,4]<-mean(cover.rcpp.33[,4],na.rm=TRUE)
table3.rcpp[3,5]<-mean(b1.rcpp.33,na.rm=TRUE)
table3.rcpp[3,6]<-mean(b1.rcpp.sd.33,na.rm=TRUE)
table3.rcpp[3,7]<-sd(b1.rcpp.33,na.rm=TRUE)
table3.rcpp[3,8]<-mean(cover.rcpp.33[,8],na.rm=TRUE)

#### t_0=3 & c=50% ####
b0.isw1.35 = c()
b0.isw1.sd.35 = c()
b1.isw1.35 = c()
b1.isw1.sd.35 = c()
cover.isw1.35=matrix(NA,2000,8)
colnames(cover.isw1.35)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.rcpp.35 = c()
b0.rcpp.sd.35 = c()
b1.rcpp.35 = c()
b1.rcpp.sd.35 = c()
cover.rcpp.35=matrix(NA,2000,8)
colnames(cover.rcpp.35)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.5)
  # method 1 : Our method (ISMB + rcpp's weight + weight-in + dfsane)
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw1.35[i] = ismb.fit[1,1]
    b0.isw1.sd.35[i] = ismb.fit[1,2]
    b1.isw1.35[i] = ismb.fit[2,1]
    b1.isw1.sd.35[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.35[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.35[i,2] = ismb.fit[1,1]
    cover.isw1.35[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.35[i,4] = ind(1.040613, cover.isw1.35[i,1], cover.isw1.35[i,3])
    cover.isw1.35[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.35[i,6] = ismb.fit[2,1]
    cover.isw1.35[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.35[i,8] = ind(1.0174711, cover.isw1.35[i,5], cover.isw1.35[i,7])}
    , error=function(e){
      b0.isw1.35[i] = NA
      b0.isw1.sd.35[i] = NA
      b1.isw1.35[i] = NA
      b1.isw1.sd.35[i] = NA
      # Coverage
      cover.isw1.35[i,1] = NA
      cover.isw1.35[i,2] = NA
      cover.isw1.35[i,3] = NA
      cover.isw1.35[i,4] = NA
      cover.isw1.35[i,5] = NA
      cover.isw1.35[i,6] = NA
      cover.isw1.35[i,7] = NA
      cover.isw1.35[i,8] = NA
    })

  # method 2 : rcpp's method (rcpp's weight + weight-in + dfsane)
  tryCatch({
    rcpp.fit = rcpp.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.rcpp.35[i] = rcpp.fit[1,1]
    b0.rcpp.sd.35[i] = rcpp.fit[1,2]
    b1.rcpp.35[i] = rcpp.fit[2,1]
    b1.rcpp.sd.35[i] = rcpp.fit[2,2]
    # Coverage
    cover.rcpp.35[i,1] = rcpp.fit[1,1]-1.96*rcpp.fit[1,2]
    cover.rcpp.35[i,2] = rcpp.fit[1,1]
    cover.rcpp.35[i,3] = rcpp.fit[1,1]+1.96*rcpp.fit[1,2]
    cover.rcpp.35[i,4] = ind(1.040613, cover.rcpp.35[i,1], cover.rcpp.35[i,3])
    cover.rcpp.35[i,5] = rcpp.fit[2,1]-1.96*rcpp.fit[2,2]
    cover.rcpp.35[i,6] = rcpp.fit[2,1]
    cover.rcpp.35[i,7] = rcpp.fit[2,1]+1.96*rcpp.fit[2,2]
    cover.rcpp.35[i,8] = ind(1.0174711, cover.rcpp.35[i,5], cover.rcpp.35[i,7])}
    , error=function(e){
      b0.rcpp.35[i] = NA
      b0.rcpp.sd.35[i] = NA
      b1.rcpp.35[i] = NA
      b1.rcpp.sd.35[i] = NA
      # Coverage
      cover.rcpp.35[i,1] = NA
      cover.rcpp.35[i,2] = NA
      cover.rcpp.35[i,3] = NA
      cover.rcpp.35[i,4] = NA
      cover.rcpp.35[i,5] = NA
      cover.rcpp.35[i,6] = NA
      cover.rcpp.35[i,7] = NA
      cover.rcpp.35[i,8] = NA
    })
}
# IS beta table
table3.isw1[4,1]<-mean(b0.isw1.35,na.rm=TRUE)
table3.isw1[4,2]<-mean(b0.isw1.sd.35,na.rm=TRUE)
table3.isw1[4,3]<-sd(b0.isw1.35,na.rm=TRUE)
table3.isw1[4,4]<-mean(cover.isw1.35[,4],na.rm=TRUE)
table3.isw1[4,5]<-mean(b1.isw1.35,na.rm=TRUE)
table3.isw1[4,6]<-mean(b1.isw1.sd.35,na.rm=TRUE)
table3.isw1[4,7]<-sd(b1.isw1.35,na.rm=TRUE)
table3.isw1[4,8]<-mean(cover.isw1.35[,8],na.rm=TRUE)

table3.rcpp[4,1]<-mean(b0.rcpp.35,na.rm=TRUE)
table3.rcpp[4,2]<-mean(b0.rcpp.sd.35,na.rm=TRUE)
table3.rcpp[4,3]<-sd(b0.rcpp.35,na.rm=TRUE)
table3.rcpp[4,4]<-mean(cover.rcpp.35[,4],na.rm=TRUE)
table3.rcpp[4,5]<-mean(b1.rcpp.35,na.rm=TRUE)
table3.rcpp[4,6]<-mean(b1.rcpp.sd.35,na.rm=TRUE)
table3.rcpp[4,7]<-sd(b1.rcpp.35,na.rm=TRUE)
table3.rcpp[4,8]<-mean(cover.rcpp.35[,8],na.rm=TRUE)

#### t_0=3 & c=70% ####
b0.isw1.37 = c()
b0.isw1.sd.37 = c()
b1.isw1.37 = c()
b1.isw1.sd.37 = c()
cover.isw1.37=matrix(NA,2000,8)
colnames(cover.isw1.37)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

b0.rcpp.37 = c()
b0.rcpp.sd.37 = c()
b1.rcpp.37 = c()
b1.rcpp.sd.37 = c()
cover.rcpp.37=matrix(NA,2000,8)
colnames(cover.rcpp.37)=c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.7)
  # method 1 : Our method (ISMB + rcpp's weight + weight-in + dfsane)
  tryCatch({
    ismb.fit = ismbw1.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.isw1.37[i] = ismb.fit[1,1]
    b0.isw1.sd.37[i] = ismb.fit[1,2]
    b1.isw1.37[i] = ismb.fit[2,1]
    b1.isw1.sd.37[i] = ismb.fit[2,2]
    # Coverage
    cover.isw1.37[i,1] = ismb.fit[1,1]-1.96*ismb.fit[1,2]
    cover.isw1.37[i,2] = ismb.fit[1,1]
    cover.isw1.37[i,3] = ismb.fit[1,1]+1.96*ismb.fit[1,2]
    cover.isw1.37[i,4] = ind(1.040613, cover.isw1.37[i,1], cover.isw1.37[i,3])
    cover.isw1.37[i,5] = ismb.fit[2,1]-1.96*ismb.fit[2,2]
    cover.isw1.37[i,6] = ismb.fit[2,1]
    cover.isw1.37[i,7] = ismb.fit[2,1]+1.96*ismb.fit[2,2]
    cover.isw1.37[i,8] = ind(1.0174711, cover.isw1.37[i,5], cover.isw1.37[i,7])}
    , error=function(e){
      b0.isw1.37[i] = NA
      b0.isw1.sd.37[i] = NA
      b1.isw1.37[i] = NA
      b1.isw1.sd.37[i] = NA
      # Coverage
      cover.isw1.37[i,1] = NA
      cover.isw1.37[i,2] = NA
      cover.isw1.37[i,3] = NA
      cover.isw1.37[i,4] = NA
      cover.isw1.37[i,5] = NA
      cover.isw1.37[i,6] = NA
      cover.isw1.37[i,7] = NA
      cover.isw1.37[i,8] = NA
    })

  # method 2 : rcpp's method (rcpp's weight + weight-in + dfsane)
  tryCatch({
    rcpp.fit = rcpp.est(a[,3], 1, a[,4], a[,5], 3, 0.5, 200)
    b0.rcpp.37[i] = rcpp.fit[1,1]
    b0.rcpp.sd.37[i] = rcpp.fit[1,2]
    b1.rcpp.37[i] = rcpp.fit[2,1]
    b1.rcpp.sd.37[i] = rcpp.fit[2,2]
    # Coverage
    cover.rcpp.37[i,1] = rcpp.fit[1,1]-1.96*rcpp.fit[1,2]
    cover.rcpp.37[i,2] = rcpp.fit[1,1]
    cover.rcpp.37[i,3] = rcpp.fit[1,1]+1.96*rcpp.fit[1,2]
    cover.rcpp.37[i,4] = ind(1.040613, cover.rcpp.37[i,1], cover.rcpp.37[i,3])
    cover.rcpp.37[i,5] = rcpp.fit[2,1]-1.96*rcpp.fit[2,2]
    cover.rcpp.37[i,6] = rcpp.fit[2,1]
    cover.rcpp.37[i,7] = rcpp.fit[2,1]+1.96*rcpp.fit[2,2]
    cover.rcpp.37[i,8] = ind(1.0174711, cover.rcpp.37[i,5], cover.rcpp.37[i,7])}
    , error=function(e){
      b0.rcpp.37[i] = NA
      b0.rcpp.sd.37[i] = NA
      b1.rcpp.37[i] = NA
      b1.rcpp.sd.37[i] = NA
      # Coverage
      cover.rcpp.37[i,1] = NA
      cover.rcpp.37[i,2] = NA
      cover.rcpp.37[i,3] = NA
      cover.rcpp.37[i,4] = NA
      cover.rcpp.37[i,5] = NA
      cover.rcpp.37[i,6] = NA
      cover.rcpp.37[i,7] = NA
      cover.rcpp.37[i,8] = NA
    })
}
# IS beta table
table3.isw1[5,1]<-mean(b0.isw1.37,na.rm=TRUE)
table3.isw1[5,2]<-mean(b0.isw1.sd.37,na.rm=TRUE)
table3.isw1[5,3]<-sd(b0.isw1.37,na.rm=TRUE)
table3.isw1[5,4]<-mean(cover.isw1.37[,4],na.rm=TRUE)
table3.isw1[5,5]<-mean(b1.isw1.37,na.rm=TRUE)
table3.isw1[5,6]<-mean(b1.isw1.sd.37,na.rm=TRUE)
table3.isw1[5,7]<-sd(b1.isw1.37,na.rm=TRUE)
table3.isw1[5,8]<-mean(cover.isw1.37[,8],na.rm=TRUE)

table3.rcpp[5,1]<-mean(b0.rcpp.37,na.rm=TRUE)
table3.rcpp[5,2]<-mean(b0.rcpp.sd.37,na.rm=TRUE)
table3.rcpp[5,3]<-sd(b0.rcpp.37,na.rm=TRUE)
table3.rcpp[5,4]<-mean(cover.rcpp.37[,4],na.rm=TRUE)
table3.rcpp[5,5]<-mean(b1.rcpp.37,na.rm=TRUE)
table3.rcpp[5,6]<-mean(b1.rcpp.sd.37,na.rm=TRUE)
table3.rcpp[5,7]<-sd(b1.rcpp.37,na.rm=TRUE)
table3.rcpp[5,8]<-mean(cover.rcpp.37[,8],na.rm=TRUE)
