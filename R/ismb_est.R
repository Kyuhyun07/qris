#' Quantile regression estimator using induced smoothing approach
#'
#' ismb_est function will calculate qunatile regression parameters "beta" and standard error of "beta" using induced smoothing approach.
#'
#' @param Z A vector of subjects' event time
#' @param nc Number of covariate to estimate
#' @param covariate A matrix of covariate (dimension number of subjects(n) x number of covariate(nc))
#' @param D A vector of subjects' censored indicator (1 = not_censored data, 0 = censored data)
#' @param t_0 Followup time or base time adjustment (minimum t_0 is 0, and it means no base time adjustment)
#' @param Q Quantile that user want to estimate (scale : 0 ~ 1)
#' @param ne Number of eta generation when estimate standard error of estimator
#' @return A nc by 2 matrix (First column is beta and second column is standard error of betas)
#' @examples
#'     covar = as.matrix(data[,1:12])
#'     ismb_est(data$survTime, 12, covar, data$event, 2, 0.5 ,100)
#' @export
ismb_est = function(Z, nc, covariate, D, t_0, Q, ne){
  n = length(Z)
  data = matrix(NA, n, nc+6)
  data[,1] = Z
  data[,2] = log(Z-t_0)
  data[,3] = as.numeric(Z>t_0)
  data[,4:(nc+3)] = covariate
  data[,(nc+4)] = D
  data[is.nan(data[,2]),2]=-10
  data = as.data.frame(data)

  # Weight
  # Kaplan-Meier estimator for censoring
  # survfit weight
  # fit = survfit(Surv(data[,1],(1-data[,(nc+4)])) ~ 1)
  # for (i in 1:length(fit$surv))
  # {
  #   data[data[,1]==fit$time[i],(nc+5)] = fit$surv[i]
  # }
  # WKM weight
  fit = WKM(data[,1],  1-data[,(nc+4)], zc = rep(1,n), w = rep(1,n))
  for (i in 1:length(fit$surv))
  {
    data[data[,1]==fit$times[i],(nc+5)] = fit$surv[i]
  }
  data[,(nc+6)] = data[,(nc+4)]/data[,(nc+5)]
  m = nrow(data)
  if (data[m,(nc+4)]==0){
    data[m,(nc+6)]=0
  } else {
    max(data[1:(m-1),(nc+6)])
  }
  colnames(data)[1:3]=c("Z", "log(Z-t_0)","I[T>=C]")
  for (j in 4:(nc+3)){
    colnames(data)[j]="covariate"
  }
  colnames(data)[(nc+4):(nc+6)] = c("delta","G_KM","Weight")

  # Objective equation
  objectF = function(beta){
    beta = as.matrix(beta)
    result = t(X*I*W) %*% {Q - (pnorm((X%*%beta-logT)/sqrt(diag(X %*% H %*% t(X)))))}
  }

  #### revised object equation ####
  rev.objectF = function(beta){
    beta = as.matrix(beta)
    result = t(eta*X*I*W) %*% {Q - (pnorm((X%*%beta-logT)/sqrt(diag(X %*% H %*% t(X)))))}
  }

  # Covariate setting (1 covariate)
  X = as.matrix(cbind(c(rep(1,n)),data[,4:(nc+3)]))
  W = data[,(nc+6)]
  logT = data[,2]
  I = data[,3]
  H = diag(1/n, nc+1, nc+1)

  # Change betastart when real data analysis c(1,rep(1,nc))
  betastart = c(-2,rep(0,nc))
  is.fit = nleqslv(betastart,objectF, control=list(ftol=1e-5))
  if (all(is.fit$fvec<1e-5&is.fit$fvec>-1e-5)){
    solbeta = is.fit$x
    # Variance estimation : ISMB
    result.ismb=c()
    for (j in 1:ne){
      eta = rexp(n,1)
      result = t(eta*X*I*W) %*% {Q - W*(pnorm((X%*%solbeta-logT)/sqrt(diag(X %*% H %*% t(X)))))}
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
