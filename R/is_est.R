#' Quantile regression estimator using induced smoothing approach
#'
#' is_est function will calculate qunatile regression parameters "beta" using quantile regression idea and induced smoothing approach.
#' (Kim et al. 2012, equation (9)) + induced smoothing = my manuscript equation(3)
#' It directly find solution from non-linear equation solver (nleqslv)
#'
#' @param Z A vector of subjects' event time
#' @param nc Number of covariate to estimate
#' @param covariate A matrix of covariate (dimension number of subjects(n) x number of covariate(nc))
#' @param D A vector of subjects' censored indicator (1 = not_censored data, 0 = censored data)
#' @param t_0 Followup time or base time adjustment (minimum t_0 is 0, and it means no base time adjustment)
#' @param Q Quantile that user want to estimate (scale : 0 ~ 1)
#' @param W A vector of weights
#' @param ne Number of eta generation when estimate standard error of estimator
#' @return A nc by 2 matrix (First column is beta and second column is standard error of betas)
#' @examples
#'     covar = as.matrix(data[,1:12])
#'     ismb_est(data$survTime, 12, covar, data$event, 2, 0.5 ,100)
#' @export
is_est = function(Z, nc, covariate, D, W, t_0, Q, ne){
  # n = number of subject
  n = length(Z)
  data = matrix(NA, n, nc+5)
  # event time = minimum of event time and censored time of subject
  data[,1] = Z
  # Take log for event time - adjusted followup time
  data[,2] = log(Z-t_0)
  # Indicator of event time > followup time
  data[,3] = as.numeric(Z>t_0)
  # covariates
  data[,4:(nc+3)] = covariate
  # Change log(negative number)=NA to -10 (any number is possible)
  data[is.na(data[,2]),2]=-10
  # censoring indicator (1=uncensored data, 0=censored data)
  data[,(nc+4)] = D
  # change last data to uncensored data (to improve performance of estimation)
  data[n,(nc+4)] = 1
  data[,(nc+5)] = W
  data = as.data.frame(data)

  # Naming data
  colnames(data)[1:3]=c("Z", "log(Z-t_0)","I[Z>t_0]")
  for (j in 4:(nc+3)){
    colnames(data)[j]="covariate"
  }
  colnames(data)[(nc+4):(nc+5)] = c("delta","Weight")

  # Covariate setting (1 covariate)
  # Making covariate matrix (intercept coefficient + covariate)
  X = as.matrix(cbind(c(rep(1,n)),data[,4:(nc+3)]))
  # response variable
  logT = data[,2]
  # Indicator of event time > followup time
  I = data[,3]
  # matrix which is used for induced smoothing
  H = diag(1/n, nc+1, nc+1)

  # Estimating equation for estimaing beta
  is_objectF = function(beta){
    beta = as.matrix(beta)
    result = t(X*I*W) %*% {(pnorm((X%*%beta-logT)/sqrt(diag(X %*% H %*% t(X)))))-Q}
  }

  # initial beta for using non-linear equation solver
  betastart = c(1,rep(1,nc))
  is.fit = nleqslv(betastart, is_objectF, control=list(ftol=1e-5))
  # When solver find converged solution, we take that value
  if (is.fit$termcd == 1){
    solbeta = is.fit$x
    print(solbeta)
    # When solver makes any error, we print out NA
  } else {
    solbeta = c(NA, NA)
    print(solbeta)
  }
}
