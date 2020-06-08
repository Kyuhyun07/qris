#' Quantile regression estimator (Kim et al & optimization & weight out)
#'
#' rq_optim_est function will calculate qunatile regression parameters "beta" using quantile regression idea. (Kim et al. 2012, function (8))
#' It find solution from optimizer (optim)
#'
#' @param Z A vector of subjects' event time
#' @param nc Number of covariate to estimate
#' @param covariate A matrix of covariate (dimension number of subjects(n) x number of covariate(nc))
#' @param D A vector of subjects' censored indicator (1 = not_censored data, 0 = censored data)
#' @param t_0 Followup time or base time adjustment (minimum t_0 is 0, and it means no base time adjustment)
#' @param Q Quantile that user want to estimate (scale : 0 ~ 1)
#' @param W A vector of weights
#' @return A nc by 2 matrix (First column is beta and second column is standard error of betas)
#' @examples
#'     covar = as.matrix(data[,1:12])
#'     ismb_est(data$survTime, 12, covar, data$event, 2, 0.5 ,100)
#' @export
rq_optim_est = function(Z, nc, covariate, D, t_0, Q, W){
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

  # Estimating function to minimize result
  rq_optim_objectF = function(beta){
    beta = as.matrix(beta)
    result = t(W*I*(logT-(X%*%beta)))%*%(Q-ifelse(logT-(X%*%beta)<=0,1,0))
  }

  # initial beta for using optimizer
  rq.optim.fit = optim(par = c(1,rep(0,nc)) , fn = rq_optim_objectF)
  # When optimizer find converged solution, we take that value
  if (rq.optim.fit$convergence == 0){
    solbeta = rq.optim.fit$par
    print(solbeta)
    # When optimizer makes any error, we print out NA
  } else {
    solbeta = c(NA, NA)
    print(solbeta)
  }
}
