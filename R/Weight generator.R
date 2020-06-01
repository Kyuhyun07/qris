#### Weight generator ####
weight_generator = function(Z, t_0, nc, covariate, D){
  n = length(Z)
  data = matrix(NA, n, nc+6)
  data[,1] = Z
  data[,2] = log(Z-t_0)
  data[,3] = as.numeric(Z>t_0)
  data[,4:(nc+3)] = covariate
  data[is.na(data[,2]),2]=-10
  data[,(nc+4)] = D
  data[n,(nc+4)] = 1
  data = as.data.frame(data)

  # Weight
  # Kaplan-Meier estimator for censoring survfit weight
  # fit = survfit(Surv(data[,1],1-(data[,(nc+4)])) ~ 1)
  # for (i in 1:length(fit$surv)){
  #   data[data[,1]==fit$time[i],(nc+5)] = fit$surv[i]
  #   }
  # data[,(nc+6)] = data[,(nc+4)]/data[,(nc+5)]
  # m = nrow(data)
  # if (data[m,(nc+6)]==Inf){
  #   data[m,(nc+6)]=max(data[1:(m-1),(nc+6)])
  #   }
  # if (data[m,(nc+4)]==0){
  #   data[m,(nc+6)]=0
  # }

  # WKM surv weight
  fit = WKM(data[,1],  1-data[,(nc+4)], zc = rep(1,n), w = rep(1,n))
  for (i in 1:length(fit$surv)){
    data[data[,1]==fit$times[i],(nc+5)] = fit$surv[i]
  }
  data[,(nc+6)] = data[,(nc+4)]/data[,(nc+5)]
  m = nrow(data)
  if (data[m,(nc+6)]==Inf){
    data[m,(nc+6)]=max(data[1:(m-1),(nc+6)])
  }
  if (data[m,(nc+4)]==0){
    data[m,(nc+6)]=0
  }

  # weight using WKM jump
  # fit = WKM(data[,1], (data[,(nc+4)]), zc = rep(1,n), w = rep(1,n))
  # data[,(nc+6)] = fit$jump
  # m = nrow(data)
  colnames(data)[1:3]=c("Z", "log(Z-t_0)","I[Z>t_0]")
  for (j in 4:(nc+3)){
    colnames(data)[j]="covariate"
  }
  colnames(data)[(nc+4):(nc+6)] = c("delta","G_KM","Weight")
  return(data[,(nc+6)])
}
