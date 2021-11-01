#### 1. simulation data ####
# a single binary covariate with success probability 0.5.
# We generate failure time, T, from a Weibull distribution.
# Potential censoring time, C, from uniform (0,c) independently from T.
# Given condition r=rho_0, k=2, exp(beta_0)=5, exp(beta_1)=2

# 1.1 data generation function
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

# 1.2 Given value
exp.beta.initial.0=5
exp.beta.initial.1=10
k=2
r.initial.0=(log(10/5))^(1/k)/exp.beta.initial.0
r.initial.1=(log(10/5))^(1/k)/exp.beta.initial.1
# c value of unif(0,c)
c=24.35
a<-data.gen(200,c)
Z=a[,3]
covariate=as.matrix(a[,4])
D=a[,5]
t_0=1
Q=0.5
ne=200

fit1 <- qrismb(Z, covariate, D, t_0, Q, ne, "rq", "smooth")
fit2 <- qrismb(Z, covariate, D, t_0, Q, ne, "random", "nonsmooth")

data(cancer, package="survival")
lung_rev <- lung[,c(2,3,4,5,7,8,9)]
# delta : censored = 0, death = 1
lung_rev$status <- as.numeric(lung$status) - 1
# covariate 1 : age (39~82)
# covariate 2 : sex(male = 0, female = 1)
lung_rev$sex <- as.numeric(lung$sex) - 1
# covariate 3 : ph.karno: Karnofsky performance score (bad=0-good=100) rated by physician
# covariate 4 : pat.karno: Karnofsky performance score as rated by patient
# covariate 5 : meal.cal: Calories consumed at meals
# covariate 6 : wt.loss: Weight loss in last six months (pounds)
Z <- lung_rev$time
covariate <- as.matrix(lung_rev[,c(3:6)])
D <- lung_rev$status
t_0 <- 5
Q <- 0.5
ne <- 200

fit1 <- qrismb(Z, covariate, D, t_0, Q, ne, "rq", "smooth")
fit2 <- qrismb(Z, covariate, D, t_0, Q, ne, "random", "nonsmooth")

coef(fit1)
coef(fit2)
