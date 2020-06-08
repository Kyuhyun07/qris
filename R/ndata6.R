#### ndata6(200608) : Condition ####
# data size = 200
# only beta0 effective
# Quantile 75%
# simulation 2000
# eta = 100
# WKM Weight outside, G(Z_i)

library(quantreg)
library(survival)
library(nleqslv)
library(xtable)
library(emplik)
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
    sim[q,1]={{-log(1-unif[q])}^(1/k)}/r.initial.0
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
k=2
r.initial.0=(log(10/5))^(1/k)/exp.beta.initial.0


#### Example parameter ####
# c.0=5000000
# c.1=123.69
# c.3=41.27
# c.5=23.55
# c.7=14.2
# a<-data.gen(200,c.3)
# Z=a[,3]
# nc=1
# covariate=a[,4]
# D=a[,5]
# t_0=0
# Q=0.75
# ne=100

#### Make table for Beta estimation and variance estimation and Coverage ####
table0<-matrix(NA,5,12)
rownames(table0)<-c(0,10,30,50,70)
colnames(table0)<-c("is_b0","is_b1","is.NA","is_optim_b0","is_optim_b1","is_optim.na","rq_b0","rq_b1","rq.NA","rq_optim_b0","rq_optim_b1","rq_optim.NA")

table1<-matrix(NA,5,12)
rownames(table1)<-c(0,10,30,50,70)
colnames(table1)<-c("is_b0","is_b1","is.NA","is_optim_b0","is_optim_b1","is_optim.na","rq_b0","rq_b1","rq.NA","rq_optim_b0","rq_optim_b1","rq_optim.NA")

table2<-matrix(NA,5,12)
rownames(table2)<-c(0,10,30,50,70)
colnames(table2)<-c("is_b0","is_b1","is.NA","is_optim_b0","is_optim_b1","is_optim.na","rq_b0","rq_b1","rq.NA","rq_optim_b0","rq_optim_b1","rq_optim.NA")

table3<-matrix(NA,5,12)
rownames(table3)<-c(0,10,30,50,70)
colnames(table3)<-c("is_b0","is_b1","is.NA","is_optim_b0","is_optim_b1","is_optim.na","rq_b0","rq_b1","rq.NA","rq_optim_b0","rq_optim_b1","rq_optim.NA")

#### censoring point at t_0=0 ####
c.0=5000000
c.1=36.98
c.3=12.53
c.5=7.47
c.7=4.81
#### t_0=0 & c=0% ####
b0.is.00 = c()
b1.is.00 = c()
b0.is.optim.00 = c()
b1.is.optim.00 = c()
b0.rq.00 = c()
b1.rq.00 = c()
b0.rq.optim.00 = c()
b1.rq.optim.00 = c()

set.seed(1)
for (i in 1:2000){
  a = data.gen(200,c.0)
  tic()
  is.fit = is_est(a[,3], 1, a[,4], a[,5], 0, 0.75)
  toc()
  b0.is.00[i] = is.fit[1]
  b1.is.00[i] = is.fit[2]
  is.optim.fit = is_optim_est(a[,3], 1, a[,4], a[,5], 0, 0.75)
  b0.is.optim.00[i] = is.optim.fit[1]
  b1.is.optim.00[i] = is.optim.fit[2]
  rq.fit = rq_est(a[,3], 1, a[,4], a[,5], 0, 0.75)
  b0.rq.00[i] = rq.fit[1]
  b1.rq.00[i] = rq.fit[2]
  rq.optim.fit = rq_optim_est(a[,3], 1, a[,4], a[,5], 0, 0.75)
  b0.rq.optim.00[i] = rq.optim.fit[1]
  b1.rq.optim.00[i] = rq.optim.fit[2]
}
table0[1,1] = mean(b0.is.00,na.rm=TRUE)
table0[1,2] = mean(b1.is.00,na.rm=TRUE)
table0[1,3] = sum(is.na(b0.is.00))
table0[1,4] = mean(b0.is.optim.00,na.rm=TRUE)
table0[1,5] = mean(b1.is.optim.00,na.rm=TRUE)
table0[1,6] = sum(is.na(b0.is.optim.00))
table0[1,7] = mean(b0.rq.00,na.rm=TRUE)
table0[1,8] = mean(b1.rq.00,na.rm=TRUE)
table0[1,9] = sum(is.na(b0.rq.00))
table0[1,10] = mean(b0.rq.optim.00,na.rm=TRUE)
table0[1,11] = mean(b1.rq.optim.00,na.rm=TRUE)
table0[1,12] = sum(is.na(b0.rq.optim.00))

#### t_0=0 & c=10% ####
b0.is.01 = c()
b1.is.01 = c()
b0.is.optim.01 = c()
b1.is.optim.01 = c()
b0.rq.01 = c()
b1.rq.01 = c()
b0.rq.optim.01 = c()
b1.rq.optim.01 = c()

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.1)
  is.fit = is_est(a[,3], 1, a[,4], a[,5], 0, 0.75)
  b0.is.01[i] = is.fit[1]
  b1.is.01[i] = is.fit[2]
  is.optim.fit = is_optim_est(a[,3], 1, a[,4], a[,5], 0, 0.75)
  b0.is.optim.01[i] = is.optim.fit[1]
  b1.is.optim.01[i] = is.optim.fit[2]
  rq.fit = rq_est(a[,3], 1, a[,4], a[,5], 0, 0.75)
  b0.rq.01[i] = rq.fit[1]
  b1.rq.01[i] = rq.fit[2]
  rq.optim.fit = rq_optim_est(a[,3], 1, a[,4], a[,5], 0, 0.75)
  b0.rq.optim.01[i] = rq.optim.fit[1]
  b1.rq.optim.01[i] = rq.optim.fit[2]
}
table0[2,1] = mean(b0.is.01,na.rm=TRUE)
table0[2,2] = mean(b1.is.01,na.rm=TRUE)
table0[2,3] = sum(is.na(b0.is.01))
table0[2,4] = mean(b0.is.optim.01,na.rm=TRUE)
table0[2,5] = mean(b1.is.optim.01,na.rm=TRUE)
table0[2,6] = sum(is.na(b0.is.optim.01))
table0[2,7] = mean(b0.rq.01,na.rm=TRUE)
table0[2,8] = mean(b1.rq.01,na.rm=TRUE)
table0[2,9] = sum(is.na(b0.rq.01))
table0[2,10] = mean(b0.rq.optim.01,na.rm=TRUE)
table0[2,11] = mean(b1.rq.optim.01,na.rm=TRUE)
table0[2,12] = sum(is.na(b0.rq.optim.01))

#### t_0=0 & c=30% ####
b0.is.03 = c()
b1.is.03 = c()
b0.is.optim.03 = c()
b1.is.optim.03 = c()
b0.rq.03 = c()
b1.rq.03 = c()
b0.rq.optim.03 = c()
b1.rq.optim.03 = c()

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.3)
  is.fit = is_est(a[,3], 1, a[,4], a[,5], 0, 0.75)
  b0.is.03[i] = is.fit[1]
  b1.is.03[i] = is.fit[2]
  is.optim.fit = is_optim_est(a[,3], 1, a[,4], a[,5], 0, 0.75)
  b0.is.optim.03[i] = is.optim.fit[1]
  b1.is.optim.03[i] = is.optim.fit[2]
  rq.fit = rq_est(a[,3], 1, a[,4], a[,5], 0, 0.75)
  b0.rq.03[i] = rq.fit[1]
  b1.rq.03[i] = rq.fit[2]
  rq.optim.fit = rq_optim_est(a[,3], 1, a[,4], a[,5], 0, 0.75)
  b0.rq.optim.03[i] = rq.optim.fit[1]
  b1.rq.optim.03[i] = rq.optim.fit[2]
}
table0[3,1] = mean(b0.is.03,na.rm=TRUE)
table0[3,2] = mean(b1.is.03,na.rm=TRUE)
table0[3,3] = sum(is.na(b0.is.03))
table0[3,4] = mean(b0.is.optim.03,na.rm=TRUE)
table0[3,5] = mean(b1.is.optim.03,na.rm=TRUE)
table0[3,6] = sum(is.na(b0.is.optim.03))
table0[3,7] = mean(b0.rq.03,na.rm=TRUE)
table0[3,8] = mean(b1.rq.03,na.rm=TRUE)
table0[3,9] = sum(is.na(b0.rq.03))
table0[3,10] = mean(b0.rq.optim.03,na.rm=TRUE)
table0[3,11] = mean(b1.rq.optim.03,na.rm=TRUE)
table0[3,12] = sum(is.na(b0.rq.optim.03))

#### t_0=0 & c=50% ####
b0.is.05 = c()
b1.is.05 = c()
b0.is.optim.05 = c()
b1.is.optim.05 = c()
b0.rq.05 = c()
b1.rq.05 = c()
b0.rq.optim.05 = c()
b1.rq.optim.05 = c()

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.5)
  is.fit = is_est(a[,3], 1, a[,4], a[,5], 0, 0.75)
  b0.is.05[i] = is.fit[1]
  b1.is.05[i] = is.fit[2]
  is.optim.fit = is_optim_est(a[,3], 1, a[,4], a[,5], 0, 0.75)
  b0.is.optim.05[i] = is.optim.fit[1]
  b1.is.optim.05[i] = is.optim.fit[2]
  rq.fit = rq_est(a[,3], 1, a[,4], a[,5], 0, 0.75)
  b0.rq.05[i] = rq.fit[1]
  b1.rq.05[i] = rq.fit[2]
  rq.optim.fit = rq_optim_est(a[,3], 1, a[,4], a[,5], 0, 0.75)
  b0.rq.optim.05[i] = rq.optim.fit[1]
  b1.rq.optim.05[i] = rq.optim.fit[2]
}
table0[4,1] = mean(b0.is.05,na.rm=TRUE)
table0[4,2] = mean(b1.is.05,na.rm=TRUE)
table0[4,3] = sum(is.na(b0.is.05))
table0[4,4] = mean(b0.is.optim.05,na.rm=TRUE)
table0[4,5] = mean(b1.is.optim.05,na.rm=TRUE)
table0[4,6] = sum(is.na(b0.is.optim.05))
table0[4,7] = mean(b0.rq.05,na.rm=TRUE)
table0[4,8] = mean(b1.rq.05,na.rm=TRUE)
table0[4,9] = sum(is.na(b0.rq.05))
table0[4,10] = mean(b0.rq.optim.05,na.rm=TRUE)
table0[4,11] = mean(b1.rq.optim.05,na.rm=TRUE)
table0[4,12] = sum(is.na(b0.rq.optim.05))

#### t_0=0 & c=70% ####
b0.is.07 = c()
b1.is.07 = c()
b0.is.optim.07 = c()
b1.is.optim.07 = c()
b0.rq.07 = c()
b1.rq.07 = c()
b0.rq.optim.07 = c()
b1.rq.optim.07 = c()

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.7)
  is.fit = is_est(a[,3], 1, a[,4], a[,5], 0, 0.75)
  b0.is.07[i] = is.fit[1]
  b1.is.07[i] = is.fit[2]
  is.optim.fit = is_optim_est(a[,3], 1, a[,4], a[,5], 0, 0.75)
  b0.is.optim.07[i] = is.optim.fit[1]
  b1.is.optim.07[i] = is.optim.fit[2]
  rq.fit = rq_est(a[,3], 1, a[,4], a[,5], 0, 0.75)
  b0.rq.07[i] = rq.fit[1]
  b1.rq.07[i] = rq.fit[2]
  rq.optim.fit = rq_optim_est(a[,3], 1, a[,4], a[,5], 0, 0.75)
  b0.rq.optim.07[i] = rq.optim.fit[1]
  b1.rq.optim.07[i] = rq.optim.fit[2]
}
table0[5,1] = mean(b0.is.07,na.rm=TRUE)
table0[5,2] = mean(b1.is.07,na.rm=TRUE)
table0[5,3] = sum(is.na(b0.is.07))
table0[5,4] = mean(b0.is.optim.07,na.rm=TRUE)
table0[5,5] = mean(b1.is.optim.07,na.rm=TRUE)
table0[5,6] = sum(is.na(b0.is.optim.07))
table0[5,7] = mean(b0.rq.07,na.rm=TRUE)
table0[5,8] = mean(b1.rq.07,na.rm=TRUE)
table0[5,9] = sum(is.na(b0.rq.07))
table0[5,10] = mean(b0.rq.optim.07,na.rm=TRUE)
table0[5,11] = mean(b1.rq.optim.07,na.rm=TRUE)
table0[5,12] = sum(is.na(b0.rq.optim.07))

#### censoring point at t_0=1 ####
c.0=5000000
c.1=70.39
c.3=24.35
c.5=14.07
c.7=8.49

#### t_0=1 & c=0% ####
b0.is.10 = c()
b1.is.10 = c()
b0.is.optim.10 = c()
b1.is.optim.10 = c()
b0.rq.10 = c()
b1.rq.10 = c()
b0.rq.optim.10 = c()
b1.rq.optim.10 = c()

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.0)
  is.fit = is_est(a[,3], 1, a[,4], a[,5],1, 0.75)
  b0.is.10[i] = is.fit[1]
  b1.is.10[i] = is.fit[2]
  is.optim.fit = is_optim_est(a[,3], 1, a[,4], a[,5], 1, 0.75)
  b0.is.optim.10[i] = is.optim.fit[1]
  b1.is.optim.10[i] = is.optim.fit[2]
  rq.fit = rq_est(a[,3], 1, a[,4], a[,5], 1, 0.75)
  b0.rq.10[i] = rq.fit[1]
  b1.rq.10[i] = rq.fit[2]
  rq.optim.fit = rq_optim_est(a[,3], 1, a[,4], a[,5], 1, 0.75)
  b0.rq.optim.10[i] = rq.optim.fit[1]
  b1.rq.optim.10[i] = rq.optim.fit[2]
}
table1[1,1] = mean(b0.is.10,na.rm=TRUE)
table1[1,2] = mean(b1.is.10,na.rm=TRUE)
table1[1,3] = sum(is.na(b0.is.10))
table1[1,4] = mean(b0.is.optim.10,na.rm=TRUE)
table1[1,5] = mean(b1.is.optim.10,na.rm=TRUE)
table1[1,6] = sum(is.na(b0.is.optim.10))
table1[1,7] = mean(b0.rq.10,na.rm=TRUE)
table1[1,8] = mean(b1.rq.10,na.rm=TRUE)
table1[1,9] = sum(is.na(b0.rq.10))
table1[1,10] = mean(b0.rq.optim.10,na.rm=TRUE)
table1[1,11] = mean(b1.rq.optim.10,na.rm=TRUE)
table1[1,12] = sum(is.na(b0.rq.optim.10))

#### t_0=1 & c=10% ####
b0.is.11 = c()
b1.is.11 = c()
b0.is.optim.11 = c()
b1.is.optim.11 = c()
b0.rq.11 = c()
b1.rq.11 = c()
b0.rq.optim.11 = c()
b1.rq.optim.11 = c()

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.1)
  is.fit = is_est(a[,3], 1, a[,4], a[,5], 1, 0.75)
  b0.is.11[i] = is.fit[1]
  b1.is.11[i] = is.fit[2]
  is.optim.fit = is_optim_est(a[,3], 1, a[,4], a[,5], 1, 0.75)
  b0.is.optim.11[i] = is.optim.fit[1]
  b1.is.optim.11[i] = is.optim.fit[2]
  rq.fit = rq_est(a[,3], 1, a[,4], a[,5], 1, 0.75)
  b0.rq.11[i] = rq.fit[1]
  b1.rq.11[i] = rq.fit[2]
  rq.optim.fit = rq_optim_est(a[,3], 1, a[,4], a[,5], 1, 0.75)
  b0.rq.optim.11[i] = rq.optim.fit[1]
  b1.rq.optim.11[i] = rq.optim.fit[2]
}
table1[2,1] = mean(b0.is.11,na.rm=TRUE)
table1[2,2] = mean(b1.is.11,na.rm=TRUE)
table1[2,3] = sum(is.na(b0.is.11))
table1[2,4] = mean(b0.is.optim.11,na.rm=TRUE)
table1[2,5] = mean(b1.is.optim.11,na.rm=TRUE)
table1[2,6] = sum(is.na(b0.is.optim.11))
table1[2,7] = mean(b0.rq.11,na.rm=TRUE)
table1[2,8] = mean(b1.rq.11,na.rm=TRUE)
table1[2,9] = sum(is.na(b0.rq.11))
table1[2,10] = mean(b0.rq.optim.11,na.rm=TRUE)
table1[2,11] = mean(b1.rq.optim.11,na.rm=TRUE)
table1[2,12] = sum(is.na(b0.rq.optim.11))


#### t_0=1 & c=30% ####
b0.is.13 = c()
b1.is.13 = c()
b0.is.optim.13 = c()
b1.is.optim.13 = c()
b0.rq.13 = c()
b1.rq.13 = c()
b0.rq.optim.13 = c()
b1.rq.optim.13 = c()

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.3)
  is.fit = is_est(a[,3], 1, a[,4], a[,5], 1, 0.75)
  b0.is.13[i] = is.fit[1]
  b1.is.13[i] = is.fit[2]
  is.optim.fit = is_optim_est(a[,3], 1, a[,4], a[,5], 1, 0.75)
  b0.is.optim.13[i] = is.optim.fit[1]
  b1.is.optim.13[i] = is.optim.fit[2]
  rq.fit = rq_est(a[,3], 1, a[,4], a[,5], 1, 0.75)
  b0.rq.13[i] = rq.fit[1]
  b1.rq.13[i] = rq.fit[2]
  rq.optim.fit = rq_optim_est(a[,3], 1, a[,4], a[,5], 1, 0.75)
  b0.rq.optim.13[i] = rq.optim.fit[1]
  b1.rq.optim.13[i] = rq.optim.fit[2]
}
table1[3,1] = mean(b0.is.13,na.rm=TRUE)
table1[3,2] = mean(b1.is.13,na.rm=TRUE)
table1[3,3] = sum(is.na(b0.is.13))
table1[3,4] = mean(b0.is.optim.13,na.rm=TRUE)
table1[3,5] = mean(b1.is.optim.13,na.rm=TRUE)
table1[3,6] = sum(is.na(b0.is.optim.13))
table1[3,7] = mean(b0.rq.13,na.rm=TRUE)
table1[3,8] = mean(b1.rq.13,na.rm=TRUE)
table1[3,9] = sum(is.na(b0.rq.13))
table1[3,10] = mean(b0.rq.optim.13,na.rm=TRUE)
table1[3,11] = mean(b1.rq.optim.13,na.rm=TRUE)
table1[3,12] = sum(is.na(b0.rq.optim.13))

#### t_0=1 & c=50% ####
b0.is.15 = c()
b1.is.15 = c()
b0.is.optim.15 = c()
b1.is.optim.15 = c()
b0.rq.15 = c()
b1.rq.15 = c()
b0.rq.optim.15 = c()
b1.rq.optim.15 = c()

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.5)
  is.fit = is_est(a[,3], 1, a[,4], a[,5], 1, 0.75)
  b0.is.15[i] = is.fit[1]
  b1.is.15[i] = is.fit[2]
  is.optim.fit = is_optim_est(a[,3], 1, a[,4], a[,5], 1, 0.75)
  b0.is.optim.15[i] = is.optim.fit[1]
  b1.is.optim.15[i] = is.optim.fit[2]
  rq.fit = rq_est(a[,3], 1, a[,4], a[,5], 1, 0.75)
  b0.rq.15[i] = rq.fit[1]
  b1.rq.15[i] = rq.fit[2]
  rq.optim.fit = rq_optim_est(a[,3], 1, a[,4], a[,5], 1, 0.75)
  b0.rq.optim.15[i] = rq.optim.fit[1]
  b1.rq.optim.15[i] = rq.optim.fit[2]
}
table1[4,1] = mean(b0.is.15,na.rm=TRUE)
table1[4,2] = mean(b1.is.15,na.rm=TRUE)
table1[4,3] = sum(is.na(b0.is.15))
table1[4,4] = mean(b0.is.optim.15,na.rm=TRUE)
table1[4,5] = mean(b1.is.optim.15,na.rm=TRUE)
table1[4,6] = sum(is.na(b0.is.optim.15))
table1[4,7] = mean(b0.rq.15,na.rm=TRUE)
table1[4,8] = mean(b1.rq.15,na.rm=TRUE)
table1[4,9] = sum(is.na(b0.rq.15))
table1[4,10] = mean(b0.rq.optim.15,na.rm=TRUE)
table1[4,11] = mean(b1.rq.optim.15,na.rm=TRUE)
table1[4,12] = sum(is.na(b0.rq.optim.15))

#### t_0=1 & c=70% ####
b0.is.17 = c()
b1.is.17 = c()
b0.is.optim.17 = c()
b1.is.optim.17 = c()
b0.rq.17 = c()
b1.rq.17 = c()
b0.rq.optim.17 = c()
b1.rq.optim.17 = c()

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.7)
  is.fit = is_est(a[,3], 1, a[,4], a[,5], 1, 0.75)
  b0.is.17[i] = is.fit[1]
  b1.is.17[i] = is.fit[2]
  is.optim.fit = is_optim_est(a[,3], 1, a[,4], a[,5], 1, 0.75)
  b0.is.optim.17[i] = is.optim.fit[1]
  b1.is.optim.17[i] = is.optim.fit[2]
  rq.fit = rq_est(a[,3], 1, a[,4], a[,5], 1, 0.75)
  b0.rq.17[i] = rq.fit[1]
  b1.rq.17[i] = rq.fit[2]
  rq.optim.fit = rq_optim_est(a[,3], 1, a[,4], a[,5], 1, 0.75)
  b0.rq.optim.17[i] = rq.optim.fit[1]
  b1.rq.optim.17[i] = rq.optim.fit[2]
}
table1[5,1] = mean(b0.is.17,na.rm=TRUE)
table1[5,2] = mean(b1.is.17,na.rm=TRUE)
table1[5,3] = sum(is.na(b0.is.17))
table1[5,4] = mean(b0.is.optim.17,na.rm=TRUE)
table1[5,5] = mean(b1.is.optim.17,na.rm=TRUE)
table1[5,6] = sum(is.na(b0.is.optim.17))
table1[5,7] = mean(b0.rq.17,na.rm=TRUE)
table1[5,8] = mean(b1.rq.17,na.rm=TRUE)
table1[5,9] = sum(is.na(b0.rq.17))
table1[5,10] = mean(b0.rq.optim.17,na.rm=TRUE)
table1[5,11] = mean(b1.rq.optim.17,na.rm=TRUE)
table1[5,12] = sum(is.na(b0.rq.optim.17))

#### censoring point at t_0=2 ####
c.0=5000000
c.1=64.86
c.3=23.34
c.5=13.62
c.7=8.36
#### t_0=2 & c=0% ####
b0.is.20 = c()
b1.is.20 = c()
b0.is.optim.20 = c()
b1.is.optim.20 = c()
b0.rq.20 = c()
b1.rq.20 = c()
b0.rq.optim.20 = c()
b1.rq.optim.20 = c()

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.0)
  is.fit = is_est(a[,3], 1, a[,4], a[,5], 2, 0.75)
  b0.is.20[i] = is.fit[1]
  b1.is.20[i] = is.fit[2]
  is.optim.fit = is_optim_est(a[,3], 1, a[,4], a[,5], 2, 0.75)
  b0.is.optim.20[i] = is.optim.fit[1]
  b1.is.optim.20[i] = is.optim.fit[2]
  rq.fit = rq_est(a[,3], 1, a[,4], a[,5], 2, 0.75)
  b0.rq.20[i] = rq.fit[1]
  b1.rq.20[i] = rq.fit[2]
  rq.optim.fit = rq_optim_est(a[,3], 1, a[,4], a[,5], 2, 0.75)
  b0.rq.optim.20[i] = rq.optim.fit[1]
  b1.rq.optim.20[i] = rq.optim.fit[2]
}
table2[1,1] = mean(b0.is.20,na.rm=TRUE)
table2[1,2] = mean(b1.is.20,na.rm=TRUE)
table2[1,3] = sum(is.na(b0.is.20))
table2[1,4] = mean(b0.is.optim.20,na.rm=TRUE)
table2[1,5] = mean(b1.is.optim.20,na.rm=TRUE)
table2[1,6] = sum(is.na(b0.is.optim.20))
table2[1,7] = mean(b0.rq.20,na.rm=TRUE)
table2[1,8] = mean(b1.rq.20,na.rm=TRUE)
table2[1,9] = sum(is.na(b0.rq.20))
table2[1,10] = mean(b0.rq.optim.20,na.rm=TRUE)
table2[1,11] = mean(b1.rq.optim.20,na.rm=TRUE)
table2[1,12] = sum(is.na(b0.rq.optim.20))

#### t_0=2 & c=10% ####
b0.is.21 = c()
b1.is.21 = c()
b0.is.optim.21 = c()
b1.is.optim.21 = c()
b0.rq.21 = c()
b1.rq.21 = c()
b0.rq.optim.21 = c()
b1.rq.optim.21 = c()

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.1)
  is.fit = is_est(a[,3], 1, a[,4], a[,5], 2, 0.75)
  b0.is.21[i] = is.fit[1]
  b1.is.21[i] = is.fit[2]
  is.optim.fit = is_optim_est(a[,3], 1, a[,4], a[,5], 2, 0.75)
  b0.is.optim.21[i] = is.optim.fit[1]
  b1.is.optim.21[i] = is.optim.fit[2]
  rq.fit = rq_est(a[,3], 1, a[,4], a[,5], 2, 0.75)
  b0.rq.21[i] = rq.fit[1]
  b1.rq.21[i] = rq.fit[2]
  rq.optim.fit = rq_optim_est(a[,3], 1, a[,4], a[,5], 2, 0.75)
  b0.rq.optim.21[i] = rq.optim.fit[1]
  b1.rq.optim.21[i] = rq.optim.fit[2]
}
table2[2,1] = mean(b0.is.21,na.rm=TRUE)
table2[2,2] = mean(b1.is.21,na.rm=TRUE)
table2[2,3] = sum(is.na(b0.is.21))
table2[2,4] = mean(b0.is.optim.21,na.rm=TRUE)
table2[2,5] = mean(b1.is.optim.21,na.rm=TRUE)
table2[2,6] = sum(is.na(b0.is.optim.21))
table2[2,7] = mean(b0.rq.21,na.rm=TRUE)
table2[2,8] = mean(b1.rq.21,na.rm=TRUE)
table2[2,9] = sum(is.na(b0.rq.21))
table2[2,10] = mean(b0.rq.optim.21,na.rm=TRUE)
table2[2,11] = mean(b1.rq.optim.21,na.rm=TRUE)
table2[2,12] = sum(is.na(b0.rq.optim.21))


#### t_0=2 & c=30% ####
b0.is.23 = c()
b1.is.23 = c()
b0.is.optim.23 = c()
b1.is.optim.23 = c()
b0.rq.23 = c()
b1.rq.23 = c()
b0.rq.optim.23 = c()
b1.rq.optim.23 = c()

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.3)
  is.fit = is_est(a[,3], 1, a[,4], a[,5], 2, 0.75)
  b0.is.23[i] = is.fit[1]
  b1.is.23[i] = is.fit[2]
  is.optim.fit = is_optim_est(a[,3], 1, a[,4], a[,5], 2, 0.75)
  b0.is.optim.23[i] = is.optim.fit[1]
  b1.is.optim.23[i] = is.optim.fit[2]
  rq.fit = rq_est(a[,3], 1, a[,4], a[,5], 2, 0.75)
  b0.rq.23[i] = rq.fit[1]
  b1.rq.23[i] = rq.fit[2]
  rq.optim.fit = rq_optim_est(a[,3], 1, a[,4], a[,5], 2, 0.75)
  b0.rq.optim.23[i] = rq.optim.fit[1]
  b1.rq.optim.23[i] = rq.optim.fit[2]
}
table2[3,1] = mean(b0.is.23,na.rm=TRUE)
table2[3,2] = mean(b1.is.23,na.rm=TRUE)
table2[3,3] = sum(is.na(b0.is.23))
table2[3,4] = mean(b0.is.optim.23,na.rm=TRUE)
table2[3,5] = mean(b1.is.optim.23,na.rm=TRUE)
table2[3,6] = sum(is.na(b0.is.optim.23))
table2[3,7] = mean(b0.rq.23,na.rm=TRUE)
table2[3,8] = mean(b1.rq.23,na.rm=TRUE)
table2[3,9] = sum(is.na(b0.rq.23))
table2[3,10] = mean(b0.rq.optim.23,na.rm=TRUE)
table2[3,11] = mean(b1.rq.optim.23,na.rm=TRUE)
table2[3,12] = sum(is.na(b0.rq.optim.23))

#### t_0=2 & c=50% ####
b0.is.25 = c()
b1.is.25 = c()
b0.is.optim.25 = c()
b1.is.optim.25 = c()
b0.rq.25 = c()
b1.rq.25 = c()
b0.rq.optim.25 = c()
b1.rq.optim.25 = c()

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.5)
  is.fit = is_est(a[,3], 1, a[,4], a[,5], 2, 0.75)
  b0.is.25[i] = is.fit[1]
  b1.is.25[i] = is.fit[2]
  is.optim.fit = is_optim_est(a[,3], 1, a[,4], a[,5], 2, 0.75)
  b0.is.optim.25[i] = is.optim.fit[1]
  b1.is.optim.25[i] = is.optim.fit[2]
  rq.fit = rq_est(a[,3], 1, a[,4], a[,5], 2, 0.75)
  b0.rq.25[i] = rq.fit[1]
  b1.rq.25[i] = rq.fit[2]
  rq.optim.fit = rq_optim_est(a[,3], 1, a[,4], a[,5], 2, 0.75)
  b0.rq.optim.25[i] = rq.optim.fit[1]
  b1.rq.optim.25[i] = rq.optim.fit[2]
}
table2[4,1] = mean(b0.is.25,na.rm=TRUE)
table2[4,2] = mean(b1.is.25,na.rm=TRUE)
table2[4,3] = sum(is.na(b0.is.25))
table2[4,4] = mean(b0.is.optim.25,na.rm=TRUE)
table2[4,5] = mean(b1.is.optim.25,na.rm=TRUE)
table2[4,6] = sum(is.na(b0.is.optim.25))
table2[4,7] = mean(b0.rq.25,na.rm=TRUE)
table2[4,8] = mean(b1.rq.25,na.rm=TRUE)
table2[4,9] = sum(is.na(b0.rq.25))
table2[4,10] = mean(b0.rq.optim.25,na.rm=TRUE)
table2[4,11] = mean(b1.rq.optim.25,na.rm=TRUE)
table2[4,12] = sum(is.na(b0.rq.optim.25))

#### t_0=2 & c=70% ####
b0.is.27 = c()
b1.is.27 = c()
b0.is.optim.27 = c()
b1.is.optim.27 = c()
b0.rq.27 = c()
b1.rq.27 = c()
b0.rq.optim.27 = c()
b1.rq.optim.27 = c()

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.7)
  is.fit = is_est(a[,3], 1, a[,4], a[,5], 2, 0.75)
  b0.is.27[i] = is.fit[1]
  b1.is.27[i] = is.fit[2]
  is.optim.fit = is_optim_est(a[,3], 1, a[,4], a[,5], 2, 0.75)
  b0.is.optim.27[i] = is.optim.fit[1]
  b1.is.optim.27[i] = is.optim.fit[2]
  rq.fit = rq_est(a[,3], 1, a[,4], a[,5], 2, 0.75)
  b0.rq.27[i] = rq.fit[1]
  b1.rq.27[i] = rq.fit[2]
  rq.optim.fit = rq_optim_est(a[,3], 1, a[,4], a[,5], 2, 0.75)
  b0.rq.optim.27[i] = rq.optim.fit[1]
  b1.rq.optim.27[i] = rq.optim.fit[2]
}
table2[5,1] = mean(b0.is.27,na.rm=TRUE)
table2[5,2] = mean(b1.is.27,na.rm=TRUE)
table2[5,3] = sum(is.na(b0.is.27))
table2[5,4] = mean(b0.is.optim.27,na.rm=TRUE)
table2[5,5] = mean(b1.is.optim.27,na.rm=TRUE)
table2[5,6] = sum(is.na(b0.is.optim.27))
table2[5,7] = mean(b0.rq.27,na.rm=TRUE)
table2[5,8] = mean(b1.rq.27,na.rm=TRUE)
table2[5,9] = sum(is.na(b0.rq.27))
table2[5,10] = mean(b0.rq.optim.27,na.rm=TRUE)
table2[5,11] = mean(b1.rq.optim.27,na.rm=TRUE)
table2[5,12] = sum(is.na(b0.rq.optim.27))

#### censoring point at t_0=3 ####
c.0=5000000
c.1=61.17
c.3=22.53
c.5=13.43
c.7=8.5
#### t_0=3 & c=0% ####
b0.is.30 = c()
b1.is.30 = c()
b0.is.optim.30 = c()
b1.is.optim.30 = c()
b0.rq.30 = c()
b1.rq.30 = c()
b0.rq.optim.30 = c()
b1.rq.optim.30 = c()

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.0)
  is.fit = is_est(a[,3], 1, a[,4], a[,5], 3, 0.75)
  b0.is.30[i] = is.fit[1]
  b1.is.30[i] = is.fit[2]
  is.optim.fit = is_optim_est(a[,3], 1, a[,4], a[,5], 3, 0.75)
  b0.is.optim.30[i] = is.optim.fit[1]
  b1.is.optim.30[i] = is.optim.fit[2]
  rq.fit = rq_est(a[,3], 1, a[,4], a[,5], 3, 0.75)
  b0.rq.30[i] = rq.fit[1]
  b1.rq.30[i] = rq.fit[2]
  rq.optim.fit = rq_optim_est(a[,3], 1, a[,4], a[,5], 3, 0.75)
  b0.rq.optim.30[i] = rq.optim.fit[1]
  b1.rq.optim.30[i] = rq.optim.fit[2]
}
table3[1,1] = mean(b0.is.30,na.rm=TRUE)
table3[1,2] = mean(b1.is.30,na.rm=TRUE)
table3[1,3] = sum(is.na(b0.is.30))
table3[1,4] = mean(b0.is.optim.30,na.rm=TRUE)
table3[1,5] = mean(b1.is.optim.30,na.rm=TRUE)
table3[1,6] = sum(is.na(b0.is.optim.30))
table3[1,7] = mean(b0.rq.30,na.rm=TRUE)
table3[1,8] = mean(b1.rq.30,na.rm=TRUE)
table3[1,9] = sum(is.na(b0.rq.30))
table3[1,10] = mean(b0.rq.optim.30,na.rm=TRUE)
table3[1,11] = mean(b1.rq.optim.30,na.rm=TRUE)
table3[1,12] = sum(is.na(b0.rq.optim.30))

#### t_0=3 & c=10% ####
b0.is.31 = c()
b1.is.31 = c()
b0.is.optim.31 = c()
b1.is.optim.31 = c()
b0.rq.31 = c()
b1.rq.31 = c()
b0.rq.optim.31 = c()
b1.rq.optim.31 = c()

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.1)
  is.fit = is_est(a[,3], 1, a[,4], a[,5], 3, 0.75)
  b0.is.31[i] = is.fit[1]
  b1.is.31[i] = is.fit[2]
  is.optim.fit = is_optim_est(a[,3], 1, a[,4], a[,5], 3, 0.75)
  b0.is.optim.31[i] = is.optim.fit[1]
  b1.is.optim.31[i] = is.optim.fit[2]
  rq.fit = rq_est(a[,3], 1, a[,4], a[,5], 3, 0.75)
  b0.rq.31[i] = rq.fit[1]
  b1.rq.31[i] = rq.fit[2]
  rq.optim.fit = rq_optim_est(a[,3], 1, a[,4], a[,5], 3, 0.75)
  b0.rq.optim.31[i] = rq.optim.fit[1]
  b1.rq.optim.31[i] = rq.optim.fit[2]
}
table3[2,1] = mean(b0.is.31,na.rm=TRUE)
table3[2,2] = mean(b1.is.31,na.rm=TRUE)
table3[2,3] = sum(is.na(b0.is.31))
table3[2,4] = mean(b0.is.optim.31,na.rm=TRUE)
table3[2,5] = mean(b1.is.optim.31,na.rm=TRUE)
table3[2,6] = sum(is.na(b0.is.optim.31))
table3[2,7] = mean(b0.rq.31,na.rm=TRUE)
table3[2,8] = mean(b1.rq.31,na.rm=TRUE)
table3[2,9] = sum(is.na(b0.rq.31))
table3[2,10] = mean(b0.rq.optim.31,na.rm=TRUE)
table3[2,11] = mean(b1.rq.optim.31,na.rm=TRUE)
table3[2,12] = sum(is.na(b0.rq.optim.31))


#### t_0=3 & c=30% ####
b0.is.33 = c()
b1.is.33 = c()
b0.is.optim.33 = c()
b1.is.optim.33 = c()
b0.rq.33 = c()
b1.rq.33 = c()
b0.rq.optim.33 = c()
b1.rq.optim.33 = c()

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.3)
  is.fit = is_est(a[,3], 1, a[,4], a[,5], 3, 0.75)
  b0.is.33[i] = is.fit[1]
  b1.is.33[i] = is.fit[2]
  is.optim.fit = is_optim_est(a[,3], 1, a[,4], a[,5], 3, 0.75)
  b0.is.optim.33[i] = is.optim.fit[1]
  b1.is.optim.33[i] = is.optim.fit[2]
  rq.fit = rq_est(a[,3], 1, a[,4], a[,5], 3, 0.75)
  b0.rq.33[i] = rq.fit[1]
  b1.rq.33[i] = rq.fit[2]
  rq.optim.fit = rq_optim_est(a[,3], 1, a[,4], a[,5], 3, 0.75)
  b0.rq.optim.33[i] = rq.optim.fit[1]
  b1.rq.optim.33[i] = rq.optim.fit[2]
}
table3[3,1] = mean(b0.is.33,na.rm=TRUE)
table3[3,2] = mean(b1.is.33,na.rm=TRUE)
table3[3,3] = sum(is.na(b0.is.33))
table3[3,4] = mean(b0.is.optim.33,na.rm=TRUE)
table3[3,5] = mean(b1.is.optim.33,na.rm=TRUE)
table3[3,6] = sum(is.na(b0.is.optim.33))
table3[3,7] = mean(b0.rq.33,na.rm=TRUE)
table3[3,8] = mean(b1.rq.33,na.rm=TRUE)
table3[3,9] = sum(is.na(b0.rq.33))
table3[3,10] = mean(b0.rq.optim.33,na.rm=TRUE)
table3[3,11] = mean(b1.rq.optim.33,na.rm=TRUE)
table3[3,12] = sum(is.na(b0.rq.optim.33))

#### t_0=3 & c=50% ####
b0.is.35 = c()
b1.is.35 = c()
b0.is.optim.35 = c()
b1.is.optim.35 = c()
b0.rq.35 = c()
b1.rq.35 = c()
b0.rq.optim.35 = c()
b1.rq.optim.35 = c()

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.5)
  is.fit = is_est(a[,3], 1, a[,4], a[,5], 3, 0.75)
  b0.is.35[i] = is.fit[1]
  b1.is.35[i] = is.fit[2]
  is.optim.fit = is_optim_est(a[,3], 1, a[,4], a[,5], 3, 0.75)
  b0.is.optim.35[i] = is.optim.fit[1]
  b1.is.optim.35[i] = is.optim.fit[2]
  rq.fit = rq_est(a[,3], 1, a[,4], a[,5], 3, 0.75)
  b0.rq.35[i] = rq.fit[1]
  b1.rq.35[i] = rq.fit[2]
  rq.optim.fit = rq_optim_est(a[,3], 1, a[,4], a[,5], 3, 0.75)
  b0.rq.optim.35[i] = rq.optim.fit[1]
  b1.rq.optim.35[i] = rq.optim.fit[2]
}
table3[4,1] = mean(b0.is.35,na.rm=TRUE)
table3[4,2] = mean(b1.is.35,na.rm=TRUE)
table3[4,3] = sum(is.na(b0.is.35))
table3[4,4] = mean(b0.is.optim.35,na.rm=TRUE)
table3[4,5] = mean(b1.is.optim.35,na.rm=TRUE)
table3[4,6] = sum(is.na(b0.is.optim.35))
table3[4,7] = mean(b0.rq.35,na.rm=TRUE)
table3[4,8] = mean(b1.rq.35,na.rm=TRUE)
table3[4,9] = sum(is.na(b0.rq.35))
table3[4,10] = mean(b0.rq.optim.35,na.rm=TRUE)
table3[4,11] = mean(b1.rq.optim.35,na.rm=TRUE)
table3[4,12] = sum(is.na(b0.rq.optim.35))

#### t_0=3 & c=70% ####
b0.is.37 = c()
b1.is.37 = c()
b0.is.optim.37 = c()
b1.is.optim.37 = c()
b0.rq.37 = c()
b1.rq.37 = c()
b0.rq.optim.37 = c()
b1.rq.optim.37 = c()

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.7)
  is.fit = is_est(a[,3], 1, a[,4], a[,5], 3, 0.75)
  b0.is.37[i] = is.fit[1]
  b1.is.37[i] = is.fit[2]
  is.optim.fit = is_optim_est(a[,3], 1, a[,4], a[,5], 3, 0.75)
  b0.is.optim.37[i] = is.optim.fit[1]
  b1.is.optim.37[i] = is.optim.fit[2]
  rq.fit = rq_est(a[,3], 1, a[,4], a[,5], 3, 0.75)
  b0.rq.37[i] = rq.fit[1]
  b1.rq.37[i] = rq.fit[2]
  rq.optim.fit = rq_optim_est(a[,3], 1, a[,4], a[,5], 3, 0.75)
  b0.rq.optim.37[i] = rq.optim.fit[1]
  b1.rq.optim.37[i] = rq.optim.fit[2]
}
table3[5,1] = mean(b0.is.37,na.rm=TRUE)
table3[5,2] = mean(b1.is.37,na.rm=TRUE)
table3[5,3] = sum(is.na(b0.is.37))
table3[5,4] = mean(b0.is.optim.37,na.rm=TRUE)
table3[5,5] = mean(b1.is.optim.37,na.rm=TRUE)
table3[5,6] = sum(is.na(b0.is.optim.37))
table3[5,7] = mean(b0.rq.37,na.rm=TRUE)
table3[5,8] = mean(b1.rq.37,na.rm=TRUE)
table3[5,9] = sum(is.na(b0.rq.37))
table3[5,10] = mean(b0.rq.optim.37,na.rm=TRUE)
table3[5,11] = mean(b1.rq.optim.37,na.rm=TRUE)
table3[5,12] = sum(is.na(b0.rq.optim.37))
