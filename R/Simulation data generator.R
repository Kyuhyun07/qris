## Simulation data generator
######## Data Generation function ########
## We refer simulation setting from Jung et al. 2009
## samplesize = size of dataset
## We generate censoring time by uniform(0, censor)
## And we empirically find censor values to keep censoring proportion 0%~70%
data.gen<-function(samplesize, censor){
    sim=matrix(NA,samplesize,5)
    colnames(sim) = c("T","C","Z","X","censored")
    ## Generate C_i
    sim[,2] = runif(samplesize,0,censor)
    ## Covariates (Control=0, Treatment=1)
    sim[,4] = rbinom(samplesize,size=1,p=0.5)
    ## Generate T_i (Given Condition r=rho_0, k=2, exp(beta_0)=5, exp(beta_1)=2))
    unif = runif(n=samplesize ,min = 0,max = 1)
    for (q in 1:samplesize){
        sim[q,1]={{-log(1-unif[q])}^(1/k)}/r.initial.0
    }
    ## Generate Y_i (min(T,C))
    sim[,3] = apply(sim[,1:2], 1, FUN=min)
    ## Censoring indicator (Censored=0, Not censored=1)
    sim[,5]=I(sim[,1]<sim[,2])
    ## Ordering
    sim = sim[order(sim[,3]),]
    n = nrow(sim)
    sim = as.data.frame(sim)
    return(sim)
}


## Assumed necessary parameters
######## ndata1-3, both beta0, 1  effective ########
exp.beta.initial.0=5
exp.beta.initial.1=10
k=2
r.initial.0=(log(10/7.5))^(1/k)/exp.beta.initial.0
r.initial.1=(log(10/7.5))^(1/k)/exp.beta.initial.1

######## ndata4-6, Only beta0 effective ########
exp.beta.initial.0=5
k=2
r.initial.0=(log(10/7.5))^(1/k)/exp.beta.initial.0

######## Censor information ########
######## ndata1 : 25% quantile ########
## t_0 = 0
c.0=5000000
c.1=123.69
c.3=41.27
c.5=23.55
c.7=14.2

## t_0 = 1
c.0=5000000
c.1=115.24
c.3=39.17
c.5=22.5
c.7=13.55

## t_0 = 2
c.0=5000000
c.1=106.73
c.3=37.34
c.5=21.61
c.7=13.07

## t_0 = 3
c.0=5000000
c.1=101.62
c.3=35.91
c.5=21.19
c.7=13.04

######## ndata2 : 50% quantile ########
## t_0 = 0
c.0=5000000
c.1=78.11
c.3=26.36
c.5=15.08
c.7=9.09

## t_0 = 1
c.0=5000000
c.1=70.39
c.3=24.35
c.5=14.07
c.7=8.49

## t_0 = 2
c.0=5000000
c.1=64.86
c.3=23.34
c.5=13.62
c.7=8.36

## t_0 = 3
c.0=5000000
c.1=61.17
c.3=22.53
c.5=13.43
c.7=8.5

######## ndata3 : 75% quantile ########
## t_0 = 0
c.0=5000000
c.1=55.68
c.3=18.79
c.5=10.69
c.7=6.49

## t_0 = 1
c.0=5000000
c.1=48.34
c.3=16.92
c.5=9.86
c.7=5.96

## t_0 = 2
c.0=5000000
c.1=43.86
c.3=16.08
c.5=9.57
c.7=6.03

## t_0 = 3
c.0=5000000
c.1=41.23
c.3=15.8
c.5=9.8
c.7=6.44

######## ndata4 : 25% quantile ########
## t_0 = 0
c.0=5000000
c.1=81.63
c.3=27.42
c.5=16.3
c.7=10.48

## t_0 = 1
c.0=5000000
c.1=73.44
c.3=25.32
c.5=15.34
c.7=9.92

## t_0 = 2
c.0=5000000
c.1=66.69
c.3=23.86
c.5=14.72
c.7=9.64

## t_0 = 3
c.0=5000000
c.1=61.86
c.3=22.66
c.5=14.33
c.7=9.58
######## ndata5 : 50% quantile ########
## t_0 = 0
c.0=5000000
c.1=52.15
c.3=17.76
c.5=10.5
c.7=6.8

## t_0 = 1
c.0=5000000
c.1=44.78
c.3=15.85
c.5=9.63
c.7=6.26

## t_0 = 2
c.0=5000000
c.1=39.03
c.3=14.55
c.5=9.22
c.7=6.2

## t_0 = 3
c.0=5000000
c.1=35
c.3=13.78
c.5=9.06
c.7=6.36

######## ndata6 : 75% quantile ########
## t_0 = 0
c.0=5000000
c.1=36.98
c.3=12.53
c.5=7.47
c.7=4.81

## t_0 = 1
c.0=5000000
c.1=30.1
c.3=10.82
c.5=6.68
c.7=4.61

## t_0 = 2
c.0=5000000
c.1=25.65
c.3=9.89
c.5=6.44
c.7=4.5

## t_0 = 3
c.0=5000000
c.1=21.92
c.3=9.46
c.5=6.59
c.7=4.41


## Example
source("allcodes.R")
a = data.gen(200,c.3)
W = weight_generator(a[,3], 0, 1, a[,4], a[,5])

library(emplik)
library(survival)
head(a)

foo <- WKM(a[,1], 1 - a[,5])
foo2 <- survfit(Surv(T, 1 - censored) ~ 1, data = a)

str(foo)
str(foo2)
