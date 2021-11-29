
<!-- README.md is generated from README.Rmd. Please edit that file -->

# qrismb

<!-- badges: start -->

<!-- badges: end -->

The **qrismb** package implements estimation procedures for a regression
model of the quantiles of residual life, remaining lifetime at a
specific time, subject to right censoring. For estimation of regression
parameters, we consider an induced smoothed method that solves smoothed
weighted estimating equations. We also consider the estimation method
that solves the original non-smooth weighted estimating equations via a
L1 minimization method. To handle data subject to right censoring,
inverse probabilities of censoring are incorporated as weights. For
standard errors estimation, a robust sandwich-type covariance estimator
aided by an efficient resampling method, and a full multiplier bootstrap
approach are considered for the induced smoothed estimator (“smooth”)
and non-smooth estimator (“nonsmooth”), respectively. Furthermore, an
iterative procedure that simultaneously estimates regression parameters
and their standard errors is implemented.

## Installation

You can install the released version of qrismb from
[GitHub](https://github.com/Kyuhyun07/qrismb) with:

``` r
> # install.packages("devtools")
> devtools::install_github("Kyuhyun07/qrismb")
> library(qrismb)
```

## Example

There are two examples

``` r
> data.gen <- function(n) {
+     r0 <- .2 * sqrt(log(2))
+     r1 <- .1 * sqrt(log(2))
+     dat <- data.frame(censoring = runif(n, 0, 24.35),
+                       Time0 = sqrt(-log(1 - runif(n))),
+                       X = rbinom(n, 1, .5))
+     dat$Time0 <- ifelse(dat$X > 0, dat$Time0 / r1, dat$Time0 / r0)
+     dat$Time <- pmin(dat$Time0, dat$censoring)
+     dat$status <- 1 * (dat$Time0 < dat$censoring)
+     subset(dat, select = c(Time, status, X))
+ }
> 
> set.seed(1)
> dat <- data.gen(200)
> 
> library(qrismb)
> fm <- Surv(Time, status) ~ X
> fit1 <- qrismb(fm, data = dat, t0 = 1, Q = 0.5, ne = 200, "smooth", "rq")
> fit2 <- qrismb(fm, data = dat, t0 = 1, Q = 0.5, ne = 200, "nonsmooth", "noeffect")
> fit3 <- qrismb(fm, data = dat, t0 = 1, Q = 0.5, ne = 200, "iterative", "userinput", c(2,1))
> 
> coef(fit1)
(Intercept)           X 
   1.252776    0.810027 
> summary(fit2)
Call:
qrismb(formula = fm, data = dat, t0 = 1, Q = 0.5, ne = 200, method = "nonsmooth", 
    init = "noeffect")

qrismb Estimator
            estimate std.Error z.value   p.value    
(Intercept)   1.2395    0.0861  14.401 < 2.2e-16 ***
X             0.8525    0.1191   7.157 < 2.2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
> vcov(fit3)
            (Intercept)           X
(Intercept)  0.01722485 -0.07026059
X           -0.07026059 27.25328957
```

``` r
> ## 2. real data example
> ## Load "retinopathy" data from R survival package
> library(survival)
Warning: package 'survival' was built under R version 3.6.2
> ## Real data application
> data(cancer, package = "survival")
> lung2 <- subset(lung, select = c(time, status, age, sex))
> ## tidy up the data
> lung2$status <- lung2$status - 1
> lung2$sex <- lung2$sex - 1
> 
> library(qrismb)
> set.seed(1)
> fm <- Surv(time, status) ~ age + sex
> fit1 <- qrismb(fm, data = lung2, t0 = 0, Q = 0.5, ne = 200, "iterative")
> fit2 <- qrismb(fm, data = lung2, t0 = 30, Q = 0.5, ne = 200, "nonsmooth", "userinput", c(1,0,1))
> fit3 <- qrismb(fm, data = lung2, t0 = 100, Q = 0.5, ne = 200,"smooth", "noeffect")
> 
> coef(fit1)
(Intercept)         age         sex 
5.206477586 0.006945976 0.182889952 
> summary(fit2)
Call:
qrismb(formula = fm, data = lung2, t0 = 30, Q = 0.5, ne = 200, 
    method = "nonsmooth", init = "userinput", userinit = c(1, 
        0, 1))

qrismb Estimator
            estimate std.Error z.value p.value    
(Intercept)   9.1260    2.2833   3.997  0.0001 ***
age          -0.0654    0.0363  -1.805  0.0711 .  
sex           1.7389    0.6270   2.773  0.0055 ** 
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
> vcov(fit3)
            (Intercept)           age           sex
(Intercept)  0.76950264 -0.0112458777  0.0200614670
age         -0.01124588  0.0001669553 -0.0004766159
sex          0.02006147 -0.0004766159  0.0434881028
```

## Reference

Chiou, S., Kang, S., and Yan, J. (2014). Fitting accelerated failure
time model in routine survival analysis with R package aftgee. *Journal
of Statistical Software*, **61**(11): 1–23.

Li, R., Huang, X., & Cortes, J. (2016). Quantile residual life
regression with longitudinal biomarker measurements for dynamic
prediction. *Journal of the Royal Statistical Society*. **Series C**
(Applied Statistics), 755-773.

Jung, S. H., Jeong, J. H., & Bandos, H. (2009). Regression on quantile
residual life. *Biometrics*, **65**(4), 1203-1212.
