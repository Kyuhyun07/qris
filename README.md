
<!-- README.md is generated from README.Rmd. Please edit that file -->

# qrismb

<!-- badges: start -->

<!-- badges: end -->

The **qrismb** package implements a regression modeling of the quantiles
of residual life, remaining lifetime at a specific time. For estimation
of regression parameters, we applies an induced smoothed version of the
existing non-smooth estimating equations approaches, and also uses
L1-minimization approach as well. Furthermore, iterative procedure is
implemented. For estimation of regression parameter, “smooth” method
uses a robust sandwich-type covariance estimator of regression
estimators with resampling method, and “nonsmooth” method uses full
multiplier bootstrap approach. To handle data subject to right
censoring, inverse probabilities of censoring are incorporated as
weights.

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
+   r0 <- .2 * sqrt(log(2))
+   r1 <- .1 * sqrt(log(2))
+   dat <- data.frame(censoring = runif(n, 0, 24.35),
+                     Time0 = sqrt(-log(1 - runif(n))),
+                     X = rbinom(n, 1, .5))
+   dat$Time0 <- ifelse(dat$X > 0, dat$Time0 / r1, dat$Time0 / r0)
+   dat$Time <- pmin(dat$Time0, dat$censoring)
+   dat$status <- 1 * (dat$Time0 < dat$censoring)
+   subset(dat, select = c(Time, status, X))
+ }
> set.seed(1)
> dat <- data.gen(200)
> 
> library(qrismb)
> fm <- Surv(Time, status) ~ X
> fit1 <- qrismb(fm, data = dat, t0 = 1, Q = 0.5, ne = 200, "rq", "smooth")
> fit2 <- qrismb(fm, data = dat, t0 = 1, Q = 0.5, ne = 200, "one", "nonsmooth")
> fit3 <- qrismb(fm, data = dat, t0 = 1, Q = 0.5, ne = 200, "random", "iterative")
> 
> coef(fit1)
(Intercept)           X 
  1.2395248   0.8525343 
> summary(fit2)
Call:
qrismb(formula = fm, data = dat, t0 = 1, Q = 0.5, ne = 200, init = "one", 
    method = "nonsmooth")

qrismb Estimator
            estimate std.Error z.value   p.value    
(Intercept)   1.2528    0.1002  12.503 < 2.2e-16 ***
X             0.8100    0.1339   6.049 < 2.2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
> vcov(fit3)
            (Intercept)          X
(Intercept)   0.1255075 -0.1245413
X            -0.1245413  0.1581005
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
> fit1 <- qrismb(fm, data = lung2, t0 = 0, Q = 0.5, ne = 200, "random", "smooth")
> fit2 <- qrismb(fm, data = lung2, t0 = 30, Q = 0.5, ne = 200, "one", "nonsmooth")
> fit3 <- qrismb(fm, data = lung2, t0 = 100, Q = 0.5, ne = 200, "rq", "iterative")
> 
> coef(fit1)
(Intercept)         age         sex 
 9.41907305 -0.06826185  1.69483630 
> summary(fit2)
Call:
qrismb(formula = fm, data = lung2, t0 = 30, Q = 0.5, ne = 200, 
    init = "one", method = "nonsmooth")

qrismb Estimator
            estimate std.Error z.value p.value    
(Intercept)   5.6362    0.8095   6.962  <2e-16 ***
age          -0.0015    0.0120  -0.122  0.9033    
sex           0.4489    0.1819   2.468  0.0136 *  
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
> vcov(fit3)
             (Intercept)           age          sex
(Intercept)  0.655922283 -0.0092565369 -0.071841802
age         -0.009256537  0.0001320656  0.000887068
sex         -0.071841802  0.0008870680  0.031495760
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
