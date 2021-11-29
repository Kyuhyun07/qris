## #########################################
## Simulated data
## #########################################
data.gen <- function(n) {
    r0 <- .2 * sqrt(log(2))
    r1 <- .1 * sqrt(log(2))
    dat <- data.frame(censoring = runif(n, 0, 24.35),
                      Time0 = sqrt(-log(1 - runif(n))),
                      X = rbinom(n, 1, .5))
    dat$Time0 <- ifelse(dat$X > 0, dat$Time0 / r1, dat$Time0 / r0)
    dat$Time <- pmin(dat$Time0, dat$censoring)
    dat$status <- 1 * (dat$Time0 < dat$censoring)
    subset(dat, select = c(Time, status, X))
}

dat <- data.gen(200)
fm <- Surv(Time, status) ~ X
fit1 <- qrismb(fm, data = dat, t0 = 1, Q = 0.5, ne = 200, "smooth", "rq")
fit2 <- qrismb(fm, data = dat, t0 = 1, Q = 0.5, ne = 200, "nonsmooth", "noeffect")
fit3 <- qrismb(fm, data = dat, t0 = 1, Q = 0.5, ne = 200, "iterative", c(2, 1))

summary(fit1)
summary(fit2)
summary(fit3)

## #########################################
## Real data application
## #########################################
data(cancer, package = "survival")
lung2 <- subset(lung, select = c(time, status, age, sex))
## tidy up the data
lung2$status <- lung2$status - 1
lung2$sex <- lung2$sex - 1

fm <- Surv(time, status) ~ age + sex
fit1 <- qrismb(fm, data = lung2, t0 = 0, Q = 0.5, ne = 200, "iterative")
fit2 <- qrismb(fm, data = lung2, t0 = 30, Q = 0.5, ne = 200, "nonsmooth", c(1, 0, 1))
fit3 <- qrismb(fm, data = lung2, t0 = 100, Q = 0.5, ne = 200,"smooth", "noeffect")

summary(fit1)
summary(fit2)
summary(fit3)
