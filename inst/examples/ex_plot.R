data(cancer, package = "survival")
lung2 <- subset(lung, select = c(time, status, age, sex))
## tidy up the data
lung2$status <- lung2$status - 1
lung2$sex <- lung2$sex - 1

fm <- Surv(time, status) ~ age + sex
fit1 <- qrismb(fm, data = lung2, t0 = 0, Q = 0.5, ne = 200, "iterative", "pmb", "rq")
fit2 <- qrismb(fm, data = lung2, t0 = 30, Q = 0.5, ne = 200, "nonsmooth", "fmb", c(1, 0, 1))

summary(fit1)
summary(fit2)

## Give a error message when method = 'iterative'
plot(fit1)

## Plot with default values; Qs <- 1:9 / 10 and t0s = fit2$para$t0 (in this case 30)
plot(fit2)

## Plot with default values and without 95% CI;
## much faster, especially with multiples values of Qs and t0s
plot(fit2, ne = 0)
## Other specifications
plot(fit2, t0s = 1:6 * 10, ne = 0)
plot(fit2, Qs = 3:6 / 10, ne = 0)
plot(fit2, Qs = 3:6 / 10, t0s = 1:6 * 10, ne = 0)
plot(fit2, Qs = 3:6 / 10, t0s = 1:6 * 10, ne = 0, byQs = TRUE)

plot(fit2, t0s = 1:5 * 10)
plot(fit2, Qs = 3:6 / 10)
plot(fit2, Qs = 3:6 / 10, t0s = 1:6 * 10)
plot(fit2, Qs = 3:6 / 10, t0s = 1:6 * 10, byQs = TRUE)

## Choose what variables to plot
plot(fit2, t0s = 1:6 * 10)
plot(fit2, t0s = 1:6 * 10, vari = c("sex", "age"))
plot(fit2, t0s = 1:6 * 10, vari = "sex")

plot(fit2, Qs = 3:6 / 10, t0s = 1:6 * 10)
plot(fit2, Qs = 3:6 / 10, t0s = 1:6 * 10, vari = c("sex", "age"))
plot(fit2, Qs = 3:6 / 10, t0s = 1:6 * 10, vari = "sex")

## Export data feature:
## this allows users to create their own ggplot if they don't like our graphical parameters
## this is also useful if one wants to change byQs or vari
fit2 <- plot(fit2, Qs = 3:6 / 10, t0s = 1:6 * 10, exportDat = TRUE)
str(fit2$ggdat)
plot(fit2, byQs = FALSE)
plot(fit2, byQs = TRUE)

plot(fit2, byQs = FALSE, vari = c("sex", "age"))
plot(fit2, byQs = TRUE, vari = c("sex", "age"))
plot(fit2, byQs = FALSE, vari = "sex")
plot(fit2, byQs = TRUE, vari = "sex")

plot(fit2, byQs = FALSE, vari = fit2$varNames[1:2])
plot(fit2, byQs = TRUE, vari = fit2$varNames[1:2])
