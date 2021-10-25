data("retinopathy", package = 'survival')
reti <- retinopathy
reti_rev <- reti
reti_rev$risk <- scales::rescale(reti$risk, to = c(0,1), from = c(0,12))
reti_rev$type <- as.numeric(reti$type) - 1
reti_rev$typextrt <- reti_rev$type * reti$trt
Z <- reti_rev$futime
nc <- 5
covariate <- as.matrix(reti_rev[,c(9,4,5,6,10)])
D <- reti_rev$status
t_0 <- 1
Q <- 0.25
ne <- 200

fit1 <- qrismb(Z, covariate, D, t_0, Q, ne, "rq", "smooth")
fit2 <- qrismb(Z, covariate, D, t_0, Q, ne, "random", "nonsmooth")

coef(fit1)
coef(fit2)
