test_that("Test Smooth PE", {
    data(cancer, package = "survival")
    lung2 <- subset(lung, select = c(time, status, age, sex))
    ## tidy up the data
    lung2$status <- lung2$status - 1
    lung2$sex <- lung2$sex - 1    
    fm <- Surv(time, status) ~ age + sex
    set.seed(1)
    fit1 <- qris(fm, data = lung2, t0 = 100, Q = 0.5, nB = 0, "smooth", "pmb", "rq")
    expect_equal(round(fit1$coef, 3), c(8.628, -0.060, 1.812), tolerance = .01)
    set.seed(1)
    fit2 <- qris(fm, data = lung2, t0 = 100, Q = 0.5, nB = 100, "iterative", "pmb", "rq")
    expect_equal(round(fit2$coef, 3), c(5.978, -0.010, 0.546), tolerance = .01)
    set.seed(1)
    fit3 <- qris(fm, data = lung2, t0 = 100, Q = 0.5, nB = 0, "nonsmooth", "pmb", "rq")
    expect_equal(round(fit3$coef, 3), c(5.995, -0.010, 0.540), tolerance = .01)
})

test_that("Test pmb", {
    data(cancer, package = "survival")
    lung2 <- subset(lung, select = c(time, status, age, sex))
    ## tidy up the data
    lung2$status <- lung2$status - 1
    lung2$sex <- lung2$sex - 1    
    fm <- Surv(time, status) ~ age + sex
    set.seed(1)
    fit1 <- qris(fm, data = lung2, t0 = 100, Q = 0.5, nB = 100, "smooth", "pmb", "rq")
    expect_equal(round(fit1$stderr, 3), c(2.404, 0.039, 0.637), tolerance = .01)
    set.seed(1)
    fit2 <- qris(fm, data = lung2, t0 = 100, Q = 0.5, nB = 100, "iterative", "pmb", "rq")
    expect_equal(round(fit2$stderr, 3), c(0.464, 0.007, 0.133), tolerance = .01)
})

test_that("Test fmb", {
    data(cancer, package = "survival")
    lung2 <- subset(lung, select = c(time, status, age, sex))
    ## tidy up the data
    lung2$status <- lung2$status - 1
    lung2$sex <- lung2$sex - 1    
    fm <- Surv(time, status) ~ age + sex
    set.seed(1)
    fit1 <- qris(fm, data = lung2, t0 = 100, Q = 0.5, nB = 100, "smooth", "fmb", "rq")
    expect_equal(round(fit1$stderr, 3), c(2.751, 0.044, 0.818), tolerance = .01)
    set.seed(1)
    fit2 <- qris(fm, data = lung2, t0 = 100, Q = 0.5, nB = 100, "iterative", "fmb", "rq")
    expect_equal(round(fit2$stderr, 3), c(0.737, 0.011, 0.192), tolerance = .01)
    set.seed(1)
    fit3 <- qris(fm, data = lung2, t0 = 100, Q = 0.5, nB = 100, "nonsmooth", "fmb", "rq")
    expect_equal(round(fit3$stderr, 3), c(1.203, 0.019, 0.263), tolerance = .01)
})

