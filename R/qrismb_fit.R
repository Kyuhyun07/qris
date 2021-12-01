#' A wrapper for the different qrismb fit
#'
#' This version is more flexible for future extension
#' 
#' @noRd
qrismb.fit <- function(info, method) {
  switch(method,
         smooth = qrismb.smooth(info),
         iterative = qrismb.iter(info),
         nonsmooth = qrismb.nonsmooth(info))
}

qrismb.nonsmooth <- function(info) {
  ## 1. : L1-minimization : Estimating equation for estimaing beta (using rq)
  out <- with(info, {
    M <- 1e6
    pseudo1 <- -colSums(X * I * W)
    pseudo2 <- 2 * colSums(X * I * Q)
    Y.reg <- c(data[,2], M, M)
    X.reg <- rbind(X, rbind(pseudo1, pseudo2))
    W.reg <- c(I * W, 1, 1)
    Li.fit <- rq.wfit(X.reg, Y.reg, weights = W.reg)
    coefficient <- as.vector(Li.fit$coefficients)
    ## Full multiplier bootstrap
    fb_result <- c()
    ## if (Li.fit$code == 1 | Li.fit$code == 2){
    if (all(Li.fit$coefficients <= 10)){
      for (j in 1:ne){
        ## generating perturbation variable
        eta <- rexp(n, 1)
        if (all(data[, 4] == rep(1, n))){
          W_star <- rep(1, n)
        } else {
          Gest <- ghat(data[, 1], 1 - data[, 4], eta)
          ghatstart0 <- 1
          if (t0 > Gest$deathtime[1]) ghatstart0 <- Gest$survp[min(which(Gest$deathtime>t0)) - 1]
          W_star <- data[,4] / Gest$survp[findInterval(data[,1] , Gest$deathtime)] * ghatstart0
          W_star[is.na(W_star)] <- max(W_star, na.rm = TRUE)
        }
        rev_pseudo1 <- -colSums(X * I * W_star * eta)
        rev_pseudo2 <- 2 * colSums(X * I * eta * Q)
        rev_Y.reg <- c(data[,2], M, M)
        rev_X.reg <- rbind(X,rbind(rev_pseudo1, rev_pseudo2))
        rev_W.reg <- c(W_star * I * eta, 1, 1)
        full_boot <- rq.wfit(rev_X.reg, rev_Y.reg, weights = rev_W.reg)
        if (all(full_boot$coef <= 10)){
          fb_result <- cbind(fb_result, as.vector(full_boot$coef))
        } else {
          fb_result <- cbind(fb_result, rep(NA, nc))
        }
      }
      sigma <- try(cov(t(fb_result), use = "complete.obs"), silent = T)
      if(class(sigma)[1] == "try-error" | sum(!is.na(colSums(fb_result))) <= 1){
        stop("More resampling iterations are necessary")
      } else {
        se <- sqrt(diag(sigma))
        out <- list(coefficient = coefficient, stderr = se, vcov = sigma)
      }
    } else {
      coefficient <- se <- rep(NA, nc)
      vcov <- matrix(NA, nc, nc)
      out <- list(coefficient = coefficient, stderr = se, vcov = vcov)
    }})
  out
}

qrismb.iter <- function(info) {
  out <- with(info, {
    iter_beta_result <- old_beta <- new_beta <- betastart
    new_sigma <- old_sigma <- diag(1, nc, nc)
    new_h <- old_h <- H
    iter_SE_result <- sqrt(diag(new_h))
    iter_norm_result <- c()
    for (k in 1:10){
      old_beta <- new_beta
      old_sigma <- new_sigma
      old_h <- new_h
      slope_olda <- Amat(old_beta, X, W, old_h, I, logZ, Q)
      ## Step 1 : Update beta()
      new_beta <- old_beta - qr.solve(slope_olda) %*% (isObj(old_beta, X, W, old_sigma, I, logZ, Q) / n)
      iter_beta_result <- rbind(iter_beta_result, t(new_beta))
      ## Step 2 : Update Sigma()
      result.ismb <- c()
      for (j in 1:ne){
        ## generating perturbation variable
        eta <- rexp(n, 1)
        if (all(data[, 4] == rep(1, n))){
          W_star <- rep(1, n)
        } else {
          Gest <- ghat(data[, 1], 1 - data[, 4], eta)
          ghatstart0 <- 1
          if (t0 > Gest$deathtime[1]) ghatstart0 <- Gest$survp[min(which(Gest$deathtime>t0))-1]
          W_star <- data[,4] / Gest$survp[findInterval(data[,1], Gest$deathtime)] * ghatstart0
          W_star[is.na(W_star)] <- max(W_star, na.rm = TRUE)
        }
        result <- rev_isObj(new_beta, X, W_star, old_h, eta, I, logZ, Q)
        result.ismb <- cbind(result.ismb, result)
      }
      new_V <- cov(t(result.ismb))
      slope_newa <- Amat(new_beta, X, W, old_h, I, logZ, Q)
      new_sigma <- t(qr.solve(slope_newa)) %*% new_V %*% qr.solve(slope_newa)
      new_h <- new_sigma
      iter_SE_result <- rbind(iter_SE_result , sqrt(diag(new_sigma)))
      iter_norm_result <- rbind(iter_norm_result , norm(new_beta-old_beta, "F"))
      if(norm(new_beta-old_beta, "i") < 1e-4) break
    } ## end for loop
    ## Last iteration
    old_beta <- new_beta
    old_sigma <- new_sigma
    old_h <- new_h
    slope_olda <- Amat(old_beta, X, W, old_h, I, logZ, Q)
    new_beta <- old_beta - qr.solve(slope_olda) %*% (isObj(old_beta, X, W, old_sigma, I, logZ, Q) / n)
    iter_beta_result <- rbind(iter_beta_result, t(new_beta))
    result.ismb <- c()
    for (j in 1:ne){
      ## generating perturbation variable
      eta <- rexp(n, 1)
      if (all(data[, 4] == rep(1, n))){
        W_star <- rep(1,n)
      } else {
        Gest <- ghat(data[,1],1-data[,4],eta)
        ghatstart0 <- 1        
        if (t0 > Gest$deathtime[1]) ghatstart0 <- Gest$survp[min(which(Gest$deathtime > t0)) - 1]
        W_star <- data[,4] / Gest$survp[findInterval(data[,1], Gest$deathtime)] * ghatstart0
        W_star[is.na(W_star)] <- max(W_star, na.rm = TRUE)
      }
      result <- rev_isObj(new_beta, X, W_star, old_h, eta, I, logZ, Q)
      result.ismb <- cbind(result.ismb,result)
    }
    new_V <- cov(t(result.ismb))
    slope_newa <- Amat(new_beta, X, W, old_h, I, logZ, Q)
    new_sigma <- t(qr.solve(slope_newa)) %*% new_V %*% qr.solve(slope_newa)
    new_h <- new_sigma
    iter_SE_result <- rbind(iter_SE_result , sqrt(diag(new_sigma)))
    iter_norm_result <- rbind(iter_norm_result , norm(new_beta - old_beta, "F"))
    out <- list(coefficient = tail(iter_beta_result, n = 1),
                coefficient_result = iter_beta_result,
                stderr = tail(iter_SE_result, n = 1),
                stderr_result = iter_SE_result,
                vcov = new_sigma, iterno = k + 1,
                norm = iter_norm_result)
  })
  out
}

qrismb.smooth <- function(info) {
  out <- with(info, {
    rcpp.fit <- nleqslv(betastart, function(b) isObj(b, X, W, H, I, logZ, Q))
    if (rcpp.fit$termcd == 1 | rcpp.fit$termcd == 2) {
      coefficient <- rcpp.fit$x
      ## Variance estimation : ISMB
      rcpp.result.ismb <- c()
      for (j in 1:ne){
        ## generating perturbation variable
        eta <- rexp(n,1)
        if (all(data[, 4] == rep(1, n))){
          W_star <- rep(1, n)
        } else {
          Gest <- ghat(data[,1],1-data[,4],eta)
          ghatstart0 <- 1
          if (t0 > Gest$deathtime[1]) ghatstart0 <- Gest$survp[min(which(Gest$deathtime>t0))-1]
          W_star <- data[,4] /
            Gest$survp[findInterval(data[,1] , Gest$deathtime)] * ghatstart0
          W_star[is.na(W_star)] <- max(W_star, na.rm = TRUE)
        }
        result <- rev_isObj(coefficient, X, W_star, H, eta, I, logZ, Q) / n
        rcpp.result.ismb <- cbind(rcpp.result.ismb,result)
      }
      v <- cov(t(rcpp.result.ismb))
      rcpp.a <- Amat(coefficient, X, W, H, I, logZ, Q) / n
      inva <- try(solve(rcpp.a))
      if(class(inva)[1] == "try-error") {
        se <- rep(NA,nc+1)
        stop("Slope matrix is singular matrix")
      } else {
        sigma <- t(inva) %*% v %*% inva
        se <- sqrt(diag(sigma))
      }
      out <- list(coefficient=coefficient, stderr = se, vcov = sigma)
    } else {
      coefficient <- se <- rep(NA, nc)
      vcov <- matrix(NA, nc, nc)
      out <- list(coefficient=coefficient, stderr = se, vcov = vcov)
    }
  })
  out
}
