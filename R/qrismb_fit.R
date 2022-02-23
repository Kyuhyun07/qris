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
  ## 1. : L1-minimization : Estimating equation for estimating beta (using rq)
  out <- with(info, {
    M <- 1e6
    pseudo1 <- -colSums(X*I*W)
    pseudo2 <- 2 * colSums(X*I*Q)
    Y.reg <- c(data[,2], M, M)
    X.reg <- rbind(X, rbind(pseudo1, pseudo2))
    W.reg <- c(I * W, 1, 1)
    Li.fit <- rq.wfit(X.reg, Y.reg, weights = W.reg)
    coefficient <- as.vector(Li.fit$coefficients)
    if (all(Li.fit$coefficients <= 10)){
      if (ne <= 1) {
        se <- rep(NA, nc)
        vcov <- matrix(NA, nc, nc)
        return(list(coefficient = coefficient, stderr = se, vcov = vcov))
      }
      if (se == "fmb"){
        fmb.result <- c()
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
          rev_pseudo1 <- -colSums(X * W_star * eta)
          rev_pseudo2 <- 2 * colSums(X * I * eta * Q)
          rev_Y.reg <- c(data[,2], M, M)
          rev_X.reg <- rbind(X,rbind(rev_pseudo1, rev_pseudo2))
          rev_W.reg <- c(W_star * I * eta, 1, 1)
          fmb.fit <- rq.wfit(rev_X.reg, rev_Y.reg, weights = rev_W.reg)
          if (all(fmb.fit$coef <= 10)){
            fmb.result <- cbind(fmb.result, as.vector(fmb.fit$coef))
          } else {
            fmb.result <- cbind(fmb.result, rep(NA, nc))
          }
        }
        fmb.sigma <- try(cov(t(fmb.result), use = "complete.obs"), silent = T)
        if(class(fmb.sigma)[1] == "try-error" | sum(!is.na(colSums(fmb.result))) <= 1){
          stop("More resampling iterations are necessary")
        } else {
          fmb.se <- sqrt(diag(fmb.sigma))
          out <- list(coefficient = coefficient, stderr = fmb.se, vcov = fmb.sigma)
        }
      } else {
        stop("Only full multiplier bootstrapping (fmb) is available for nonsmooth estimating equation approach")
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
    new_h <- old_h <- H
    new_sigma <- old_sigma <- H*n
    iter_SE_result <- sqrt(diag(new_sigma))
    iter_norm_result <- c()
    if (se == "fmb") {
      for (k in 1:control$maxiter){
        old_beta <- new_beta
        old_sigma <- new_sigma
        old_h <- new_h
        slope_a <- Amat(old_beta, X, W, old_h, I, logZ, Q)/n
        ## Step 1 : Update beta()
        ## Singular matrix 'a' error message and break
        if (class(try(qr.solve(slope_a),silent=TRUE))[1]=="try-error") {
          warning("‘A’ matrix is singular during iteration. Please try the non-iterative method.")
          break
        } else {
          new_beta <- old_beta + qr.solve(slope_a) %*% (isObj(old_beta, X, W, old_h, I, logZ, Q)/n)
          iter_beta_result <- rbind(iter_beta_result, t(new_beta))
          ## Step 2 : Update Sigma()
          result.fmb <- c()
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
            fmb.fit <- nleqslv(old_beta, function(b) rev_isObj(b, X, W_star, old_h, eta, I, logZ, Q)/n)
            if (fmb.fit$termcd == 1 | fmb.fit$termcd == 2) {
              result.fmb <- cbind(result.fmb,fmb.fit$x)
            } else {
              result.fmb <- cbind(result.fmb,rep(NA, length(fmb.fit$x)))
            }
          }
          new_sigma <- try(cov(t(result.fmb), use = "complete.obs"), silent = T)
          ## Trace the result
          if (control$trace) {
            cat("\n Step:", k)
            cat("\n beta:", as.numeric(new_beta))
            cat("\n se:", as.numeric(sqrt(diag(new_sigma))), "\n")
          }
          if (class(try(new_sigma,silent=TRUE))[1]=="try-error") {
            warning("Futher fmb method is inapplicable to this dataset. Please try other estimation methods")
            new_sigma <- old_sigma
            break
          } else {
            new_h <- new_sigma
            iter_SE_result <- rbind(iter_SE_result , sqrt(diag(new_sigma)))
            iter_norm_result <- rbind(iter_norm_result , norm(new_beta-old_beta, "F"))
            if(iter_norm_result[k]>=100*iter_norm_result[1]) {
              warning("Point estimation failed to converge.")
              break
            }
            if(norm(new_beta-old_beta, "i") < control$tol) break
          }
        }
      } ## end for loop
      
      out <- list(coefficient = tail(iter_beta_result, n = 1),
                  trace.coefficient = iter_beta_result,
                  stderr = tail(iter_SE_result, n = 1),
                  trace.stderr = iter_SE_result,
                  vcov = new_sigma, iterno = k,
                  norm = iter_norm_result)
    } else {
      for (k in 1:control$maxiter){
        old_beta <- new_beta
        old_sigma <- new_sigma
        old_h <- new_h
        slope_a <- Amat(old_beta, X, W, old_h, I, logZ, Q)/n
        ## Step 1 : Update beta()
        ## Singular matrix 'a' error message and break
        if (class(try(qr.solve(slope_a),silent=TRUE))[1]=="try-error") {
          warning("‘A’ matrix is singular during iteration. Please try the non-iterative method.")
          break
        } else {
          new_beta <- old_beta + qr.solve(slope_a) %*% (isObj(old_beta, X, W, old_h, I, logZ, Q)/n)
          iter_beta_result <- rbind(iter_beta_result, t(new_beta))
          ## Step 2 : Update Sigma()
          result.pmb <- c()
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
            pmb.eval <- rev_isObj(old_beta, X, W_star, old_h, eta, I, logZ, Q)/n
            result.pmb <- cbind(result.pmb, pmb.eval)
          }
          new_V <- cov(t(result.pmb), use = "complete.obs")
          new_sigma <- t(qr.solve(slope_a)) %*% new_V %*% qr.solve(slope_a)
          new_h <- new_sigma
          ## Trace the result
          if (control$trace) {
            cat("\n beta:", as.numeric(new_beta), "\n")
            cat("\n se:", as.numeric(sqrt(diag(new_sigma))), "\n")
          }
          iter_SE_result <- rbind(iter_SE_result , sqrt(diag(new_sigma)))
          iter_norm_result <- rbind(iter_norm_result , norm(new_beta-old_beta, "F"))
          if(iter_norm_result[k]>=100*iter_norm_result[1]) {
            warning("Point estimation result is diverging")
            break
            }
          if(norm(new_beta-old_beta, "i") < control$tol) break
        }
      } ## end for loop
      
      # ## Last iteration
      # old_beta <- new_beta
      # old_sigma <- new_sigma
      # old_h <- new_h
      # slope_a <- Amat(old_beta, X, W, old_h, I, logZ, Q)/n
      # ## Step 1 : Update beta()
      # new_beta <- old_beta + qr.solve(slope_a) %*% (isObj(old_beta, X, W, old_h, I, logZ, Q)/n)
      # iter_beta_result <- rbind(iter_beta_result, t(new_beta))
      # ## Step 2 : Update Sigma()
      # result.pmb <- c()
      # for (j in 1:ne){
      #   ## generating perturbation variable
      #   eta <- rexp(n, 1)
      #   if (all(data[, 4] == rep(1, n))){
      #     W_star <- rep(1, n)
      #   } else {
      #     Gest <- ghat(data[, 1], 1 - data[, 4], eta)
      #     ghatstart0 <- 1
      #     if (t0 > Gest$deathtime[1]) ghatstart0 <- Gest$survp[min(which(Gest$deathtime>t0))-1]
      #     W_star <- data[,4] / Gest$survp[findInterval(data[,1], Gest$deathtime)] * ghatstart0
      #     W_star[is.na(W_star)] <- max(W_star, na.rm = TRUE)
      #   }
      #   pmb.eval <- rev_isObj(old_beta, X, W_star, old_h, eta, I, logZ, Q)/n
      #   result.pmb <- cbind(result.pmb, pmb.eval)
      # }
      # new_V <- cov(t(result.pmb), use = "complete.obs")
      # new_sigma <- t(qr.solve(slope_a)) %*% new_V %*% qr.solve(slope_a)
      # new_h <- new_sigma
      # iter_SE_result <- rbind(iter_SE_result , sqrt(diag(new_sigma)))
      # iter_norm_result <- rbind(iter_norm_result , norm(new_beta-old_beta, "F"))
      out <- list(coefficient = tail(iter_beta_result, n = 1),
                  trace.coefficient = iter_beta_result,
                  stderr = tail(iter_SE_result, n = 1),
                  trace.stderr = iter_SE_result,
                  vcov = new_sigma, iterno = k,
                  norm = iter_norm_result)
    }})
  out
}

qrismb.smooth <- function(info) {
  out <- with(info, {
    smooth.fit <- nleqslv(betastart, function(b) isObj(b, X, W, H, I, logZ, Q))
    if (smooth.fit$termcd == 1 | smooth.fit$termcd == 2) {
      coefficient <- smooth.fit$x
      if (ne <= 1) {
        se <- rep(NA, nc)
        vcov <- matrix(NA, nc, nc)
        return(list(coefficient = coefficient, stderr = se, vcov = vcov))
      }
      if (se == "fmb"){
        ## Full multiplier bootstrap
        smooth.fmb.result <- c()
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
          fmb.fit <- nleqslv(coefficient, function(b) rev_isObj(b, X, W_star, H, eta, I, logZ, Q))
          if (fmb.fit$termcd == 1 | fmb.fit$termcd == 2) {
            smooth.fmb.result <- cbind(smooth.fmb.result,fmb.fit$x)
          } else {
            smooth.fmb.result <- cbind(smooth.fmb.result,rep(NA, length(fmb.fit$x)))
          }
        }
        fmb.sigma <- try(cov(t(smooth.fmb.result), use = "complete.obs"), silent = T)
        if(class(fmb.sigma)[1] == "try-error" | sum(!is.na(colSums(smooth.fmb.result))) <= 1){
          stop("More resampling iterations are needed")
        } else {
          fmb.se <- sqrt(diag(fmb.sigma))
          out <- list(coefficient = coefficient, stderr = fmb.se, vcov = fmb.sigma)
        }
      } else if (se == "pmb") {
        ## Partial Multiplier Bootstrap
        smooth.pmb.result <- c()
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
          pmb.eval <- rev_isObj(coefficient, X, W_star, H, eta, I, logZ, Q) / n
          smooth.pmb.result <- cbind(smooth.pmb.result,pmb.eval)
        }
        pmb.v <- try(cov(t(smooth.pmb.result), use = "complete.obs"), silent = T)
        pmb.a <- Amat(coefficient, X, W, H, I, logZ, Q) / n
        ## Singular matrix 'a' error message and break
        if (class(try(qr.solve(pmb.a),silent=TRUE))[1]=="try-error") {
          pmb.se <- rep(NA,nc+1)
          stop("‘A’ matrix is singular during iteration. Please try the nonsmooth method." )
        } else {
          pmb.inva <- qr.solve(pmb.a)
          pmb.sigma <- t(pmb.inva) %*% pmb.v %*% pmb.inva
          pmb.se <- sqrt(diag(pmb.sigma))
        }
        out <- list(coefficient=coefficient, stderr = pmb.se, vcov = pmb.sigma)
      } else {
        print ("Select either 'fmb' or 'pmb'")
      }
    } else {
      coefficient <- se <- rep(NA, nc)
      vcov <- matrix(NA, nc, nc)
      out <- list(coefficient=coefficient, stderr = se, vcov = vcov)
    }})
  out
}
