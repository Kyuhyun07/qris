#' @exportS3Method coef qris
#' @importFrom stats model.matrix na.omit printCoefmat
coef.qris <- function(object, ...) {
  coef <- as.numeric(object$coefficient)
  names(coef) <- object$varNames
  coef
}

#' @exportS3Method vcov qris
vcov.qris <- function(object, ...) {
  vcov <- object$vcov
  colnames(vcov) <- rownames(vcov) <- object$varNames
  vcov
}

#' @exportS3Method print qris
print.qris <- function(x, ...) {
  cat("Call: \n")
  dput(x$call)
  mat <- rbind(x$varNames, format(x$coefficient, digits = 5))
  prmatrix(mat, rowlab = rep("", nrow(mat)),
           collab = rep("", ncol(mat)), quote = FALSE)
}

#' @exportS3Method summary qris
summary.qris <- function(object, ...) {
  if (class(object) != "qris"){
    stop("Must be qris class")
  }
  ans <- object["call"]
  est.qris <- object$coefficient
  if (is.null(object$stderr)) se.qris <- rep(NaN, length(est.qris))
  else se.qris <- object$stderr
  z.qris <- as.numeric(est.qris) / as.numeric(se.qris)
  TAB <- data.frame(estimate = round(drop(est.qris), 4),
                    std.Error = round(drop(se.qris), 4),
                    z.value = round(z.qris, 3),
                    p.value = round(2 * pnorm(-abs(z.qris)), 4))
  rownames(TAB) <- object$varNames
  out <- list(call = object$call, coefficients = TAB)
  class(out) <- "summary.qris"
  out
}


#' @exportS3Method print summary.qris
print.summary.qris <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("qris Estimator")
  cat("\n")
  printCoefmat(as.matrix(x$coefficients), P.values = TRUE, has.Pvalue = TRUE)
}

#' @exportS3Method confint qris
#' @importFrom stats qnorm
confint.qris <- function(object, level = 0.95, ...) {
  cf <- coef(object)
  pnames <- names(cf)
  p <- (1 - level)/2
  p <- c(p, 1 - p)
  prange <- qnorm(p)
  pct <- paste(format(100 * p, trim = TRUE, scientific = FALSE, digits = 3),"%")
  ci <- array(NA, dim = c(length(pnames), 2L),
              dimnames = list(pnames, pct))
  ses <- object$stderr
  ci[] <- cf + ses %o% prange
  ci
}

