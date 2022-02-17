#' @exportS3Method coef qrismb
#' @importFrom stats model.matrix na.omit printCoefmat
coef.qrismb <- function(object, ...) {
  coef <- as.numeric(object$coefficient)
  names(coef) <- object$varNames
  coef
}

#' @exportS3Method vcov qrismb
vcov.qrismb <- function(object, ...) {
  vcov <- object$vcov
  colnames(vcov) <- rownames(vcov) <- object$varNames
  vcov
}

#' @exportS3Method print qrismb
print.qrismb <- function(x, ...) {
  cat("Call: \n")
  dput(x$call)
  mat <- rbind(x$varNames, format(x$coefficient, digits = 5))
  prmatrix(mat, rowlab = rep("", nrow(mat)),
           collab = rep("", ncol(mat)), quote = FALSE)
}

#' @exportS3Method summary qrismb
summary.qrismb <- function(object, ...) {
  if (class(object) != "qrismb"){
    stop("Must be qrismb class")
  }
  ans <- object["call"]
  est.qrismb <- object$coefficient
  if (is.null(object$stderr)) se.qrismb <- rep(NaN, length(est.qrismb))
  else se.qrismb <- object$stderr
  z.qrismb <- as.numeric(est.qrismb) / as.numeric(se.qrismb)
  TAB <- data.frame(estimate = round(drop(est.qrismb), 4),
                    std.Error = round(drop(se.qrismb), 4),
                    z.value = round(z.qrismb, 3),
                    p.value = round(2 * pnorm(-abs(z.qrismb)), 4))
  rownames(TAB) <- object$varNames
  out <- list(call = object$call, coefficients = TAB)
  class(out) <- "summary.qrismb"
  out
}


#' @exportS3Method print summary.qrismb
print.summary.qrismb <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("qrismb Estimator")
  cat("\n")
  printCoefmat(as.matrix(x$coefficients), P.values = TRUE, has.Pvalue = TRUE)
}
