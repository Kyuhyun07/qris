#' @exportS3Method coef qrismb
coef.qrismb <- function(x, ...) {
    coef <- as.numeric(x$coefficient)
    names(coef) <- x$varNames
    coef
}

#' @exportS3Method vcov qrismb
vcov.qrismb <- function(x, ...) {
    vcov <- x$vcov
    colnames(vcov) <- rownames(vcov) <- x$varNames
    vcov
}

#' @exportS3Method print qrismb
print.qrismb <- function(x, ...) {
    cat("Call: \n")
    dput(x$call)
    mat <- rbind(x$varNames, format(x$coefficient, digits = 5))
    prmatrix(mat, rowlab = rep("", nrow(mat)), collab = rep("", ncol(mat)), quote = FALSE)
}

#' @exportS3Method summary qrismb
summary.qrismb <- function(x, ...) {
    if (class(x) != "qrismb"){
        stop("Must be qrismb class")
    }
    ans <- x["call"]
    est.qrismb <- x$coefficient
    if (is.null(x$stderr)) {
        se.qrismb <- rep(NaN, length(est.qrismb))
    } else {
        se.qrismb <- x$stderr
    }
    z.qrismb <- as.numeric(est.qrismb)/as.numeric(se.qrismb)
    TAB <- cbind(estimate = round(est.qrismb, 4),
                 std.Error = round(se.qrismb, 4),
                 z.value = round(z.qrismb, 3),
                 p.value = round(2 * pnorm(-abs(z.qrismb)), 4))
    rownames(TAB) <- x$varNames
    out <- list(call = x$call, coefficients=TAB)
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
