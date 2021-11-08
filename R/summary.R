#' @exportS3Method coef qrismb
coef.qrismb <- function(x, ...) {
    as.numeric(x$coefficient)
}

#' @exportS3Method print qrismb
print.qrismb <- function(x, ...) {
    cat("Call: \n")
    dput(x$call)
    mat <- rbind(x$varNames, format(x$coefficient, digits = 5))
    prmatrix(mat, rowlab = rep("", nrow(mat)), collab = rep("", ncol(mat)), quote = FALSE)
}
