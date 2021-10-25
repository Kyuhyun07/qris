#' @exportS3Method coef qrismb
coef.qrismb <- function(object, ...) {
    as.numeric(object$coefficient)
}
