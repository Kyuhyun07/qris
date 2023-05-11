#' Prediction for Quantile Regression Model Fitted on Residual life
#' 
#' Prediction based on fitted quantile regression model
#'
#' @param object is a qris object
#' @param newdata is a data frame for an optional new data to do predictions.
#' If omitted, the fitted values based on the original data and fit will be returned.
#' @param ... for future extension
#'
#' @method predict qris
#'
#' @export
#' @return A vector of prediction
predict.qris <- function(object, newdata, ...) {
    if (missing(newdata)) 
        X <- model.matrix(formula(object$call[[2]]), dat = object$data)
    else
        X <- model.matrix(formula(object$call[[2]]), dat = newdata)
    exp(drop(X %*% object$coef)) + object$para$t0
}


#' Residuals for Quantile Regression Model Fitted on Residual life
#'
#' Residual based on fitted quantile regression model
#' 
#' @method residuals qris
#'
#' @export
#' @return A vector of residual
residuals.qris <- function(object, newdata, ...) {
    X <- model.matrix(formula(object$call[[2]]), dat = object$data)
    exp(drop(X %*% object$coef)) + object$para$t0 - object$data[,1]    
}
