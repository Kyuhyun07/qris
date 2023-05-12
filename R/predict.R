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
        X <- model.matrix(object$formula, dat = object$data)
    else {
        object$formula[[2]] <- NULL
        X <- model.matrix(object$formula, dat = newdata)
    }
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
    if (missing(newdata)) {
        X <- model.matrix(object$formula, dat = object$data)
        return(exp(drop(X %*% object$coef)) + object$para$t0 - object$data[,1])
    }
    newY <- drop(newdata[all.vars(object$formula[[2]])[1]])
    object$formula[[2]] <- NULL
    newX <- model.matrix(object$formula, dat = newdata)
    exp(drop(newX %*% object$coef)) + object$para$t0 - newY        
}
