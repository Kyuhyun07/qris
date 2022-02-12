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

#' Draw 95% confidence interval by a quantile regression estimator of residual lifetime from survival data
#'
#' @exportS3Method plot qrismb
#' @importFrom ggplot2 facet_wrap geom_line ggplot xlab ylab
#' @importFrom reshape melt
#' @importFrom ggpubr ggarrange
#'
#' @param object is an qrismb object
#' @param type is a indicator that which parameter designated to x value (1 = t0, 2 = quantile); when not specified, the default value is 2.
#' @param t0s is a vector of range of t0 to plot; when not specified, the default value is from 0 to presently defined t_0
#' @param Qs  is a vector of range of Q to plot; when not specified, the default value is from 5% to presently defined Q
#' @param ne is the number of multiplier bootstrapping
#' when not specified, multiplier bootstrap will be carried out with the \code{ne} specified in \code{object};
#' when ne = 0, only the point estimates will be plotted;
#' when ne > 1, both the point estimates and the 95% Wald CI will be plotted.
#### plot function ####
#### plot function ####
qrplot.qrismb <- function(object, type = 2, t0s = c(0,t0), Qs = c(0.05,Q), ne = 100, xlab, ylab, ...) {
  # Draw 95% CI plot of point estimate against t0 at different quantiles
  if (type == 1) {
    t0s <- seq(t0s[1], t0s[2], length.out=4)
    if (t0s[1] == t0s[4]){
      t0s <- t0s[1]
    }
    t0s <- round(t0s,2)
    Qs <- seq(Qs[1], Qs[2], by = 0.05)
    n <- length(Qs)
    result_coef <- result_se <- plots <- list()
    for (i in 1:n){
      q.now <- Qs[i]
      d <- as.data.frame(do.call(cbind, lapply(t0s, function(t) summary(update(fit1, t0 = t , Q = q.now, ne = 100))[[2]][,1:2])))
      d <- as.data.frame(t(d))
      d$t <- rep(t0s,each=2)
      # Seperate coef and StdErr data.frame (1=coef, 0=StdErr)
      d$row <- seq_len(nrow(d)) %% 2
      d_coef <- d[d$row == 1, ]
      d_se <- d[d$row == 0, ]
      d_coef$row <- NULL
      d_se$row <- NULL
      result_coef[[i]] <- melt(d_coef, id = "t")
      se <- melt(d_se, id = "t")$value
      result_coef[[i]] <- cbind(result_coef[[i]], se)
      plots[[i]] <- ggplot(result_coef[[i]], aes(x = t, y = value)) +
        geom_ribbon(aes(ymax = value + 1.96*se, ymin = value - 1.96*se),fill = "slategray3") +
        geom_line(color = "firebrick",size = 1) +
        facet_wrap(~ variable, scales = "free") +
        scale_x_continuous(limits=c(0, t0), breaks = t0s) +
        scale_y_continuous(limits=c(min(result_coef[[i]]$value - 4*result_coef[[i]]$se), max(result_coef[[i]]$value + 4*result_coef[[i]]$se))) +
        xlab(expression(t[0])) + ylab(expression(beta))
    }
    Qs <- round(Qs,2)
    Qs.label <- paste("Q = ", Qs, sep="")
    ggarrange(plotlist=plots, ncol = 1, nrow = n, labels = Qs.label, hjust = 0, vjust = 3, font.label = list(size = 9, color = "blue", face = "bold.italic"))
    # Draw 95% CI plot of point estimate against quantile at different t0
  } else if (type == 2) {
    t0s <- seq(t0s[1], t0s[2], length.out=4)
    if (t0s[1] == t0s[4]){
      t0s <- t0s[1]
    }
    Qs <- seq(Qs[1], Qs[2], by = 0.05)
    n <- length(t0s)
    result_coef <- result_se <- plots <- list()
    for (i in 1:n){
      t0.now <- t0s[i]
      d <- as.data.frame(do.call(cbind, lapply(Qs, function(q) summary(update(fit1, t0 = t0.now , Q = q, ne = 100))[[2]][,1:2])))
      d <- as.data.frame(t(d))
      d$Q <- rep(Qs,each=2)
      # Seperate coef and StdErr data.frame (1=coef, 0=StdErr)
      d$row <- seq_len(nrow(d)) %% 2
      d_coef <- d[d$row == 1, ]
      d_se <- d[d$row == 0, ]
      d_coef$row <- NULL
      d_se$row <- NULL
      result_coef[[i]] <- melt(d_coef, id = "Q")
      se <- melt(d_se, id = "Q")$value
      result_coef[[i]] <- cbind(result_coef[[i]], se)
      plots[[i]] <- ggplot(result_coef[[i]], aes(x = Q, y = value)) +
        geom_ribbon(aes(ymax = value + 1.96*se, ymin = value - 1.96*se),fill = "slategray3") +
        geom_line(color = "firebrick",size = 1) +
        facet_wrap(~ variable, scales = "free") +
        scale_x_continuous(limits=c(0.05, Q), breaks = Qs) +
        scale_y_continuous(limits=c(min(result_coef[[i]]$value - 4*result_coef[[i]]$se), max(result_coef[[i]]$value + 4*result_coef[[i]]$se))) +
        xlab(expression(tau)) + ylab(expression(beta))
    }
    t0s <- round(t0s,2)
    t0s.label <- paste("t0 = ", t0s, sep="")
    ggarrange(plotlist=plots, ncol = 1, nrow = n, labels = t0s.label, hjust = 0, vjust = 3, font.label = list(size = 9, color = "blue", face = "bold.italic"))
  } else {
    stop("Please choose x variable of graph (either 1 = t0 or 2 = quantile")
  }
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
