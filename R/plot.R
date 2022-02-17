#' Draw 95\% confidence interval by a quantile regression estimator of residual lifetime from survival data
#'
#' @param object is an qrismb object
#' @param type is a indicator that which parameter designated to x value (1 = t0, 2 = quantile); when not specified, the default value is 2.
#' @param t0s is a vector of range of t0 to plot; when not specified, the default value is from 0 to presently defined \eqn{t_0}
#' @param Qs  is a vector of range of Q to plot; when not specified, the default value is from 5\% to presently defined \eqn{Q}
#' @param ne is the number of multiplier bootstrapping
#' when not specified, multiplier bootstrap will be carried out with the \code{ne} specified in \code{object};
#' when \eqn{ne = 0}, only the point estimates will be plotted;
#' when \eqn{ne > 1}, both the point estimates and the 95\% Wald CI will be plotted.
#'
#' @importFrom ggplot2 facet_wrap geom_line ggplot xlab ylab aes geom_ribbon
#' @importFrom ggplot2 ggtitle scale_x_continuous scale_y_continuous theme element_text
#' @importFrom reshape melt
#' @importFrom ggpubr ggarrange
#' @export
#' @method plot qrismb
plot.qrismb <- function(object, type = 2, t0s = c(0,t0), Qs = c(0.05,Q), ne = 100, ...) {
  ## Draw 95% CI plot of point estimate against t0 at different quantiles
  if (type == 1) {
    t0s <- seq(t0s[1], t0s[2], length.out = 4)
    if (t0s[1] == t0s[4]) t0s <- t0s[1]
    t0s <- round(t0s,2)
    Qs <- seq(Qs[1], Qs[2], by = 0.05)
    n <- length(Qs)
    result_coef <- result_se <- plots <- list()
    for (i in 1:n){
      q.now <- Qs[i]
      d <- as.data.frame(do.call(cbind, lapply(t0s, function(t)
        summary(update(fit1, t0 = t , Q = q.now, ne = 100))[[2]][,1:2])))
      d <- as.data.frame(t(d))
      d$t <- rep(t0s,each=2)
      ## Seperate coef and StdErr data.frame (1=coef, 0=StdErr)
      d$row <- seq_len(nrow(d)) %% 2
      d_coef <- d[d$row == 1, ]
      d_se <- d[d$row == 0, ]
      d_coef$row <- NULL
      d_se$row <- NULL
      result_coef[[i]] <- melt(d_coef, id = "t")
      se <- melt(d_se, id = "t")$value
      result_coef[[i]] <- cbind(result_coef[[i]], se)
      plots[[i]] <- ggplot(result_coef[[i]], aes(x = t, y = value)) +
        geom_ribbon(aes(ymax = value + 1.96 * se, ymin = value - 1.96 * se),fill = "slategray3") +
        geom_line(color = "firebrick",size = 1) +
        facet_wrap(~ variable, scales = "free") +
        scale_x_continuous(limits = c(0, max(result_coef[[i]]$value)), breaks = t0s) +
        scale_y_continuous(limits = c(min(result_coef[[i]]$value - 4 * result_coef[[i]]$se),
                                      max(result_coef[[i]]$value + 4 * result_coef[[i]]$se))) +
        xlab(expression(t[0])) + ylab(expression(beta)) +
        ggtitle(paste0("Q = ", round(q.now, 2))) +
        theme(plot.title = element_text(hjust = 0.5))
    }
    ggarrange(plotlist=plots, ncol = 1, nrow = n)
    ## Draw 95% CI plot of point estimate against quantile at different t0
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
      d <- as.data.frame(do.call(cbind, lapply(Qs, function(q)
        summary(update(fit1, t0 = t0.now , Q = q, ne = 100))[[2]][,1:2])))
      d <- as.data.frame(t(d))
      d$Q <- rep(Qs,each=2)
      ## Seperate coef and StdErr data.frame (1=coef, 0=StdErr)
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
        scale_x_continuous(limits = c(0.05, max(result_coef[[i]]$Q)), breaks = Qs) +
        scale_y_continuous(limits = c(min(result_coef[[i]]$value - 4 * result_coef[[i]]$se),
                                      max(result_coef[[i]]$value + 4 * result_coef[[i]]$se))) +
        xlab(expression(tau)) + ylab(expression(beta)) +
        ggtitle(paste0("t0 = ", round(t0.now, 2))) +
        theme(plot.title = element_text(hjust = 0.5))
    }
    ggarrange(plotlist = plots, ncol = 1, nrow = n)
  } else {
    stop("Please choose x variable of graph (either 1 = t0 or 2 = quantile")
  }
}
