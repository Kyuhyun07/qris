## #########################################
## Real data application
## #########################################
library(qrismb)

data(cancer, package = "survival")
lung2 <- subset(lung, select = c(time, status, age, sex))
## tidy up the data
lung2$status <- lung2$status - 1
lung2$sex <- lung2$sex - 1

fm <- Surv(time, status) ~ age + sex
fit1 <- qrismb(fm, data = lung2, t0 = 0, Q = 0.5, ne = 200, "iterative")

library(ggplot2)
library(reshape)
## By Q (I wanted to set ne = 0 but couldn't. see #1 in to-do)
Qs <- 1:19 / 20
d <- as.data.frame(do.call(rbind, lapply(Qs, function(q) coef(update(fit1, Q = q, ne = 2)))))
d$Q <- Qs
d <- melt(d, id = "Q")

## effect plot without 95% CI
ggplot(d, aes(x = Q, y = value)) + geom_line() +
  facet_wrap(~ variable, scales = "free") +
  xlab(expression(tau)) + ylab(expression(beta))

## By t0
t0s <- 0:20 / 10
d <- as.data.frame(do.call(rbind, lapply(Qs, function(q) coef(update(fit1, t0 = q, ne = 2)))))
d$Q <- Qs
d <- melt(d, id = "Q")

## effect plot without 95% CI
ggplot(d, aes(x = Q, y = value)) + geom_line() +
  facet_wrap(~ variable, scales = "free") +
  xlab(expression(t[0])) + ylab(expression(beta))

## To-do
## 1. To save computing time, we might want to default \code{ne = 0} when making the plots.
##    To make this happen, we need to modify \code{qrismb()} so that bootstrap is not called when ne = 0.
## 2. If user wants to see the 95% wald CI in coefficient plots,
##    then allow them to specify an ne > 1 or use the ne from the original model fit.
## 3. When 95% CI is enabled, need to modify the ggplot function to alllow extra geom_line()
##    or geom_ribbon()
## 4. Allow user to plot by Q or by t0
## 5. I see the error "singular matrix 'a' in solve", but this randomly pops up. Why?
## 6. Might want to change the x and y labels or allow users to specify those.
## 7. I used free scales here. Should we adjust the y-scale so that beta = 0 is always visible? 
## 8. Might want the users to have an option to choose what variables to choose.
##    Can this be done at the ggplot level so we don't need to rerun the lapply()?
## 9. Need to clean the x and y ticks. Maybe we can show 5 ticks by default? 
## 10. Give some default Qs and t0s if the user is too lazy to choose those.
## 11. Implement those in a function, e.g., S3 method - plot() - for qrismb objects.
##     We could start with the following template
## 12. We will require at least one of t0 and Q to be a vector in the following for now.

#' @exportS3Method plot qrismb
#' @importFrom ggplot2 facet_wrap geom_line ggplot xlab ylab
#' @importFrom reshape melt
#'
#' @param object is an qrismb object
#' @param t0 is a vector of t0 to plot; when not specified, the default value is ...
#' @param Q  is a vector of Q to plot; when not specified, the default value is ...
#' @param ne is the number of multiplier bootstrapping;
#' when not specified, multiplier bootstrap will be carried out with the \code{ne} specified in \code{object};
#' when ne = 0, only the point estimates will be plotted;
#' when ne > 1, both the point estimates and the 95% Wald CI will be plotted.
plot.qrismb <- function(object, t0 = NULL, Q = NULL, ne = NULL, xlab, ylab, ...) {
  ...
}
  
