## #################################################################
## Simulation used in the manuscript
## Summarized from Simulation\ data\ generator.R and ndata2.R
## #################################################################

n <- 10

#' External variables
#' @param n is sample size
#' @param cen is censoring rate
#' @param sce scenario 1 assumes no treatment effect
#'
#' Internal variables
#' @param b is true parameter matrix
#' for the simulation in the manuscript it is 4 by 2; 4 t0's 2 X's

datGen <- function(n, cen, sce) {
    X <- model.matrix(~ sample(0:1, n, TRUE))
    if (sce == 1) b <- matrix(c(1.61, 1.41, 1.22, 1.04, 0, 0, 0, 0), 4)
    else b <- matrix(c(1.61, 1.41, 1.22, 1.04, 0.69, 0.80, 0.91, 1.02), 4)
    
}
