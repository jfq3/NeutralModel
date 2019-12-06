#' pval_NeutralModel
#'
#' Makes a test of signifcance for the neutral model.
#'
#' @usage pval_NeutralModel(nm.rslt, n.iter=99999, plot=TRUE)
#' @param nm.rslt Result from neutral_model function.
#' @param n.iter Number of iterations.
#' @param plot A logical; TRUE plots simulated distribution, with red line denoting observed statistic.
#'
#' @return A list of two items: the observed RMS error and associated p.values.
#' @export
#'
#' @details This test is not really valid. Non-linear models can only be tested against each other
#'   to see which is better. All are signficant by this method.
#' @author Arvind Venkataraman and John Quensen
#' @examples
#' data(rslt)
#' pval_NeutralModel(rslt)
#'
pval_NeutralModel <- function(nm.rslt, n.iter=99999, plot=TRUE){
  obs.rmse <- sqrt(mean(nm.rslt$bestneutralmatrix[ , 5]))
  obs.freq <- nm.rslt$obs[ , 1]
  pred.freq <- nm.rslt$bestneutralmatrix[ , 4]

  rand.rmse.fun <- function(obs.freq, pred.freq) {
    sim.freq <- sample(pred.freq)
    sim.rmse <- sqrt(mean((obs.freq-sim.freq)^2))
    sim.rmse
  }

  rmse.dist <- replicate(n.iter, rand.rmse.fun(obs.freq, pred.freq))

  if (plot) {
    m <- 0.75*(min(obs.rmse, min(rmse.dist)))
    n <- 1.25*(max(obs.rmse, max(rmse.dist)))
    graphics::hist(rmse.dist, xlim=c(m,n), xlab="Root Mean Square Error",
         main="Simulated RMSE Distribution")
    graphics::abline(v=obs.rmse, lwd=2, col="red")
    graphics::box(which="plot")
  }

  p.val <- (rank(c(obs.rmse, rmse.dist))[1])/(n.iter+1)
  rslt <- list(obs.rmse=obs.rmse, p.val=p.val)
  class(rslt) <- "NCMpval"
  return(rslt)
}
