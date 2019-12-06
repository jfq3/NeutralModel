#' Plot method for class NCMpval
#'
#' Prints a summary for a class NCMpval object.
#'
#' @param x Output of the function pval_NeutralModel.
#' @param ... Other parameters passed to the generic print function.
#'
#' @return Printed summary.
#' @export
#' @author John Quensen
#'
print.NCMpval <- function(x, ...) {
  cat(paste("Class: ", class(x), "\n"))
  cat("Fit = ", x$Fit, "\n")
  cat("Fit std error = ", x$Fit.std.error, "\n")
  # cat("Observed Root Mean Square error = ", x$obs.rmse, "\n")
  cat("p-value = ", x$p.val,"\n")
}
