#' Print method for class NCM
#'
#' Prints statistics for a class NCM object.
#'
#' @param x Output of the function neutral_model.
#' @param ... Other parameters passed to the generic print function.
#'
#' @return Printed summary.
#' @export
#'
#' @author John Quensen
#' @examples
#' data(rslt)
#' rslt
#'
print.NCM <- function(x, ...) {
  cat(paste("Class: ",class(x),"\n"))
  cat("Ntm = ", x$stats$Ntm,"\n")
  cat("Mean Nt = ", x$stats$mean.Nt,"\n")
  cat("m = ", x$stats$m,"\n")
}
