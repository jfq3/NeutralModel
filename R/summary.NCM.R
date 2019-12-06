#' Summary method for class NCM
#'
#' Prints a summary for an NCM class object.
#'
#'
#' @param object Output of the function neutral_model.
#' @param ... Other parameters passed to the generic print function.
#'
#' @return A printed summary of the output of the function neutral_model.
#' @export
#' @author John Quensen
#' @examples
#' data(rslt)
#' summary(rslt)
summary.NCM <- function(object, ...) {
  print(object$sum1)
  cat("\n\n")
  print(object$sum2)
}
