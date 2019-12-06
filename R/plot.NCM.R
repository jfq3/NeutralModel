#' Plot method for class NCM
#'
#' Plots an NCM object.
#'
#' @method plot NCM
#' @usage plot(x, xlab = "Relative Abundance in Source",
#'   ylab = "Det.Freq. in Target", main = NULL,
#'   xlim = c(0.001, 100), ylim = c(0, 100), ...,
#'   symb.col = c("darkgrey", "indianred", "darkgreen"))
#' @param x Output of the function neutral_model.
#' @param xlab Y-axis label.
#' @param ylab X-axis label.
#' @param main Plot title.
#' @param xlim X-axis limits.
#' @param ylim Y-axis limits
#' @param ... Other parmeters passed to the generic plot function.
#' @param symb.col Vector of symbol colors for neutral, selected against, and selected for OTUs.
#'
#' @return A plot.
#' @export
#' @author John Quensen
#' @examples
#' data(rslt)
#' plot(rslt, xlab="Relative Abundance in Mouth", ylab="Det. Freq in Lung")
#'
plot.NCM <- function(x,
                     xlab="Relative Abundance in Source",
                     ylab="Det.Freq. in Target",
                     main=NULL,
                     xlim=c(0.001, 100),
                     ylim=c(0,100), ...,
                     symb.col=c("darkgrey", "indianred", "darkgreen")) {
  graphics::plot(x$obs[,4]*100,x$bestneutralmatrix[,4]*100,type="n",col="black",lwd=2,
       xlab=xlab, ylab=ylab,log="x", frame.plot=T, xlim=xlim, ylim=ylim, main=main)
  graphics::points(jitter(x$neutral_matrix[,4]*100),x$neutral_matrix[,1]*100,col=symb.col[1],pch=16)
  if(nrow(x$against_matrix) > 0) {graphics::points(jitter(x$against_matrix[,4]*100),
                                           jitter(x$against_matrix[,1]*100),
                                           col=symb.col[2],pch=16,cex=1.25)}
  if(nrow(x$for_matrix) > 0) {graphics::points(jitter(x$for_matrix[,4]*100),
                                       jitter(x$for_matrix[,1]*100),col=symb.col[3],
                                       pch=16,cex=1.25)}
  graphics::lines(x$obs[,4]*100,x$bestneutralmatrix[,6]*100,type="l",col="black",lwd=2,lty=2)
  graphics::lines(x$obs[,4]*100,x$bestneutralmatrix[,7]*100,type="l",col="black",lwd=2,lty=2)
  graphics::lines(x$obs[,4]*100,x$bestneutralmatrix[,4]*100,type="l",col="black",lwd=2)
}
