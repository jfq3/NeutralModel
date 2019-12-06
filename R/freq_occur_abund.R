#' Calculate abundance statistics
#'
#' @param otu A community matrix.
#'
#' @return A named list.
#' @importFrom vegan decostand
#' @importFrom stats sd
#' @keywords internal
freq_occur_abund <-
function(otu) {
  pa <- decostand(otu, "pa")
  freq.occur <- apply(pa, MARGIN=2, FUN=mean)
  otu <- decostand(otu, "total")
  mean.abund <- apply(otu, MARGIN=2, mean)
  mean.sd <- apply(otu, MARGIN=2, FUN=sd)
  sqrt.n <- sqrt(nrow(otu))
  mean.se <- mean.sd/sqrt.n
  rslt <- list(freq.occur, mean.abund, mean.sd, mean.se)
  names(rslt) <- c("freq.occur", "mean.abund", "mean.abund.sd", "mean.abund.se")
  return(rslt)
}
