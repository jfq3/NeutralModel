#' Calculate Ntm parameter
#'
#' Finds the Ntm parameter by reiteration.
#'
#' @param obs The observed commmunity matrix.
#' @param Ntm A starting sequence of integers used in finding Ntm.
#' @details Calculates the expected probability of detection of each OTU in
#'  the target community if it was present merely because of neutral
#'   processes from the source community.
#'   These probabilities are calculated for all values of Ntm and then
#'   the best fit candidate is determined based upon a least sum of squares
#'   approach. The function is called recursively until Ntm is estimated
#'   to 2 decimal places.
#' @return A list consiting of Ntm and the fit neutral matrix.
#' @importFrom stats pbeta
#' @keywords internal
calc_Ntm <-
function(obs, Ntm = seq(1,100,1)){
  detlim <- obs[nrow(obs),4]
  objfun <- matrix(nrow=length(Ntm),ncol=2)
  fit_neutral_matrix <- matrix(nrow=nrow(obs),ncol=4)
  for (j in 1:length(Ntm)){
    Ntmloop <- Ntm[j]
    fit_neutral_matrix[,1] <-  Ntmloop*obs[,4]
    fit_neutral_matrix[,2] <-  Ntmloop*(1-obs[,4])
    fit_neutral_matrix[,3] <- 1-(pbeta(detlim,fit_neutral_matrix[ ,1],fit_neutral_matrix[ ,2]))
    fit_neutral_matrix[,4] <- (obs[,1]-fit_neutral_matrix[,3])^2
    objfun[j,1] <- sum(fit_neutral_matrix[,4])
    objfun[j,2] <- Ntmloop
  }

  # Determine the best minimum objective function
  minobj <- which.min(objfun[,1])
  minobjfun <- objfun[minobj[1],1]

  # Determine the best Ntm value
  bestNtm <- objfun[minobj[1],2]

  if(bestNtm == range(Ntm)[1]) {
      Ntm1 <- bestNtm/10
      Ntm2 <- 2*bestNtm
      intv <- 1/(2*(Ntm2/Ntm1))
      Ntm  <- seq(Ntm1, Ntm2, intv)
      calc_Ntm(obs, Ntm)
    } else if (bestNtm == range(Ntm)[2]) {
      Ntm1 <- bestNtm*0.75
      Ntm2 <- bestNtm*1.25
      intv <- 1/(2*(Ntm2/Ntm1))
      Ntm <- seq(Ntm1, Ntm2, intv)
      calc_Ntm(obs, Ntm)
    } else if ((nchar(bestNtm-trunc(bestNtm))>2)) {
      # return(bestNtm)
      rslt <- list(bestNtm = bestNtm, fit_neutral_matrix = fit_neutral_matrix)
      return(rslt)
    }
      else {
      Ntm1 <- bestNtm*0.75
      Ntm2 <- bestNtm*1.25
      intv <- 1/(2*(Ntm2/Ntm1))
      Ntm <- seq(Ntm1, Ntm2, intv)
      calc_Ntm(obs, Ntm)
    }

}
