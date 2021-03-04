#' @name likelihoodMDCRSF
#'
#' @title Likelihood calculation forthe RSF fitting procedure
#'
#' @description This is a function to calculate the likelihood. The function is used
#' in an iterative optimization process. At the end of the process, the likelihood, when
#' maximized, will provide thevalues for the coefficients from the RSF fitting. This is used when
#' we do NOT have a probability of detection argument (100% fix success)
#'
#' @param beta current values of the coefficients, size k1+k2
#'
#' @param X1 2-D matric of covariates for RSF model. Size should be ncells x k1
#'
#'
#' @param locations Observation matrix. The size is the number of scheduled fixes by 2
#'
#' @param k1 Number of coefficients in the RSF model
#'
#' @return This function returns a likelihood calculation
#'
#' @author John Lombardi (email: \code{jlombardi@west-inc.com})
#'
#' @export likelihoodMDCRSF

likelihoodMDCRSF <- function(beta, X1,  locations, k1){

  # -------------------------------------------------
  rsf.coefs = beta[1:k1]
  ncells = nrow(X1)

  # Number of fix attempts
  fix.attempts = length(locations)

  # Lag between observations.  Should be 1 when no locations are missing.

  ##JAL EDIT:This locations needs to be a vector
  Lag = rep(1, times=fix.attempts)
  for(i in 3:fix.attempts){
    if(is.na(locations[i-1])) Lag[i] = Lag[i-1] + 1
  }

  # Detected (Yes = 1; No = 0)
  detected = !is.na(locations)

  # evaluate movement model
  EXP = matrix(0, nrow=ncells, ncol=ncells)
  for(i in 1:k1){
    EXP = EXP + X1[,(1+(i-1)*ncells):(i*ncells)]*rsf.coefs[i]
  }
  EXP = exp(EXP)
  sumEXP = matrix(rowSums(EXP), nrow=ncells, ncol=ncells)

  # Matrix of Pr[movement]
  D = (EXP / sumEXP)

  # Individual likelihood for each observation
  prob = rep(NA, times=fix.attempts)
  for(i in 2:fix.attempts){
    # If 2 consecutive scheduled GPS fixes were successful
    if(detected[i] & Lag[i] == 1){
      prob[i] = D[locations[i-1], locations[i]]
    }

  }
  -sum(log(prob), na.rm=T)


  print(-sum(log(prob),na.rm=T))
}
