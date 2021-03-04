#' @name likelihoodMDC
#'
#' @title Likelihood calculation forthe RSF fitting procedure
#'
#' @description This is a function to calculate the likelihood. The function is used
#' in an iterative optimization process. At the end of the process, the likelihood, when
#' maximized, will provide thevalues for the coefficients from the RSF fitting.
#'
#' @param beta current values of the coefficients, size k1+k2
#'
#' @param X1 2-D matric of covariates for RSF model. Size should be ncells x k1
#'
#' @param X2 2-D matrix of covariates for the probability ofdetection. Size is
#' ncells x k2
#'
#' @param locations Observation matrix. The size is the number of scheduled fixes by 2
#'
#' @param k1 Number of coefficients in the RSF model
#'
#' @param k2 Number of coefficients in the probability of detection model
#'
#' @param maxLagArg Maximum number of matrix multiplications to try.
#'
#' @return This function returns a likelihood calculation
#'
#' @author John Lombardi (email: \code{jlombardi@west-inc.com})
#'
#' @export likelihoodMDC

likelihoodMDC <- function(beta, X1, X2, locations, k1, k2,maxLagArg){

  # -------------------------------------------------
  # Link function for probability of detection
  link <- function(x) {
    exp(x)/(1+exp(x))
  }

  rsf.coefs = beta[1:k1]
  p.coefs = t(t(beta[(k1+1):(k1+k2)]))
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
  #sumEXP = matrix(apply(EXP, 1, sum), nrow=ncells, ncol=ncells, byrow=F)
  sumEXP = matrix(rowSums(EXP), nrow=ncells, ncol=ncells)

  # browser()

  # evaluate p model
  P = link(X2%*%p.coefs)          #p is ncells X 1

  # Matrix of Pr[movement] x Pr[detection]
  D = (EXP / sumEXP) * matrix(P, nrow=ncells, ncol=ncells, byrow=T)
  # Matrix of Pr[movement] x (1 - Pr[detection])
  B = (EXP / sumEXP) * matrix(1-P, nrow=ncells, ncol=ncells, byrow=T)


  ##JAL
  ##Below is the most computational intensive part
  ##We can maybe just require Rcpp and do a function in-line call

  # Creating B^(m-2) when there are m-1 missing locations in a sequence
  if(max(Lag) > 1){
      B.s = vector("list", max(Lag-1))
      B.s[[1]] = diag(ncells)
      B.s[[2]] = B
      if(max(Lag) > 2){

        ##Use a maxlag
        for(j in 3:max(Lag-1)){

          output = B
          ##Below is for exp by squaring
          #Modify length(3:j)
          if(length(3:j) > maxLagArg)
          {
            numMults <- maxLagArg
          } else{

            numMults <- length(3:j)
          }



          B.s[[j]] <-MDC:::matrixCalcFastTwo(matrixOne = output,multTimes = numMults,
                                        matrixOfDiag = diag(1,ncells))

        }
      }

    }

    # Individual likelihood for each observation
    prob = rep(NA, times=fix.attempts)
    for(i in 2:fix.attempts){
      # If 2 consecutive scheduled GPS fixes were successful
      if(detected[i] & Lag[i] == 1){
        prob[i] = D[locations[i-1], locations[i]]
      }
      # If there is a gap
      if(detected[i] & Lag[i] > 1){
        prob[i] = B[locations[i-Lag[i]],] %*%
          B.s[[Lag[i]-1]] %*%
          D[,locations[i]]
      }
    }
    -sum(log(prob), na.rm=T)


  print(-sum(log(prob),na.rm=T))
}
