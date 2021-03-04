#' @name MDC
#'
#' @title Fits the discrete choise RSF
#'
#' @description This is a function to fit the modified
#' discrete choice RSF to animal locations
#' when there is <100% fix success, assuming fix success is a function of habitat.
#'
#' @param selection model formula for rsf parameters (the relative probability of
#' selection parameters).
#' Covariates in this model are matrices of size ncells X ncells.
#' Value of element [i,j]
#' of covariate k in this formula is the kth covariate value for the jth cell.
#' In this formulation the rsf parameters are not allowed to vary over time.  Any
#' distance related covariate (e.g, distance from previous cell) should be a matrix with
#' value of element [i,j] being the distance between the ith and jth cells.
#' No intercept for this model.
#'
#' @param p model formula for probability of detection given presence parameters.  Covariates
#' in this mocel are vectors of size ncells.  Value in element [i] of variable
#' k is value of kth covariate for ith site.  Intercept (1) must be present.
#'
#' @param locations Cell ID # where the animal was located
#' (=NA if the animal was not found) during each fix attempt. Cell ID
#' should be numeric, the minimum possible is 1, and the max = ncells.
#'
#' @param ncells number of available cells in the study area.
#'
#' @param maxLagArg Integer amount of maximum number of matrix multiplications to perform.
#'
#' @param RSFCovar Characer string for RSF covariates
#'
#' @param LogCovar Character vector for the covariates for logistic regression
#'
#' @return This function returns a MDC fit.
#'
#' @author John Lombardi (email: \code{jlombardi@west-inc.com})
#'
#' @export MDC


MDC <- function(selection, p, locations, ncells,maxLagArg,RSFCovar,
                LogCovar){

  # -------- Main code for siteocc.fn ------
  if(missing(selection)) stop("Model for habitat selection must be specified")
  if(missing(p)) stop("Model for p (encounter probabilities) must be specified")
  if(missing(locations)) stop("Animal locations must be specified")
  if(max(locations, na.rm=T) > ncells) stop("Cell ID #s are not correct")
  if(!is.numeric(locations)) stop("Cell ID #s are not numeric")

  orig.call <- match.call()
  # browser()
  # memory.limit(6.291e+6)
  rsf.mod.mat <- getModMatrix(selection, ncells)

  #Define for nlminb - these will never be null
  X.rsf <- rsf.mod.mat$X
  k.rsf <- rsf.mod.mat$n.covars
  # browser()
  #May have to check tosee if it is not null
if(is.null(p))
  {
    #Do maximizaiton here
    strt.vals <- rep(0, k.rsf)

    #   Do maximization
    out <- nlminb(start = strt.vals, objective = likelihoodMDCRSF,
                  X1=X.rsf,
                  locations=locations, k1=k.rsf)



    # Function call
    hessian = F.2nd.deriv(out$par, likelihoodMDCRSF, X1=X.rsf,
                          locations=locations, k1=k.rsf)
    SEs = sqrt(diag(solve(hessian)))

    rsf.coefs <- out$par[1:k.rsf]
    RSFDataframe <- data.frame("Covar" = RSFCovar,
                               "Coef" = rsf.coefs,
                               "SE" = SEs[1:k.rsf])

    ans <- list(    loglik = -out$objective,
                    convergence = out$convergence,
                    call = orig.call,
                    ncells = ncells,
                    n.fix.attempts = length(locations),
                    RSF = RSFDataframe,
                    aic = 2*out$objective  + 2*(k.rsf),
                    bic = 2*out$objective + (k.rsf)*log(sum(!is.na(locations))-1))
    class( ans ) <- "rsf"

    return(ans)


}

  ##IF we have a probability of detection model
if(!is.null(p)){


    p.mod.mat   <- getModMatrix(p, ncells)


    X.p   <- p.mod.mat$X
    k.p   <- p.mod.mat$n.covars  # includes intercept if present

    strt.vals <- rep(0, k.rsf+k.p)

    #   Do maximization
    out <- nlminb(start = strt.vals, objective = likelihoodMDC,
                  X1=X.rsf, X2=X.p,
                  locations=locations, k1=k.rsf, k2=k.p,
                  maxLagArg = maxLagArg)


    # Function call
    hessian = F.2nd.deriv(out$par, likelihoodMDC, X1=X.rsf, X2=X.p,
                          locations=locations, k1=k.rsf, k2=k.p,
                          maxLagArg = maxLagArg)
    SEs = sqrt(diag(solve(hessian)))


    rsf.coefs <- out$par[1:k.rsf]
    p.coefs <- out$par[(k.rsf+1):(k.rsf+k.p)]

    RSFDataframe <- data.frame("Covar" = RSFCovar,
                               "Coef" = rsf.coefs,
                               "SE" = SEs[1:k.rsf])

    LogDataframe <- data.frame("Covar" = LogCovar,
                               "Coef" =  p.coefs,
                               "SE" = SEs[(k.rsf+1):(k.rsf+k.p)])

    ans <- list(    loglik = -out$objective,
                    convergence = out$convergence,
                    call = orig.call,
                    ncells = ncells,
                    n.fix.attempts = length(locations),
                    RSF = RSFDataframe,
                    ProbOfDet = LogDataframe,
                    aic = 2*out$objective  + 2*(k.rsf+k.p),
                    bic = 2*out$objective + (k.rsf+k.p)*log(sum(!is.na(locations))-1))
    class( ans ) <- "rsf"

    return(ans)
  }
}


