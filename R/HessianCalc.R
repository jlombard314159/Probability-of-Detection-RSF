#' @name HessianCalc
#
#' @title F.2nd.deriv
#'
#' @description Computes numeric 2nd derivatives (hessian) of an
#' arbitrary multidimensional function at a particular location.
#'
#' @param x The location (a vector) where the second derivatives
#' of \code{FUN} are desired.
#'
#' @param FUN An R function to compute second derivatives for.
#' This must be a function of the form FUN <- function(x, ...){...}
#' where x is the parameters of the function (e.g., location \code{x}),
#' and ... are additional paramters needed to evaluate the function.
#' FUN must return a single value (scalar), the height of the
#' surface above \code{x}, i.e., FUN evaluated at x.
#'
#' @param ... Additional arguments
#'
#' @details This function uses the "5-point" numeric second derivative
#' method advocated in numerous numerical recipie texts.  During computation
#' of the 2nd derivative, FUN will be evaluated
#' capable of evaluating at locations within a hyper-elipsoid
#' with cardinal radii 2*\code{x}*(eps)^0.25, where eps=10e-7.
#' If those radii don't
#' work for your function, make eps a parameter to this function and
#' change it, possibly varying by dimension.
#'
#' A handy way to use this function is to call an optimization routine
#' like \code{nlminb} with FUN, then call this function with the
#' optimized values (solution) and FUN.  This will yeild the hessian
#' at the solution rather than the hessian at the previous step of the
#' optimization.
#'
#' @author Trent McDonald


F.2nd.deriv <- function(x, FUN, ...){
  FUN <- match.fun(FUN)
  d <- length(x)   # number of dimensions
  hess <- matrix(0, nrow=d, ncol=d)
  eps <- 10e-7
  h <- ifelse(x==0, eps^0.25, (eps^(0.25))*x )
  for(i in 1:d){
    ei <- rep(0,d)
    ei[i] <- 1
    # compute diagonal element
    hess[i,i] <- (-FUN(x+2*h*ei, ...) + 16*FUN(x+h*ei, ...) - 30*FUN(x, ...) +
                    16*FUN(x-h*ei, ...) - FUN(x-2*h*ei, ...)) / (12*h[i]*h[i])
    if((i+1) <= d){
      for(j in (i+1):d){
        ej <- rep(0,d)
        ej[j] <- 1
        # compute off diagonal element
        hess[i,j] <- (FUN(x+h*ei+h*ej, ...) - FUN(x+h*ei-h*ej, ...) -
                        FUN(x-h*ei+h*ej, ...) + FUN(x-h*ei-h*ej, ...)) / (4*h[i]*h[j])
        # Assume symetric
        hess[j,i] <- hess[i,j]
      }
    }
  }
  hess
}
