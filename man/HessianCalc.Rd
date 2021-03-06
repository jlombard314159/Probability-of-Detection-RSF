% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HessianCalc.R
\name{HessianCalc}
\alias{HessianCalc}
\alias{F.2nd.deriv}
\title{F.2nd.deriv}
\usage{
F.2nd.deriv(x, FUN, ...)
}
\arguments{
\item{x}{The location (a vector) where the second derivatives
of \code{FUN} are desired.}

\item{FUN}{An R function to compute second derivatives for.
This must be a function of the form FUN <- function(x, ...){...}
where x is the parameters of the function (e.g., location \code{x}),
and ... are additional paramters needed to evaluate the function.
FUN must return a single value (scalar), the height of the
surface above \code{x}, i.e., FUN evaluated at x.}

\item{...}{Additional arguments}
}
\description{
Computes numeric 2nd derivatives (hessian) of an
arbitrary multidimensional function at a particular location.
}
\details{
This function uses the "5-point" numeric second derivative
method advocated in numerous numerical recipie texts.  During computation
of the 2nd derivative, FUN will be evaluated
capable of evaluating at locations within a hyper-elipsoid
with cardinal radii 2*\code{x}*(eps)^0.25, where eps=10e-7.
If those radii don't
work for your function, make eps a parameter to this function and
change it, possibly varying by dimension.

A handy way to use this function is to call an optimization routine
like \code{nlminb} with FUN, then call this function with the
optimized values (solution) and FUN.  This will yeild the hessian
at the solution rather than the hessian at the previous step of the
optimization.
}
\author{
Trent McDonald
}
