#' @name getModMatrix
#'
#' @title Prepares data for the main MDC workhorse function
#'
#' @description Function to return the model matrix from a formula object.  I.e., convert
# ~ x1 + x2 into a usable model matrix.
#'
#' @param f  an R formula object, without a response.  I.e., of the form
#'  "~ x1 + x2 + ..."
#'
#' @param ncells number of cells in study area = row dimension of matrices
#'
#'
#' @return A matrix with variables in formula as columns.  Also included is
#' the names of variables in the formula and whether the model has an intercept.
#'
#' @author John Lombardi (email: \code{jlombardi@west-inc.com})
#'
#' @export getModMatrix



getModMatrix <- function(f, ncells){
  # browser()

  call <- match.call()
  contrasts <- NULL
  mf <- match.call(expand.dots = FALSE)
  mf$family <- mf$start <- mf$control <- mf$maxit <- NULL
  mf$model <- mf$method <- mf$x <- mf$y <- mf$contrasts <- NULL
  mf$... <- NULL
  mf$f <- mf$ncells <- NULL
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf$formula <- f
  form <- as.character(formula(mf))[2]
  intercept <- NA
  if(nchar(form) == 1 & form == "1"){
    # this is a model with intercept only
    X <- matrix(1, ncells, 1)
    intercept <- TRUE
    nx <- 1
    x.names <- character(0)
  } else {  #This else is if for more than just an intercept only model
    if(substr(form, start = 1, stop = 1) == "1"){
      int <- 1
      intercept <- TRUE
    } else {
      int <- 0
    }

    #mf is the model form - specifically the data from the formula
    #So for log sub it is a ncells x n covars
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    if(int == 0) attr(mt, "intercept") <- 0

    # ##### #Addin factor
    #
    # #Look for where there are factors in MT
    # factorsToFind <- as.data.frame(attr(mt, "dataClasses"))
    # colnames(factorsToFind) <- c("Type")
    #
    # factorsToFind$DataName <- row.names(factorsToFind)
    #
    # factorsToFind <- factorsToFind[factorsToFind$Type == "factor",]
    #
    # #Now grab the row names
    # contrastsToAdd <-factorsToFind$DataName
    #
    # #Now Mke a list or make it null
    #
    # if(identical(contrastsToAdd, character(0)))
    # {
    #   contrasts <- NULL
    # } else {
    #   # browser()
    #   contrasts <- list()
    #   for(i in 1:length(contrastsToAdd)){
    #
    #     contrasts[[i]] <- "contr.treatment"
    #   }
    #
    #   names(contrasts) <- contrastsToAdd
    #
    # }

    #Get rid of first col which is intercept
    xvars <- as.character(attr(mt, "variables"))[-1]

    #JAL: This is getting rid of the resposne variable if it is in there
    if ((yvar <- attr(mt, "response")) > 0) xvars <- xvars[-yvar]

    #No clue what this is doing - its not used anywhere after this
    xlev <- if (length(xvars) > 0) {
      xlev <- lapply(mf[xvars], levels)
      xlev[!sapply(xlev, is.null)]
    }

    X <- NA
    X <- if (length(attr(mt, "order")) != 0) {
      model.matrix(mt, mf,contrasts)

    } else {
      stop("Empty models not allowed")
    }
    assign.col <- attr(X, "assign")
    nx <- length(unique(assign.col))
    x.names <- xvars
    dimnames(X) <- NULL
  }
  ans <- list(
    X = X,
    n.covars = nx,
    intercept = intercept,
    vars = x.names)
  return(ans)
}
