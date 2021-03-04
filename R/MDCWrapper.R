#' @name MDCWrapper
#'
#' @title Prepares data for the main MDC workhorse function
#'
#' @description This function takes inputs for locations, habitat, and covariates for the RSF and probability
#' of detection. A bit of data massaging goes on to convert dataframes to a more usable format for the
#' main MDC function.
#'
#' @param habitatDF Dataframe of habitat for the area of interest. Columns should include: coordinates,
#' a unique identifier, and relevant environmental or other covariates for each unit.
#'
#' @param CellID Dataframe of fixes for one animal. NOTE: For >1 animals this function will have to be
#' subsetted to one individual animal. Columns should include:a unique identifier that matches to the
#' habitatDF, a record of whether or not a fix was successful.
#'
#' @param nCells An integer value for the number of cells in the study area.
#'
#' @param rsfCovars A vector of strings that the RSF will be regressed on.
#' The values of the vector needs to match the names of columns in
#' the \code{habitatDF}. Example inputs: \code{c("elevation","percent.sage")}, \code{c("elevation")}.
#' As said before, these need to match the \code{colnames(habitatDF)}. If no environmental covariates are wanted,
#' e.g. an intercept-only model is the RSF model, don't provide any input to rsfCovars argument.
#'
#' @param probDetCovars A vector of strings for the probability of detection. The values of the vector needs to
#' match the names of columns in the \code{habitatDF}. Example inputs: \code{c("elevation","percent.sage")},
#' \code{c("elevation")}. As said before,
#' these need to match the \code{colnames(habitatDF)}. If the probability of detection will be an intercept-only
#' model, don't provide any input to the probDetCovars argument.
#'
#' @param distColumns A vector of strings for the name of the columns that contain coordinates in the
#' \code{habitatDF}. If distance is to be included as a covariate in the RSF, an extra step is needed
#' to calculate a distance matrix.
#'
#' @param maxLagArg An integer for the number of max matrix multiplications to perform. The default is 3.
#' Anything larger than 3 will greatly increase computation time.
#'
#' @return This function returns a MDC fit.
#'
#' @author John Lombardi (email: \code{jlombardi@west-inc.com})
#'
#' @export MDCWrapper


MDCWrapper <- function(habitatDF,CellID,nCells,
                       rsfCovars = NULL,probDetCovars = NULL,
                       distColumns = NULL,
                       maxLagArg = 3)
{

  #--------------------------------------------------------------------------------------------------
  #QAQC of inputs

  #Check for DF format
  if(!is.data.frame(habitatDF))
  {
    stop(paste("habitatDF should be a data frame."))
  }


  #--------------------------------------------------------------------------------------------------
  #Converting factors to be useful
  #Per fawn: if you have a 3 level factor we are splitting up that into 2 columns
  #Should be whatever base factor level is compared to the rest
  #For example: factor with 1, 2, and 3
  #Original col is 1, 1, 1, 2, 2, 2, 3, 3, 3
  #2 cols:
  #0 0 0 1 1 1 0 0 0
  #0 0 0 0 0 0 1 1 1

  ###NOTE: Below won't work if we have multiple factors (im pretty sure)

  ##NOTE we only want to do above if it's in the thing...

  #First look in habitatDF
  #Colnames to convert to multiple things
  factorHabitat <- sapply(habitatDF,is.factor)
  factorHabitat <- factorHabitat[!(factorHabitat == FALSE)]

  toConvertHabitat <- names(factorHabitat)

  factorDF <- as.data.frame(habitatDF[,colnames(habitatDF) %in% toConvertHabitat])

  if(ncol(factorDF) > 0){
    #Expand out each column and storein a list
    factorDFList <- list()
    for(i in 1:ncol(factorDF)){

      mm <- model.matrix(~factorDF[,i] - 1, model.frame(~factorDF[,i] - 1), contrasts = FALSE)

      factorDFList[[i]] <- cbind(mm[,2:ncol(mm)])

      #Really only works for >1 columns
      colnames(factorDFList[[i]]) <- gsub(".*]","",colnames(factorDFList[[i]]))

    }

    newDF <- do.call(cbind.data.frame,factorDFList)

    #modify rsfCovars
    #This next section ensures that the old names are changed from user inputs
    if(toConvertHabitat %in% rsfCovars){
      rsfCovars <- rsfCovars[!rsfCovars == toConvertHabitat]
      rsfCovars <- c(rsfCovars,colnames(newDF))
    }

    if(toConvertHabitat %in% probDetCovars){
      probDetCovars <- probDetCovars[!probDetCovars == toConvertHabitat]
      probDetCovars <- c(probDetCovars,colnames(newDF))
    }

    habitatDF <- habitatDF[,!(colnames(habitatDF) %in% toConvertHabitat)]
    habitatDF <- cbind(habitatDF,newDF)

  }

  ####



  #--------------------------------------------------------------------------------------------------
  #Calculate distance matrix here
  distMatrix = as.matrix(x = stats::dist(cbind(habitatDF[,colnames(habitatDF) %in% distColumns[1]],
                                      habitatDF[,colnames(habitatDF) %in% distColumns[2]])))


  #First check if rsfCovars is null
  if(is.null(rsfCovars)){

    #Throw an error
    stop("A covariate must be specific for the RSF")

  }else{

  #HabitatDF has our covariates on the habitat
    #NOTE: "distance" is never in the dataframe so we make it here
    #Wedo NOT need to get rid of it
  habitatDFSub <- as.data.frame(habitatDF[,colnames(habitatDF) %in% rsfCovars])

  #Convert each column to a matrix
  matrixList <- list()
  for(i in 1:ncol(habitatDFSub)){

    matrixList[[i]] <- matrix(habitatDFSub[,i], nrow=nCells, ncol=nCells, byrow=T)

  }

  if(!is.null(distMatrix) &
     (any(c("distance","DISTANCE","dist","DIST","Distance") %in% rsfCovars)))
  {
     matrixList[[length(matrixList) + 1]] <- distMatrix
  }

  #Modify matrixList to work with MDC
  #We need a vector of each element of the list
  vectorToConvert <- c()
  for(i in 1:length(matrixList)){

    vectorToConvert[i] <- paste("matrixList[[",i,"]]",sep="")

  }

  #Now conver to formula
  formulaRest <- paste("~",paste(vectorToConvert,collapse = "+"))

  selectionFormula <- as.formula(formulaRest)
  }##End else

  #probability of detection formula is a bit tricky
  #default is just intercept only
  #Otherwise it should be a vector

  if(is.null(probDetCovars)){

    probDetFormula <- NULL

  }else{

    logisticSub <- as.data.frame(habitatDF[,colnames(habitatDF) %in% probDetCovars])

    vectorToConvert <- c()
    for(i in 1:ncol(logisticSub)){

      vectorToConvert[i] <- paste("logisticSub[,",i,"]",sep="")

    }
      finalFormula <- paste("~ 1 +", paste(vectorToConvert,collapse = "+"))
      probDetFormula <- as.formula(finalFormula)

  }


  MDCFit <- MDC(selection = selectionFormula,
                p = probDetFormula,
                locations = CellID,
                ncells=nCells,
                maxLagArg = maxLagArg,
                RSFCovar = rsfCovars,
                LogCovar = c("Intercept",probDetCovars))

  return(MDCFit)

}
