## ------------------------------------------------------------------------

#habitatDF is a dataframe that contains environmental covariates for the RSF

habitatDF <- data.frame(UTMX = c(203453,201353,205122,230131,205315),
                        UTMY = c(4534122,4134124,4342313,4314431,423434),
                        UnitID = c(1,2,3,4,5),
                        Elevation = c(2.12,2.34,2.24,2.24,2.56))

#Total number of cells in the study area
nCells = max(habitatDF$UnitID)

#locationsDF is a dataframe of fix attempts for an individual collared animal

locationsDF <- data.frame(UnitID = c(3,4,2,4,NA,NA,1),
                          FixAttempt = c(1,2,3,4,5,6,7))

##For inputs, we only need the record of fix/not fixed or the UnitID column
# test<- MDCWrapper(habitatDF = habitatDF,
#                        locationsDF = locationsDF$UnitID,
#                        nCells = nCells,
#                        rsfCovars = c("Elevation"),
#                        distColumns  = c("UTMX","UTMY"))



