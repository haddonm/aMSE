# The first version of creating a FIS will be to assume that
# each population will have a FIS and that this function is called in the projection part
# so called every year of the projection.
# Also no error as yet. Later I will change the control input files and slowly add
# all required changes throughout the code as selecting a population needs more thought.
# I have put all the FIS functions at the top so they can be easily moved if need be
#' Title GetFIS
#'
#' @description This generates the fishery independent data (FIS) required to run some of the
#'    harvest strategies. This code is still under development so do not use. At present assumes
#'    function is being called for each year of the projection and each iteration within that
#'    loop of DoTasProjection
#' @param zoneDDR
#' @param glb global variables
#' @param ctrl control settings
#' @param FISsettings FISsettings are settings required to calculate the index.
#'    Some of this should later be moved to one of csv files but presently in the first_use file.
#' @param iter replicate number
#' @param year year of the projection
#'
#' @return
#' @export
#'
#' @examples
getFIS <- function(zoneDP=zoneDP, glb=glb,ctrl=ctrl,FISsettings = FISsettings,iter=1, year=1) {
  Nclass <- glb$Nclass
  BinSize <- glb$midpts[2]-glb$midpts[1] # this is also in zone1 but big object so hopefully saving time
  Nt <- zoneDP$Nt
  MidSizeAll <- as.numeric(names(Nt[,1,1,1])) # Mid pt sizes of the observed data

  # NB NB replace this with truncated normal but first check if in mUtils as not in base or stats. this is a hack presently
  # Calculate selectivity
  knifeedge <- 0
  if(FISsettings$FISKnife == TRUE) knifeedge <- FISsettings$FISSel50[year] # Test if knifeedged
  Select <- logistic(FISsettings$FISSel50[year],FISsettings$FISSel95[year]-FISsettings$FISSel50[year],
                     MidSizeAll, knifeedge=knifeedge) # Calculate selectivity. Is knife edged if knifeedge >0
  MinIndexSizeBin <- ceiling(FISsettings$MinIndexSize/BinSize) # Min size of index (bin)
  MaxIndexSizeBin <- ceiling(FISsettings$MaxIndexSize/BinSize) # Max size of index (bin)

  # ToDo - need to extend the FIS back in time but getting it to work for proejction first

  FISSF <- data.frame(matrix(NA, nrow=glb$Nclass, ncol=glb$numpop)) # Once I understand the changes MH is making will need to change these 2 lines
  ObsFISIndices <- data.frame(matrix(NA, nrow=glb$numpop,ncol=3))

  for(pop in 1:glb$numpop) {

    FISBias <- runif(1,FISsettings$FISBiasLower,FISsettings$FISBiasUpper) # FISBias - Bias changes by year, iteration and population.
    # May need to change to truncated normal as more likely to centre around 1
    FISN <- Nt[,year,pop,iter] # Absolute size frequency without error, sampling in same year at the start of the year
    error <-runif(glb$Nclass, min=0,max=10)
    ObsFISNAct <- rowMeans(rmultinom(FISsettings$SFSampleSize,sum(FISN),FISN)) # generate a multinomial sample
    ObsFISN <- FISBias[iter]*FISsettings$qFIS*ObsFISNAct# we scale to catchability, bias
    # Note selectivity has still not been applied as worth code checking the numbers here

    #plot(glb$midpts,FISN) # Testing the multinomial, don't run
    #lines(glb$midpts,ObsFISNAct, lty=2) # Testing the multinomial, don't run
    #temp<- rmultinom(200,sum(FISN),FISN)
    #plot(glb$midpts[-1],temp[-1,1])
    #lines(glb$midpts[-1],temp[-1,2], lty=2)
    #difference <- temp[,1]-temp[,2]
    #plot(difference)


    ObsFISSFnFinal <- Select*ObsFISN # Now apply selectivity. This is the actual size-frequency the FIS will obtain
    # Need to check this as numbers strange

    ObsFISSFpFinal <- ObsFISSFnFinal/sum(ObsFISSFnFinal) # Size frequency as a proportion

    # NB I am proposing this is where we stop. If the HS needs an index then we give it ObsFISSFnFinal,
    # if the HS needs size-frequency but doesn't have an index, we give ObsFISSFpFinal
    # The HS can then cut its cloth as needed within its own code. So the below will be in the HS code rather.
    # It can arguably be here but if they need 2 indices then the person would do the above and cut in their HS. Same thing.
    # So from here on I think should be in the HS function e.g. SA HS.

    # creating a numbers index based on the cut-offs required for the HS.
    ObsFISNIndex <- ObsFISSFnFinal[MinIndexSizeBin:MaxIndexSizeBin]
    MidSize <- as.numeric(names(ObsFISNIndex)) # Mid pt sizes of the observed data
    ObsFISBIndex <- ObsFISNIndex * FISsettings$WtaObs*MidSize^(FISsettings$WtbObs) # FIS biomass per bin, NB Need to check that the sizes are commensurate to Nt
    ObsFISNumber <- sum(ObsFISSFnFinal) # Index in numbers for the selected range
    ObsFISBiomass <- sum(ObsFISBIndex) # Actual absolute weight except for s-f error
    temp <- cbind(pop,ObsFISNumber,ObsFISBiomass)
    temp2 <- cbind(pop,ObsFISSFnFinal,ObsFISSFpFinal)
    # IndexInB = q * qyear * B^gamma*exp(error) NB Still need to add gamma and exp error
    if(pop == 1) {
      ObsFISIndices <- temp
      ObsFISSF <- temp2
    } else {
      ObsFISIndices <- rbind(ObsFISIndices,temp) # put all the population data together
      ObsFISSF <- rbind(ObsFISSF, temp2)
    } # if pop

    # Still not complete at all! And not checked but now think it is time wasted as MH changing code
    # to much for this to be relevant
  } # for pop

  return(list(ObsFISIndices, ObsFISSF))
} # function getFIS

getSF <- function(iter=1,year=1){
  Nclass <- glb$Nclass
  BinSize <- glb$midpts[2]-glb$midpts[1] # this is also in zone1 but big object so hopefully saving time
  Nt <- zoneDP$Nt
  MidSizeAll <- as.numeric(names(Nt[,1,1,1])) # Mid pt sizes of the observed data
  SampleSize <- 20
  FISN2 <-matrix(glb$Nclass)
  for(size in 1:glb$Nclass) {
    FISN2[size] <- sample(Nt[size,year,pop,],SampleSize,replace=TRUE,prob = NULL)
  }


  SF <- NULL
  return(SF)
}
