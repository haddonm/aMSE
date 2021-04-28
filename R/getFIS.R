FISsettings <- list("MinIndexSize" = 1, # Lower size class in mm of the Index
                    "MaxIndexSize" = 210, # Upper size class of the index in mm, maybe relate to glb?
                    "FISKnife" = FALSE, # Knide edged selectivity (TRUE) or logistic (FALSE)
                    "FISSel50" = rep(100, times=ctrl$projection), # FIS 50% selectivity, needs to be more elegant later
                    #If knife-edged will use this number
                    "FISSel95" = rep(140, times=ctrl$projection), # FIS 95% selectivity
                    "SFSampleSize" = 10, # effective sample size for multinomial distribution to generate size-frequency, NB This one needs to be expanded
                    "qFIS" = 0.01,  # FIS catchability related to numbers not biomass (needs discussion)
                    "FISBiasLower" = 1, # Lower value of uniform random for bias, bias of 1 is no bias
                    "FISBiasUpper" = 1, # Upper value for uniform random for bias
                    "ObsErrCompsMean" = 1, # Mean observation error on size frequency composition (normal distribution)
                    "ObsErrCompsSd" = 0.001, # Sd of observation error on s-f composition (normal distribution)
                    "WtaObs" = 2.413551e-05, # Observed length-weight a used in the HS
                    "WtbObs" = 3.385185e+00) # Observed length-weight b used in HS. NB at this stage assuming 1 global l-w regression

#' Title getSF
#' @description Generates the size frequency which is drawn based on the samples size from the population
#'     of actual numbers by iteration. The more you sample, the more the size-frequency will look like
#'     the average across all the iteration. This is not the same as the actual or perfect knowledge
#'
#' @param zoneDP the dynamic components of the simulated zone
#' @param glb the general global variables
#' @param ctrl the control settings
#' @param iter the iteration number
#' @param year the year to apply
#'
#' @return FISNSF a dataframe of dimension size and population
#' @export
#'
#' @examples
getSF <- function(zoneDP=zoneDP, glb=glb,ctrl=ctrl,FISsettings=FISsettings,iter=1, year=1) {
  Nclass <- glb$Nclass
  BinSize <- glb$midpts[2]-glb$midpts[1] # this is also in zone1 but big object so hopefully saving time
  Nt <- zoneDP$Nt
  MidSizeAll <- as.numeric(names(Nt[,1,1,1])) # Mid pt sizes of the observed data
  SampleSize <- rep(FISsettings$SFSampleSize,times=glb$numpop) # NB NB will put this a FISSettinglist and make by year and population. Need to chat how best to do this
  FISNSF <-data.frame(matrix(NA, nrow=glb$Nclass, ncol=glb$numpop))
  for (pop in 1:glb$numpop) {
    for(size in 1:glb$Nclass) {
      FISNSF[size,pop] <- round(mean(sample(Nt[size,year,pop,],SampleSize[pop],replace=TRUE,prob = NULL)))
    } # for size
  } #end for pop
  return(FISNSF) # return the size frequency dataframe by size and pop
} # end function getSF

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


  # Calculate selectivity
  knifeedge <- 0
  if(FISsettings$FISKnife == TRUE) knifeedge <- FISsettings$FISSel50[year] # Test if knifeedged
  Select <- logistic(FISsettings$FISSel50[year],FISsettings$FISSel95[year]-FISsettings$FISSel50[year],
                     MidSizeAll, knifeedge=knifeedge) # Calculate selectivity. Is knife edged if knifeedge >0
  MinIndexSizeBin <- ceiling(FISsettings$MinIndexSize/BinSize) # Min size of index (bin)
  MaxIndexSizeBin <- ceiling(FISsettings$MaxIndexSize/BinSize) # Max size of index (bin)

  # ToDo - need to extend the FIS back in time but getting it to work for proejction first

  #FISSF <- data.frame(matrix(NA, nrow=glb$Nclass, ncol=glb$numpop)) # Once I understand the changes MH is making will need to change these 2 lines
  ObsFISIndices <- data.frame(matrix(NA, nrow=glb$numpop,ncol=3))
  FISSF <- getSF(zoneDP,glb,ctrl,FISsettings, 1,1) # size frequency rows by population in columns
  for(pop in 1:glb$numpop) {

    FISBias <- runif(1,FISsettings$FISBiasLower,FISsettings$FISBiasUpper) # FISBias - Bias changes by year, iteration and population.
    # May need to change to truncated normal as more likely to centre around 1


    # Preserving the multinomial code in case I come back to this, dlete later if not needed
    ###################
    #FISN <- Nt[,year,pop,iter] # Absolute size frequency without error, sampling in same year at the start of the year
    #ObsFISNAct <- rowMeans(rmultinom(FISsettings$SFSampleSize,sum(FISN),FISN)) # generate a multinomial sample
    #ObsFISN <- FISBias[iter]*FISsettings$qFIS*ObsFISNAct# we scale to catchability, bias
    # Note selectivity has still not been applied as worth code checking the numbers here

    #plot(glb$midpts,FISN) # Testing the multinomial, don't run
    #lines(glb$midpts,ObsFISNAct, lty=2) # Testing the multinomial, don't run
    #temp<- rmultinom(200,sum(FISN),FISN)
    #plot(glb$midpts[-1],temp[-1,1])
    #lines(glb$midpts[-1],temp[-1,2], lty=2)
    #difference <- temp[,1]-temp[,2]
    #plot(difference)

    ############################# Stops here, delete above

    ObsFISN <- FISBias[iter]*FISsettings$qFIS*FISSF[,pop] #Scale to catchability and bias
    ObsFISSFnFinal <- Select*ObsFISN # Now apply selectivity. This is the actual size-frequency the FIS will obtain

    # NB Need to include exp(error)
    ObsFISSFpFinal <- ObsFISSFnFinal/sum(ObsFISSFnFinal)# Size frequency as a proportion

    ObsFISSFpop <- cbind(as.data.frame(glb$midpts),ObsFISSFnFinal,ObsFISSFpFinal)

    # Creating a numbers index based on the cut-offs required for the HS.
    ObsFISNIndex <- ObsFISSFnFinal[MinIndexSizeBin:MaxIndexSizeBin]
    ObsFISNumber <- sum(ObsFISSFnFinal) # Index in numbers for the selected range

    # Creating a biomass index from the same cutoffs

    ObsFISB <- ObsFISSFpop[,"ObsFISSFnFinal"] * FISsettings$WtaObs*ObsFISSFpop[,"glb$midpts"]^(FISsettings$WtbObs) # FIS biomass per bin for all sizes
    ObsFISBiomass <- sum(ObsFISB[MinIndexSizeBin:MaxIndexSizeBin])

    temp <- cbind(iter,year,pop,ObsFISNumber,ObsFISBiomass)
    temp2 <- cbind(iter,year,pop,ObsFISSFpop)
    # IndexInB = q * qyear * B^gamma*exp(error) NB Still need to add gamma and exp error
    if(pop == 1) {
      ObsFISIndices <- temp
      ObsFISSF <- temp2
    } else {
      ObsFISIndices <- rbind(ObsFISIndices,temp) # put all the population data together
      ObsFISSF <- rbind(ObsFISSF, temp2)
    } # if pop

    # Needs a deep check. Also size frequency is more noisy, but enough?
  } # for pop

  return(list(ObsFISIndices, ObsFISSF))
} # function getFIS

test <- getFIS(zoneDP=zoneDP, glb=glb,ctrl=ctrl,FISsettings = FISsettings,iter=1, year=1)

# Note to self: ObsFISBiomass does not seem to have worked. Check this as worked before
