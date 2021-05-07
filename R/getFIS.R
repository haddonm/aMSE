FISsettings <- list("MinIndexSize" = 1, # Lower size class in mm of the Index
                    "MaxIndexSize" = 210, # Upper size class of the index in mm, maybe relate to glb?
                    "FISType" = "Rough", # whether want good size frequency or bad "Rough"
                    "FISKnife" = TRUE, # Knife edged selectivity (TRUE) or logistic (FALSE)
                    "FISSel50" = rep(100, times=ctrl$projection), # FIS 50% selectivity, needs to be more elegant later
                    #If knife-edged will use this number
                    "FISSel95" = rep(140, times=ctrl$projection), # FIS 95% selectivity
                    "SFSampleSize" = rep(3,times=numpop), # NB NB will put this a FISSettinglist and make SampleSize by year and population. Need to chat how best to do this
                    "qFIS" = 0.01,  # FIS catchability related to numbers not biomass (needs discussion)
                    "FISBiasLower" = 1, # Lower value of uniform random for bias, bias of 1 is no bias
                    "FISBiasUpper" = 1, # Upper value for uniform random for bias
                    "ObsErrCompsMean" = 1, # Mean observation error on size frequency composition (normal distribution)
                    "ObsErrCompsSd" = 0.001, # Sd of observation error on s-f composition (normal distribution)
                    "WtaObs" = 2.413551e-05, # Observed length-weight a used in the HS
                    "WtbObs" = 3.385185e+00) # Observed length-weight b used in HS. NB at this stage assuming 1 global l-w regression


logisticCD <- function(inL50,delta,lens,knifeedge=0) {
  ans <- 1/(1+exp(-log(19.0)*(lens-inL50)/(delta)))
  if (knifeedge > 0) {
    pick <- which(lens < knifeedge)
    if (length(pick) > 0) {
      ans <- rep(1, times =length(lens))
      ans[1:pick[length(pick)]] <- 0.0
    } #if length(pick) >0
  }
  return(ans)
} # function logisticCD





#' Title getSF
#' @description Generates the average size frequency for a site which is randomly selected
#'     from the population of actual numbers by iteration. This returns site specific size-frequencies
#'
#' @param zoneDP the dynamic components of the simulated zone
#' @param glb the general global variables
#' @param ctrl the control settings
#' @param iter the iteration number
#' @param year the year to apply
#'
#' @return FISNSF a dataframe of dimension size, population, site
#' @export
#'
#' @examples
getSF <- function(zoneDP=zoneDP, glb=glb,ctrl=ctrl,FISsettings=FISsettings,iter=1, year=1) {
  Nclass <- glb$Nclass
  numpop <- glb$numpop
  BinSize <- glb$midpts[2]-glb$midpts[1] # this is also in zone1 but big object so hopefully saving time
  Nt <- zoneDP$Nt
  qFIS <- FISsettings$qFIS # bring in the catchability to scale the size-frequency
  MidSizeAll <- as.numeric(names(Nt[,1,1,1])) # Mid pt sizes of the observed data

  # Calculate the selectivity
  knifeedge <- 0
  if(FISsettings$FISKnife == TRUE) knifeedge <- FISsettings$FISSel50[year] # Test if knifeedged
  Select <- logisticCD(FISsettings$FISSel50[year],FISsettings$FISSel95[year]-FISsettings$FISSel50[year],
                     MidSizeAll, knifeedge=knifeedge) # Calculate selectivity. Is knife edged if knifeedge >0

  # set up the array and the calculate the observed size-frequency
  # the below can be done more elegantly
  FISNSF <-array(NA,dim = c(Size=Nclass,Pop=numpop,Site=max(SampleSize)))
  for (pop in 1:numpop) {
    FISBias <- runif(1,FISsettings$FISBiasLower,FISsettings$FISBiasUpper) # FISBias - Bias changes by year, iteration and population.
    for(site in 1:SampleSize[pop]) {
      IterationNumber <- round(runif(1,1,glb$reps)) # Calculate the iteration that will become the observed site's frequency
      if(FISsettings$FISType=="Smooth") {
        FISNSF[,pop,site] <- qFIS*FISBias*Select*Nt[,year,pop,IterationNumber] # Create the observed size frequency for each site
      } else if (FISsettings$FISType=="Rough") {
        # else Rough
        for (size in 1: Nclass) {
          FISNSF[size,pop,site] <- qFIS*FISBias*Select[size]*Nt[size,year,pop,IterationNumber] # randomly select at the size level
        } # for size

      } else { # if missed something
        print("One option not considered in FISType")
      } # else misses
    } # for site
  } #end for pop
  round(FISNSF) # round to nearest integer as a size-frequency of individual numbers
  return(FISNSF) # return the size frequency dataframe by size, pop and site within pop
} # end function getSF



#' Title getFIS Calculates the FIS data to simulate both a size frequency distribution presently at the zone level and
#'     an index of numbers and biomass. Will change the code to producing at the SAU level as think this is what is needed
#'     by all the Australian harvest strategies.
#'
#' @param zoneDP Dynamic operating model numbers
#' @param glb General global variables
#' @param ctrl Scenario specifi settings
#' @param FISsettings FIS specific settings
#' @param iter Iteration number
#' @param year Year in which the FIS is undertaken
#'
#' @return a list with the size-frequency in numbers by site, and the FIS numbers and biomass index
#' @export
#'
#' @examples
getFIS <- function(zoneDP=zoneDP, glb=glb,ctrl=ctrl,FISsettings = FISsettings,iter=1, year=1) {
  Nclass <- glb$Nclass
  BinSize <- glb$midpts[2]-glb$midpts[1] # this is also in zone1 but big object so hopefully saving time
  Nt <- zoneDP$Nt
  MidSizeAll <- as.numeric(names(Nt[,1,1,1])) # Mid pt sizes of the observed data

   # ToDo - need to extend the FIS back in time but getting it to work for proejction first

  #FISSF <- data.frame(matrix(NA, nrow=glb$Nclass, ncol=glb$numpop)) # Once I understand the changes MH is making will need to change these 2 lines
  ObsFISIndices <- data.frame(matrix(NA, nrow=glb$numpop,ncol=3))
  FISSF <- getSF(zoneDP,glb,ctrl,FISsettings,1,1) # this is the size frequency in observed numbers
  FISSFMean <- rowMeans(FISSF) # For Victoria this is density weighted. Need to check with SA. At this stage using simple mean across sites
  ObsFISN <- FISSFMean[MinIndexSizeBin:MaxIndexSizeBin] # select size that fall itno the index
  ObsFISNIndex <- sum(ObsFISN) # Sum so get index in numbers for the selected range

  # Calculate the index in biomass
  ObsFISB <- FISSFMean * FISsettings$WtaObs*glb$midpts^(FISsettings$WtbObs) # FIS biomass per bin for all sizes
  ObsFISBIndex <- sum(ObsFISB) # need to add kappa and error here
  FISIndices <- cbind(year,iter,ObsFISNIndex,ObsFISBIndex)

  FISData <- list(FISSF,FISIndices) # Put the sife frequency and indices together as a list()

} # function getFIS


SFOnly <- getSF(zoneDP,glb,ctrl,FISsettings,1,1) # this is the size frequency in observed numbers at the site level
DataFIS <- getFIS(zoneDP=zoneDP, glb=glb,ctrl=ctrl,FISsettings=FISsettings,iter=1, year=1) # this calculates both SF and Zone level Indices

# not sure why the "rough"version is still smooth
# code needs a through check
# havent added the kappa and other errors
