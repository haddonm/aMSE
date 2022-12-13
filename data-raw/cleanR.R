

#' @title blockE13 is an abalone data-set for testing performance measures
#'
#' @description blockE13 is a fishery data-set for blacklip abalone
#'     (\emph{Haliotis rubra}) from block 13 in the eastern zone, this is
#'     Tasmania's Block 13. It constitutes three time-series of the same cpue
#'     data in different formats (see below). It is for use when testing the
#'     performance measures within Tasmania's MCDA, although it could be used
#'     for other purposes, such as illustrating the typical linear relationship
#'     between catch and CPUE. Note there will only ever be one less PM value
#'     than there are cpue data in the time-series.
#'
#' @name blockE13
#'
#' @docType data
#'
#' @format A data.frame of abalone fishery data
#' \describe{
#'   \item{year}{the year of fishing}
#'   \item{coeff}{the back-transformed coefficients from a standardization}
#'   \item{scaled}{the coefficients scaled to a mean of one}
#'   \item{cpue}{the coefficients scaled to the nominal geometric mean of the
#'       time series, to place it on the nominal scale}
#'   \item{catch}{the related catch from the eastern parts of block 13}
#' }
#'
#' @section Subjects:
#'  \itemize{
#'    \item performance measure estimation
#'    \item relationship between catch and CPUE
#'  }
#'
#' @source Mundy, C. and J.M. McAllister (2020) Tasmanian Abalone Assessment
#'     2019. IMAS, University of Tasmania.
#'
#' @examples
#' data(blockE13)
#' blockE13
NULL


#' @title lf10 contains three years of length-composition data for block 10
#'
#' @description lf10 is a data.frame of length-composition of commercial catch
#'     from block 10 on the west coast of Tasmania in a longdat format as
#'     derived from the commlf function makelongdat.This makes it useful for
#'     illustrating the use of makewidedat.
#'
#' @name lf10
#'
#' @docType data
#'
#' @format A data.frame of abalone size-composition data
#' \describe{
#'   \item{year}{the year of sampling}
#'   \item{sau}{the block or SAU}
#'   \item{length}{the measured size in mm, with 2mm size classes}
#'   \item{counts}{the counts at each length class}
#'   \item{propcounts}{the proportion of each length class of the total}
#' }
#'
#' @section Subjects:
#'  \itemize{
#'    \item selectivity
#'    \item size-based stock assessment model fitting
#'    \item size-distributions of commercial catches
#'  }
#'
#' @source Thanks to the Institute of Marine and Antarctic Science,
#'     which is part of the University of Tasmania, and especially to
#'     Dr Craig Mundy, leader of the Abalone Group, for permission to use
#'     this data.
#'
#' @examples
#'  data(lf10)
#'  head(lf10,20)
#'  mids <- seq(138,210,2)
#'  answer <- makewidedat(lf10,mids)
#'  answer[1:20,]
NULL


#' @title NOTdodepletion resets zoneD to an input depletion level
#'
#' @description Notdodepletion resets the depletion level of the whole
#'     zone and does this by searching for the harvest rate that
#'     leads to the sum of the mature biomass, across populations,
#'     divided by the sum of the effective B0 across populations is
#'     as close as possible to the desired depletion level. This means
#'     the individual populations will likely vary around the target
#'     depletion, but across the zone it will be correct. The
#'     depletion is measured relative to the effective B0 as that
#'     takes account of any larval dispersal. This function uses the
#'     production curve matrix to search for harvest rates that bound
#'     the target depletion and then re-searches across those bounds
#'     using len intervals.
#'
#' @param zoneC the constants components of the simulated zone
#' @param zoneD the dynamic components of the simulated zone
#' @param glob the general global variables
#' @param depl the target depletion proportion for the whole zone
#' @param product the production curve matrix from doproduction
#' @param len the number of intervals in the trial harvest rates
#'
#' @return a revised zoneD object
#' @export
#'
#' @examples
#' \dontrun{
#'   data(zone1)
#'   glb <- zone1$globals
#'   data(constants)
#'   ans <- makezoneC(zone1,constants)
#'   zoneC <- ans$zoneC
#'   glb <- ans$glb  # now contains the move matrix
#'   ans <- makezone(glb,zoneC)
#'   zoneC <- ans$zoneC
#'   zoneD <- ans$zoneD
#'   ans2 <- modzoneC(zoneC,zoneD,glb)
#'   zoneC <- ans2$zoneC  # now has MSY and deplMSY
#'   product <- ans2$product
#'   zoneDD <- Notdodepletion(zoneC,zoneD,glb,depl=0.3,product)
#'   sum((zoneDD$matureB[1,]/sum(zoneDD$matureB[1,]))*zoneDD$deplsB[1,])
#'   mean(zoneDD$deplsB[1,])
#' }
Notdodepletion <- function(zoneC,zoneD,glob,depl,product,len=15) {
  #  zoneC=zoneC; zoneD=zoneD; glob=glb;  product=zone$product; depl=0.2;len=15
  # use product to find bounds on H around the desired depletion level
  spb <- rowSums(product[,"MatB",])
  initH <- as.numeric(rownames(product))
  zoneDyn <- cbind(initH,spb,spb/spb[1])
  colnames(zoneDyn) <- c("harv","spb","deplet")
  regB0 <- spb[1]
  pick <- which.closest(depl,zoneDyn[,"deplet"])
  if (abs(zoneDyn[pick,"deplet"] - depl) > 0.05)
    warning("Resolution of production curves may be insufficient.")
  if (pick == nrow(product)) {
    mssg <- paste0("Initial maximum harvest rate when estimating ",
                   "productivity was too low to generate desired ",
                   "depletion level")
    stop(mssg)
  }
  lowl <- initH[pick-1]
  upl <- initH[pick+1]
  dinitH <- seq(lowl,upl,length=len)
  zoneDepl <- numeric(len)
  numpop <- glob$numpop
  for (harv in 1:len) { # harv=1
    doharv <- rep(dinitH[harv],numpop)
    zoneDD <- runthreeH(zoneC=zoneC,zoneD,inHarv=doharv,glob)
    zoneDepl[harv] <- sum((zoneDD$matureB[1,]/sum(zoneDD$matureB[1,])) *
                            zoneDD$deplsB[1,])
  }
  pick <- which.closest(depl,zoneDepl)
  pickharv <- rep(dinitH[pick],numpop)
  zoneDD <- runthreeH(zoneC=zoneC,zoneD,inHarv=pickharv,glob)
  return(zoneDD)
} # end of Notdodepletion

























