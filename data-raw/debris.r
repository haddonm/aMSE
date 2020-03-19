


#' @title makeZone the list of numpop populations defines the abalone zone
#'
#' @description makeZone the list of numpop populations defines the abalone
#'     zone. This is an S3 class names 'zone'
#'
#' @param condDat the matrix of biological constants read from the data file
#' @param uplim the upper limit on the harvest rate used when estimating
#'     the productivity of each population, defaults to 0.4; a scalar
#'
#' @return a list containing the simulated zone, the production matrix from
#'     the function doproduction, and the popdefs used to create the zone.
#' @export
#'
#' @examples
#' \dontrun{
#' data(condDat)  # replaces the lines above
#' out <- makeZone(condDat,uplim=0.38)  # Define the Zone
#' zone <- out$zone
#' str(zone[[1]],max.level=1)
#' production <- out$Production
#' popdefs <- out$popdefs
#' print(round(t(popdefs),3))
#' }
makeZone <- function(condDat,uplim=0.4) { # condDat=condDat; uplim=0.4
  glb <- condDat$globals
  for (i in 1:length(glb))
    assign(names(condDat$globals)[i],condDat$globals[[i]]) #  defines:

  nblock <- glb$nblock
  blkdef <- condDat$blkpop
  numpop <- condDat$globals$numpop
  blockI <- defineBlock(nblock,blkdef,numpop)
  blockNames <- condDat$blockNames
  midpts <- condDat$globals$midpts
  #   Nyrs <- condDat$outYear[1]
  #   glb <- condDat$globals
  set.seed(condDat$randomseed)
  projectionLML <- condDat$projLML
  historicalLML <- condDat$histLML
  if (condDat$Condition)
    projLML <- historicalLML else projLML <- projectionLML
  zone <- vector("list",(numpop))
  imph <- seq(0.01,uplim,0.01)   # harvest rate for productivity
  numrow <- length(imph)
  pops <- seq(1,numpop,1)
  columns <- c("ExB","MatB","AnnH","Catch","Deplet","RelCE")
  production <- array(0,dim=c(numpop,numrow,6),
                      dimnames=list(pops,imph,columns))
  popdefs <- definepops(nblock,blockI,condDat) # define all pops in one go
  for (pop in 1:numpop) {      # pop <- 1
    popdef <- popdefs[pop,]
    zone[[pop]] <- makeabpop(popdef,midpts,projLML)
    #  out <- doproduction(zone[[pop]],uplim=uplim)
    #  zone[[pop]]$MSY <- out$MSY
    #  zone[[pop]]$MSYDepl <- out$MSYDepl
    #   production[pop,,] <- out$Productivity
    tmpL <- oneyrgrowth(zone[[pop]],zone[[pop]]$SaM)
    zone[[pop]]$bLML <- oneyrgrowth(zone[[pop]],tmpL)
    zone[[pop]]$blockName <- blockNames[blockI[pop]]
  }
  class(zone) <- "zone"
  # ans <- list(zone,production,popdefs)
  # names(ans) <- c("zone","production","popdefs")
  ans <- list(zone,popdefs)
  names(ans) <- c("zone","popdefs")
  return(ans)
}  # End of makeZone
