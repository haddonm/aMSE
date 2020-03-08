
library(rutilsMH)
library(aMSE)
library(microbenchmark)
datadir <- "./../../rcode2/aMSE/data-raw/"

  data(condDat)

#makeZone <- function(condDat,uplim=0.4) { #
#   condDat=condDat;  uplim=0.4
  glb <- condDat$globals
  for (i in 1:length(glb)) assign(names(glb)[i],glb[[i]])
  ls()

  nblock <- glb$nblock
  blkdef <- condDat$blkpop
  numpop <- glb$numpop
  blockI <- defineBlock(nblock,blkdef,numpop)
  blockNames <- condDat$blockNames
  midpts <- glb$midpts
  #   Nyrs <- condDat$outYear[1]
  #   glb <- condDat$globals
  set.seed(condDat$randomseed)
  projectionLML <- condDat$projLML
  historicalLML <- condDat$histLML
  if (condDat$Condition)
    projLML <- historicalLML else projLML <- projectionLML
  zoneD <- vector(list,numpop)  # isnt it just a single list
  zone <- vector("list",numpop)
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

