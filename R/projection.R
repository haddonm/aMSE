

#' @title makezoneDR generates the container for the projection dynamics
#'
#' @description makezoneDR generates an object designed to hold the outputs
#'     from each replicate within a set of projections. This is identical to
#'     zoneD except it contains a repaet for each iteration.
#'
#' @param projyr the number of years of the projection
#' @param iter the number of replicates
#' @param glob the global object
#' @param inzoneD the zoneD object after any preliminary depletion
#'
#' @return a large list containing an object ready for the projection dynamics
#' @export
#'
#' @examples
#' print("wait for a considerd example")
makezoneDR <- function(projyr,iter,glob,inzoneD) {
  # projyr=projyrs; iter=reps; glob=glb; inzoneD=zoneDD
  numpop <- glob$numpop
  N <- glob$Nclass
  namedims <- list(1:projyr,1:numpop,1:iter)
  namendims <- list(glob$midpts,1:projyr,1:numpop,1:iter)
  SAU <- inzoneD$SAU #as.numeric(sapply(zoneC,"[[","SAU"))
  MatB <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  MatB[1,,] <- inzoneD$matureB[1,]
  ExplB <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  ExplB[1,,] <- inzoneD$exploitB[1,]
  Catch <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  Catch[1,,] <- inzoneD$catch[1,]
  Harvest <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  Harvest[1,,] <- inzoneD$harvestR[1,]
  cpue <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  cpue[1,,] <- inzoneD$cpue[1,]
  Recruit <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  Recruit[1,,] <- inzoneD$recruit[1,]
  deplSpB <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  deplSpB[1,,] <- inzoneD$deplsB[1,]
  deplExB <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  deplExB[1,,] <- inzoneD$depleB[1,]
  CatchN <- array(data=0,dim=c(N,projyr,numpop,iter),dimnames=namendims)
  CatchN[,1,,] <- inzoneD$catchN[,1,]
  Nt <- array(data=0,dim=c(N,projyr,numpop,iter),dimnames=namendims)
  Nt[,1,,] <- inzoneD$Nt[,1,]
  zoneDP <- list(SAU=SAU,matureB=MatB,exploitB=ExplB,catch=Catch,
                 harvestR=Harvest,cpue=cpue,recruit=Recruit,
                 deplsB=deplSpB,depleB=deplExB,catchN=CatchN,Nt=Nt)
  return(zoneDP)
} # end of makezoneDR
