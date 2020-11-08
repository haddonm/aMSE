
#' @title addrepvar adds variation to the start of a generic simulation run
#'
#' @description addrepvar conducts projyrs of simulation at the same constant
#'     harvest rate used to deplete each population to the desired level, only
#'     it does this with recruitment variability turned on, which sets up a
#'     series of replicates with different starting positions. The final step
#'     is to take the final 10 yearly values from each population and replicate
#'     and place them in the first ten years of the projection matrices within
#'     the zoneDR object
#'
#' @param zoneC the constant zone object. Its selectivity patterns should be
#'     altered to reflect the projection conditions using modzoneCSel
#' @param zoneDR the dynamic zone object expanded to include all replicates,
#'     using makezoneDR
#' @param inHt the initial harvest strategy that was used to deplete the zone
#' @param glb the globals object
#' @param ctrl the ctrl object
#' @param nyrs defaults to 25, which appears to be enough to produce the
#'     required variation
#' @param reset should the first ten rows of the each output be set to the last
#'     10 rows of the initial variation run. default = TRUE. If FALSE the
#'     reps trajectories when exposed to recruitment variation are passed out
#'
#' @return a revised zoneDR object
#' @export
#'
#' @examples
#' print("wait on suitable data")
addrepvar <- function(zoneC,zoneDR,inHt,glb,ctrl,nyrs=25,reset=TRUE) {
  # zoneC=zoneC;zoneDR=zoneDR;inHt=zoneDR$harvestR;glb=glb;ctrl=ctrl
  varzoneDR <- zoneDR
  sigmar <- ctrl$withsigR
  npop <- glb$numpop
  Ncl <- glb$Nclass
  movem <- glb$move
  reps <- ctrl$reps
  matb <- numeric(npop)
  for (iter in 1:reps) {
    for (year in 2:nyrs) {
      for (popn in 1:npop) { # year=2; iter=1; popn=1
        out <- oneyear(inpopC=zoneC[[popn]],inNt=zoneDR$Nt[,year-1,popn,iter],
                       Nclass=Ncl,inH=inHt[1,popn,iter],yr=year)
        zoneDR$exploitB[year,popn,iter] <- out$ExploitB
        zoneDR$matureB[year,popn,iter] <- out$MatureB
        zoneDR$catch[year,popn,iter] <- out$Catch
        zoneDR$harvestR[year,popn,iter] <- out$Harvest
        zoneDR$cpue[year,popn,iter] <- out$ce
        zoneDR$Nt[,year,popn,iter] <- out$Nt
        zoneDR$catchN[,year,popn,iter] <- out$CatchN
        matb[popn] <- out$MatureB
      } # pop
      steep <- getvect(zoneC,"steeph") #sapply(zoneC,"[[","popdef")["steeph",]
      r0 <- getvar(zoneC,"R0") #sapply(zoneC,"[[","R0")
      b0 <- getvar(zoneC,"B0") #sapply(zoneC,"[[","B0")
      recs <- oneyearrec(steep,r0,b0,matb,sigR=sigmar)
      newrecs <- movem %*% recs
      zoneDR$recruit[year,,iter] <- newrecs
      zoneDR$Nt[1,year,,iter] <- newrecs
      zoneDR$deplsB[year,,iter] <- zoneDR$matureB[year,,iter]/b0
      zoneDR$depleB[year,,iter] <- zoneDR$exploitB[year,,iter]/getvar(zoneC,"ExB0")
    }   # year loop        zoneDR$matureB[,,1]
  }     # rep loop
  if (reset) {
    varzoneDR$matureB[1:10,,] <- zoneDR$matureB[(nyrs-9):nyrs,,]
    varzoneDR$exploitB[1:10,,] <- zoneDR$exploitB[(nyrs-9):nyrs,,]
    varzoneDR$catch[1:10,,] <- zoneDR$catch[(nyrs-9):nyrs,,]
    varzoneDR$harvestR[1:10,,] <- zoneDR$harvestR[(nyrs-9):nyrs,,]
    varzoneDR$cpue[1:10,,] <- zoneDR$cpue[(nyrs-9):nyrs,,]
    varzoneDR$recruit[1:10,,] <- zoneDR$recruit[(nyrs-9):nyrs,,]
    varzoneDR$deplsB[1:10,,] <- zoneDR$deplsB[(nyrs-9):nyrs,,]
    varzoneDR$depleB[1:10,,] <- zoneDR$depleB[(nyrs-9):nyrs,,]
    varzoneDR$catchN[,1:10,,] <- zoneDR$catchN[,(nyrs-9):nyrs,,]
    varzoneDR$Nt[,1:10,,] <- zoneDR$Nt[,(nyrs-9):nyrs,,]
  } else {
    varzoneDR <- zoneDR
  }
  return(varzoneDR)
} # end of addrepvar

#' @title asSAU generates an output list of SAU values from the in put zoneDP
#'
#' @description asSAU takes the population based results in zoneDP and converts
#'     them into a list based around SAU
#'
#' @param zoneDP the population based results list zoneDP
#' @param sauindex the SAU index of each population
#' @param saunames the names of each SAU
#' @param b0 a vector of the B0 for each SAU
#' @param exb0 a vector of the ExB0 for each SAU
#'
#' @return a list of results based around SAU
#' @export
#'
#' @examples
#' print("wait on new data")
asSAU <- function(zoneDP,sauindex,saunames,b0,exb0) {
  #  zoneDP=zoneDP; sauindex=sauindex;saunames=saunames;b0=B0;exb0=exB0
  # dimension <- dim(zoneDP$matureB)
  # projyrs <- dimension[1]
  # reps <- dimension[3]
  # nsau <- length(saunames)
  matB <- sumpops(zoneDP$matureB,sauindex,saunames)
  expB <- sumpops(zoneDP$exploitB,sauindex,saunames)
  catS <- sumpops(zoneDP$catch,sauindex,saunames)
  recS <- sumpops(zoneDP$recruit,sauindex,saunames)
  harvS <- catS/expB
  saudeplsB <- calcsau(matB,saunames,b0)
  saudepleB <- calcsau(expB,saunames,exb0)
  ans <- list(matB=matB,expB=expB,catS=catS,harvS=harvS,recS=recS,
              saudeplsB=saudeplsB,saudepleB=saudepleB)
  return(ans)
} #end of asSAU



#' @title calcsau compares an input variable to a constant for each SAU
#'
#' @description calcsau divides an input variable, such as matureB or exploitB
#'     by their respective constant, such as B0 or ExB0, from a vector of
#'     values for each SAU. This is used to calculate the depletion levels of
#'     mature and exploitable biomass
#'
#' @param invar either zoneDP$matureB or zoneDP$exploitB
#' @param saunames the names of each SAU
#' @param ref0 either a vector of B0 or ExB0 for each SAU
#'
#' @return an array of projyrs x nSAU x reps
#' @export
#'
#' @examples
#' print("wait on new data")
calcsau <-  function(invar,saunames,ref0) {# for deplsb depleB
  dimension <- dim(invar)
  projyrs <- dimension[1]
  reps <- dimension[3]
  nsau <- length(saunames)
  res <- array(0,dim=c(projyrs,nsau,reps))
  dimnames(res) <- list(1:projyrs,saunames,1:reps)
  for (sau in 1:nsau) res[,sau,] <- invar[,sau,]/ref0[sau]
  return(res)
} # end of calcsau

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



#' @title sumpops takes the zoneDP and sums the populations within each SAU
#'
#' @description sumpops summarizes the matureB, exploitB, catch, recruit
#'     variables within the depleted zone that has undergone the replicate
#'     application of an HS, and sums the values from each population within
#'     each SAU.
#'
#' @param invar either zoneDP$matureB or exploitB, or catch, or recruit
#' @param sauindex the SAU index of each population
#' @param saunames the names of each SAU
#'
#' @return an array of projyrs x nSAU x reps
#' @export
#'
#' @examples
#' print("wait on new data")
sumpops <- function(invar,sauindex,saunames) {
  # for matureB, exploitB, catch, recruit
  dimension <- dim(invar)
  projyrs <- dimension[1]
  reps <- dimension[3]
  nsau <- length(saunames)
  res <- array(0,dim=c(projyrs,nsau,reps))
  dimnames(res) <- list(1:projyrs,saunames,1:reps)
  for (iter in 1:reps) {
    for (i in 1:nsau) { # iter=1; i = 1
      pickc <- which(sauindex == i)
      res[,i,iter] <- rowSums(invar[,pickc,iter])
    }
  }
  return(res)
} # end of sumpops

