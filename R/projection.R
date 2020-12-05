
#' @title addrecvar adds recruitment variation to start of generic simulation
#'
#' @description addrecvar conducts projyrs of simulation at the same constant
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
addrecvar <- function(zoneC,zoneDR,inHt,glb,ctrl,nyrs=25,reset=TRUE) {
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
} # end of addrecvar

#' @title asSAU generates an output list of SAU values from the in put zoneDP
#'
#' @description asSAU takes the population based results in zoneDP and converts
#'     them into a list based around SAU. saucpue is already calculated during
#'     the application of the harvest strategy.
#'
#' @param projzone the population based results list zoneDP plus the saucpue
#'     all from applymcda
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
asSAU <- function(projzone,sauindex,saunames,b0,exb0) {
# projzone=mseproj;sauindex=sauindex;saunames=saunames;B0=B0;exB0=exB0
  zoneDP <- projzone
  cpue <- zoneDP$cesau
  matB <- sumpops(zoneDP$matureB,sauindex,saunames)
  expB <- sumpops(zoneDP$exploitB,sauindex,saunames)
  catS <- sumpops(zoneDP$catch,sauindex,saunames)
  recS <- sumpops(zoneDP$recruit,sauindex,saunames)
  harvS <- catS/expB
  saudeplsB <- calcsau(matB,saunames,b0)
  saudepleB <- calcsau(expB,saunames,exb0)
  ans <- list(matB=matB,expB=expB,catS=catS,cpue=cpue,harvS=harvS,
              recS=recS,saudeplsB=saudeplsB,saudepleB=saudepleB)
  return(ans)
} #end of asSAU

#' @title aszone sums the various SAU dynamics and fishery properties
#'
#' @description aszone calculates zone wide totals by summing the mature and
#'     exploitable biomass across SAU, it also sums the catches and recruitment
#'     levels. It calculates zone wide harvest rates by dividing the catches by
#'     the available exploitable biomass. It also calculates the depletion
#'     levels of both the mature and exploitable biomass by dividing their
#'     zone totals through time by the initial zone wide B0 and ExB0. Finally,
#'     it calculates a zone-wide cpue by catch-weighting each of the SAU cpue
#'     values and then summing each set of SAU for each year and iteration.
#'
#' @param sauzone the output from the applyharvest strategy function, a large
#'     list of arrays or results
#' @param zoneCP the constant part of the projection zone
#'
#' @return a list of 8 nyrs x reps matrices, summarizing the fishery outputs at
#'     the geographical scale of the zone
#' @export
#'
#' @examples
#' print("wait on data")
aszone <- function(sauzone,zoneCP) { # sauzone=sauzoneDP; zoneCP=zoneCP; reps=ctrl$reps
  B0 <- sum(getvar(zoneCP,"B0"),na.rm=TRUE)
  exB0 <- sum(getvar(zoneCP,"ExB0"),na.rm=TRUE)
  zonesB <- apply(sauzone$matB,c(1,3),sum,na.rm=TRUE) # spawning biomass
  zoneeB <- apply(sauzone$expB,c(1,3),sum,na.rm=TRUE) # exploitable biomass
  zoneC <- apply(sauzone$catS,c(1,3),sum,na.rm=TRUE)  # catches
  zoneH <- zoneC/zoneeB                               # gross harvest rate
  zoneR <- apply(sauzone$recS,c(1,3),sum,na.rm=TRUE)  # total recruitment
  zonedeplsB <- zonesB/B0               # spawning biomass depletion
  zonedepleB <- zoneeB/exB0             # xploitable biomass depletion
  nyrs <- nrow(zonesB)
  reps <- ncol(zonesB)
  zonece <- matrix(0,nrow=nyrs,ncol=reps,dimnames=list(1:nyrs,1:reps))
  for (iter in 1:reps) { # generate SAU catch-weighted zonal cpue
    wts <- sauzone$catS[,,iter]/zoneC[,iter]
    zonece[,iter] <- rowSums((sauzone$cpue[,,iter] * wts),na.rm=TRUE)
  }
  outzone=list(zonesB=zonesB,zoneeB=zoneeB,zoneC=zoneC,zoneH=zoneH,zoneR=zoneR,
               zonece=zonece,zonedeplsB=zonedeplsB,zonedepleB=zonedepleB)
  return(outzone)
} # end of aszone


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
#'     zoneD except it contains a repeat for each iteration and now includes
#'     a cesau, which is the population-catch-weighted cpue for each SAU.
#'     Variation should be added by addrecvar.
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
#' print("Could add variation to the harvest rates so that when ")
#' print("prepareprojection was run the range of initial H values would increase ")
makezoneDR <- function(projyr,iter,glob,inzoneD) {
  # projyr=projyrs; iter=reps; glob=glb; inzoneD=zoneDD
  numpop <- glob$numpop
  nSAU <- glob$nSAU
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
  sauindex <- glob$sauindex
  cesau <- array(0,dim=c(projyr,nSAU,iter),
                 dimnames=list(1:projyr,1:nSAU,1:iter))
  catwt <- tapply(inzoneD$catch[1,],sauindex,sum)
  wts <- inzoneD$catch[1,] / catwt[sauindex]
  for (i in 1:nSAU)
    cesau[1,i,] <- sum(inzoneD$cpue[1,(sauindex==i)] * wts[(sauindex==i)])
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
                 harvestR=Harvest,cpue=cpue,cesau=cesau,recruit=Recruit,
                 deplsB=deplSpB,depleB=deplExB,catchN=CatchN,Nt=Nt)
  return(zoneDP)
} # end of makezoneDR

#' @title modprojC produces three objects used to condition the zone
#'
#' @description modprojC produces an object used to condition the zone when
#'     projecting it following any conditioning. Once the initial conditions
#'     for the projection have been attained through an initial depletion of an
#'     unfished equilibrium, or more particularly by conditioning on historical
#'     data then the projC object will be used to condition the projections. In
#'     essence, this puts the selectivity and selectivity x weight-at-size into
#'     projC
#'
#' @param zoneC the constant part of the zone
#' @param glob the globals for the simulation
#' @param inzone the zone object from readctrlzone
#'
#' @return a list of projSel and projSelWt within the projC object from inzone
#' @export
#'
#' @examples
#' print("wait on new example data")
modprojC <- function(zoneC,glob,inzone) { # zoneC=zoneC; glob=glb; inzone=zone1
  numpop <- glob$numpop
  popdefs <- getlistvar(zoneC,"popdef")
  sau <- getvar(zoneC,"SAU")
  midpts <- glob$midpts
  projSel <- matrix(0,nrow=glob$Nclass,ncol=numpop,dimnames=list(midpts,sau))
  projSelWt <- matrix(0,nrow=glob$Nclass,ncol=numpop,dimnames=list(midpts,sau))
  pLML <- inzone$projC$projLML[1] # needs development to allow variation in projLML
  for (pop in 1:numpop) {
    selL50 <- popdefs["SelP1",pop]
    selL95 <- popdefs["SelP2",pop]
    projSel[,pop] <- logistic((pLML + selL50),selL95,midpts)
    projSelWt[,pop] <- projSel[,pop] * zoneC[[pop]]$WtL
  }
  projC <- inzone$projC
  projC$Sel <- projSel
  projC$SelWt <- projSelWt
  return(projC=projC)
} # end of modprojC

#' @title modzoneCSel changes the selectivity characteristics in zoneC
#'
#' @description modzoneCSel  changes the selectivity characteristics in zoneC
#'     which is necessary when making a projection under a different LML.
#'     Currently this function is designed to allow only a single LML during any
#'     projections. It will obviously need modification if it is desired to
#'     explore the option of gradually changing an LML
#'
#' @param zoneC the constant zone object from setupzone
#' @param sel the new selectivity as a matrix of selectivity by population
#' @param selwt new selectivity x WtL as a matrix of SelWt x Population
#' @param glb the globals object
#' @param yrs the number of years of projection
#'
#' @return the zoneC object modified ready for projection
#' @export
#'
#' @examples
#' print("wait on suitable data")
modzoneCSel <- function(zoneC,sel,selwt,glb,yrs) {
  numpop <- glb$numpop
  Nclass <- glb$Nclass
  for (pop in 1:numpop){ # pop = 1
    select <- matrix(sel[,pop],nrow=Nclass,ncol=yrs)
    selectwt <- matrix(selwt[,pop],nrow=glb$Nclass,ncol=yrs)
    zoneC[[pop]]$Select <- select
    zoneC[[pop]]$SelWt <- selectwt
  }
  return(zoneC)
} # end of modzoneCSel

#' @title prepareprojection high level function that sets up a projection
#'
#' @description prepareprojection is a high level function that restructures
#'     zoneC by including the projected selectivity (in case the LML changes),
#'     and also selectivity x weight-at-size (for computational speed). It
#'     then generates a new zoneD with room for all replicates. Finally, it
#'     adds recruitment variabilkity to each of the replicates and keeps
#'     zone1$projC$inityrs of each replicates series as the start of each
#'     replicate projection.
#'
#' @param zone1 the original zone1 object from readctrlzone
#' @param zoneC the constant part of the zone
#' @param glb the global variables
#' @param zoneDep the zone after initial depletion
#' @param ctrl the ctrl object
#'
#' @return a list of the dynamic zone object as a list of arrays of projyrs x
#'     populations x replicates, plus the revised projC and revised zoneC
#' @export
#'
#' @examples
#' print("wait on data files")
prepareprojection <- function(zone1,zoneC,glb,zoneDep,ctrl) {
# zone1=zone$zone1;zoneC=zone$zoneC; glb=glb; zoneDep=zoneDD; ctrl=ctrl
  projyrs <- zone1$projC$projyrs + zone1$projC$inityrs
  projC <- modprojC(zoneC,glb,zone1) # include selectivity into projC
  zoneC <- modzoneCSel(zoneC,projC$Sel,projC$SelWt,glb,projyrs)
  zoneDR <- makezoneDR(projyrs,ctrl$reps,glb,zoneDep) # zoneDReplicates
  zoneDRp <- addrecvar(zoneC,zoneDR,zoneDR$harvestR,glb,ctrl)
  return(list(zoneDP=zoneDRp,projC=projC,zoneC=zoneC))
} # end of prepareprojection

#' @title scaleto1 scales an input vector of CPUE to a mean of one x avCE
#'
#' @description scaleto1 scales a vector of CPUE to a mean of
#'     one or avCE. The use of a mean of one means that visual comparisons
#'     between different time-series becomes visually simplified. The
#'     avCE option could be used to scale the CPUE to the average
#'     geometric mean - so as to put it on the nominal scale
#'
#' @param invect a vector of linear scale CPUE
#'
#' @return a vector of CPUE re-scaled to a mean of one
#' @export
#'
#' @examples
#'  ce <- c(0.4667187,1.2628564,0.8442146,0.9813531, 0.5554076,0.7426321)
#'  scaleto1(ce)
scaleto1 <- function(invect) {
  return(invect/mean(invect,na.rm=TRUE))
} # end of scaleto1

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
sumpops <- function(invar,sauindex,saunames) { #invar=zoneDP$matureB
  # for matureB, exploitB, catch, recruit
  dimension <- dim(invar)
  projyrs <- dimension[1]
  reps <- dimension[3]
  nsau <- length(saunames)
  res <- array(0,dim=c(projyrs,nsau,reps))
  dimnames(res) <- list(1:projyrs,saunames,1:reps)
  for (iter in 1:reps) {# iter=1; i = 1
    for (i in 1:nsau) {
      res[,i,iter] <- rowSums(invar[,(sauindex == i),iter])
    }
  }
  return(res)
} # end of sumpops

