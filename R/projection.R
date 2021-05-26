
#' @title addrecvar adds recruitment variation to end of conditioning step
#'
#' @description addrecvar replicates the historical depletion step in zoneDD
#'     reps times and copies the first Nyrs - varyrs of values into each
#'     replicate. Then it steps through each replicate adding recruitment
#'     variation to the dynamics for the last varyrs years. This is so that
#'     when the projections begin each replicate will already have recruitment
#'     variability fully developed and there will be no time lag on the
#'     recruitment and subsequent dynamics.
#'
#' @param zoneDD the dynamic zone after historical fishery data has been used
#'     to finalize the conditioning
#' @param zoneDP the empty dynamic zone object ready to receive its first year
#' @param zoneC the zone's constants object for each population
#' @param glob the object containing the global constants
#' @param condC the object containing the historical fishery data
#' @param ctrl the control object
#' @param varyrs the number of years at the end of the historical period to
#'     which recruitment variation is to be added
#' @param multTAC the TAC multiplication matrix from the HCR
#' @param calcpopC a function that takes the output from hcrfun and gernerates
#'     the actual catch per population expected in the current year.
#' @param sigR the initial recruitment variation default=1e-08
#' @param sigB the initial biomass cpuie variation default = 1e-08
#' @param lastsigR the recruitment variation to be added to the final varyrs
#'
#' @return an initialized dynamics zone object for the projecitons with the
#'     first year populated
#' @export
#'
#' @examples
#' print("wait on suitable data-sets")
addrecvar <- function(zoneDD,zoneDP,zoneC,glob,condC,ctrl,varyrs,multTAC,
                      calcpopC,sigR=1e-08,sigB=1e-08,lastsigR=0.3) {
 #  zoneDD=zoneDD;zoneDP=zoneDR;zoneC=zoneC;glob=glb
 #  condC=condC;ctrl=ctrl;varyrs=6;lastsigR=lastsigR; multTAC=multTAC
 #  sigR=1e-08; sigB=1e-08; lastsigR=0.1
  sauindex <- glob$sauindex
  histC <- condC$histCatch
  yrs <- condC$histyr[,"year"]
  nyrs <- length(yrs)
  reps <- ctrl$reps
  zoneDDR <- makezoneDP(nyrs,reps,glob)
  finalyr <- nyrs - varyrs
  r0 <- getvar(zoneC,"R0") #sapply(zoneC,"[[","R0")
  b0 <- getvar(zoneC,"B0") #sapply(zoneC,"[[","B0")
  exb0 <- getvar(zoneC,"ExB0")
  zoneDDR$matureB[1:finalyr,,] <- zoneDD$matureB[1:finalyr,]
  zoneDDR$exploitB[1:finalyr,,] <- zoneDD$exploitB[1:finalyr,]
  zoneDDR$midyexpB[1:finalyr,,] <- zoneDD$midyexpB[1:finalyr,]
  zoneDDR$catch[1:finalyr,,] <- zoneDD$catch[1:finalyr,]
  zoneDDR$harvestR[1:finalyr,,] <- zoneDD$harvestR[1:finalyr,]
  zoneDDR$cpue[1:finalyr,,] <- zoneDD$cpue[1:finalyr,]
  zoneDDR$recruit[1:finalyr,,] <- zoneDD$recruit[1:finalyr,]
  zoneDDR$catchN[,1:finalyr,,] <- zoneDD$catchN[,1:finalyr,]
  zoneDDR$Nt[,1:finalyr,,] <- zoneDD$Nt[,1:finalyr,]
  for (iter in 1:reps) { # iter = 1
    for (year in (finalyr+1):nyrs) { # year = finalyr + 1
      catchsau <- histC[year,]
      hcrout <- list(acatch=catchsau)
      popC <- calcpopC(hcrout,exb=zoneDD$exploitB[year-1,],
                             sauindex,sigmab=sigB)
    #  exb <- zoneDDR$exploitB[year-1,,iter]
      inN <- zoneDDR$Nt[,year-1,,iter]
      out <- oneyearsauC(zoneCC=zoneC,inN=inN,popC=popC,year=year,
                         Ncl=glob$Nclass,sauindex=sauindex,movem=glob$move,
                         sigmar=lastsigR,r0=r0,b0=b0,exb0=exb0)
      dyn <- out$dyn
      saudyn <- poptosauCE(dyn["catch",],dyn["cpue",],sauindex)
      zoneDDR$exploitB[year,,iter] <- dyn["exploitb",]
      zoneDDR$midyexpB[year,,iter] <- dyn["midyexpB",]
      zoneDDR$matureB[year,,iter] <- dyn["matureb",]
      zoneDDR$catch[year,,iter] <- dyn["catch",]
      zoneDDR$harvestR[year,,iter] <- dyn["catch",]/out$dyn["exploitb",]
      zoneDDR$cpue[year,,iter] <- dyn["cpue",]
      zoneDDR$cesau[year,,iter] <- saudyn$saucpue
      zoneDDR$catsau[year,,iter] <- saudyn$saucatch
      zoneDDR$recruit[year,,iter] <- dyn["recruits",]
      zoneDDR$deplsB[year,,iter] <- dyn["deplsB",]
      zoneDDR$depleB[year,,iter] <- dyn["depleB",]
      zoneDDR$Nt[,year,,iter] <- out$NaL
      zoneDDR$catchN[,year,,iter] <- out$catchN
    }
  }
  endcatch <- tapply(zoneDDR$catch[nyrs,,1],sauindex,sum,na.rm=TRUE)
  yrce <- nrow(multTAC)
  acatch <- endcatch * multTAC[yrce,]  # predicted aspirational catches
  sigmar=ctrl$withsigR
  sigmab=ctrl$withsigB
  for (iter in 1:reps) {  #  iter=1
  #  exb=zoneDDR$exploitB[nyrs,,iter]
    inN=zoneDDR$Nt[,nyrs,,iter]
    hcrout <- list(acatch=catchsau)
    popC <- calcpopC(hcrout,exb=zoneDD$exploitB[year-1,],
                           sauindex,sigmab=sigB)
    outy <- oneyearsauC(zoneCC=zoneC,inN=inN,popC=popC,year=1,
                        Ncl=glob$Nclass,sauindex=sauindex,movem=glob$move,
                        sigmar=sigmar,r0=r0,b0=b0,exb0=exb0)
    dyn <- outy$dyn
    saudyn <- poptosauCE(dyn["catch",],dyn["cpue",],sauindex)
    zoneDP$exploitB[1,,iter] <- dyn["exploitb",]
    zoneDP$midyexpB[1,,iter] <- dyn["midyexpB",]
    zoneDP$matureB[1,,iter] <- dyn["matureb",]
    zoneDP$catch[1,,iter] <- dyn["catch",]
    zoneDP$acatch[1,,iter] <- acatch
    zoneDP$catsau[1,,iter] <- saudyn$saucatch
    zoneDP$harvestR[1,,iter] <- dyn["catch",]/dyn["exploitb",]
    zoneDP$cpue[1,,iter] <- dyn["cpue",]
    zoneDP$cesau[1,,iter] <- saudyn$saucpue
    zoneDP$recruit[1,,iter] <- dyn["recruits",]
    zoneDP$deplsB[1,,iter] <- dyn["deplsB",]
    zoneDP$depleB[1,,iter] <- dyn["depleB",]
    zoneDP$Nt[,1,,iter] <- outy$NaL
    zoneDP$catchN[,1,,iter] <- outy$catchN
  }
  return(list(zoneDP=zoneDP,zoneDDR=zoneDDR))
} # end of addrecvar

#' @title asSAU generates an output list of SAU values from the in put zoneDP
#'
#' @description asSAU takes the population based results in zoneDP and converts
#'     them into a list based around SAU. saucpue is already calculated during
#'     the application of the harvest strategy.
#'
#' @param projzone the population based results list zoneDP plus the saucpue
#'     all from doprojection
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

#' @title calcprojsel generates the selectivity and selectivity x weight
#'
#' @description calcprojsel is used to produce each projection years'
#'     selectivity and selectivity x weight-at-size. projSel is an
#'     Nclass x projyrs matrix, a separate selectivity for each projection year,
#'     which allows for changing the selectivity during the projection. The
#'     projSelWt is a Nclass x projyrs x numpop array, required because each
#'     population's weight-at-size relationship will be slihtly different.
#'
#' @param zoneC the constants object
#' @param projC the projection object from zone1
#' @param glb the globals object
#'
#' @return a list of projSel and projSelWt
#' @export
#'
#' @examples
#'  print("wait on new example data")
calcprojsel <- function(zoneC,projC,glb) {
  # selP=c(,zoneC[[1]]$popdef["SelP2"]);
  # glb=glb; projC=zone1$projC
  selL50 <- zoneC[[1]]$popdef["SelP1"]
  selL95 <- zoneC[[1]]$popdef["SelP2"]
  projyrs <- projC$projyrs
  numpop <- glb$numpop
  pLML <- projC$projLML
  midpts <- glb$midpts
  Nclass <- glb$Nclass
  projSel <- matrix(0,nrow=Nclass,ncol=projyrs,dimnames=list(midpts,1:projyrs))
  projSelWt <- array(0,dim=c(Nclass,projyrs,numpop),
                     dimnames=list(midpts,1:projyrs,1:numpop))
  for (yr in 1:projyrs) {
    projSel[,yr] <- logistic((pLML[yr] + selL50),selL95,midpts)
    for (pop in 1:numpop) {
      projSelWt[,yr,pop] <- projSel[,yr] * zoneC[[pop]]$WtL
    }
  }
  return(list(projSel=projSel,projSelWt=projSelWt))
} # end of calcprojsel


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



#' @title doprojections conducts the replicate model runs for Tasmania
#'
#' @description doprojections conducts the replicate model runs for
#'     Tasmania using the mcdahcr and the input hsargs
#'
#' @param ctrl the ctrl object from readctrlfile
#' @param zoneDP the object used to contain the dynamics from the replicate
#'     model runs
#' @param zoneCP the object used to contain the constants for each population
#'     used in model dynamics
#' @param otherdata an object containing any other data used by the different
#'     functions used n the projections. For example, for Tasmania this would
#'     be pointed at the histCE data set which is appended to each year's new
#'     predicted cpue data from each replicate to calibrate the HCR within the
#'     current HS. Other jurisdictions may use other data.
#' @param glb the object containing the global constants for the given run
#' @param hcrfun the name of the harvest control rule that is used to
#'     calculate the multiplier for the previous aspirational catches (possibly
#'     for each SAU but possibly the TAC for the wholse zone) so as to
#'     estimate the aspirational catches or TAC or the following year
#' @param hsargs the constants used to define the workings of the hcr
#' @param sampleCE a function that generates the CPUE statistics
#' @param sampleFIS a function that generates the FIS statistics
#' @param sampleNaS a function that generates the Numbers-at-size samples
#' @param getdata a function that gathers all the data required by the hcrfun
#'     and combines it into an hcrdata object ready for the hcrfun
#' @param calcpopC a function that takes the output from hcrfun and gernerates
#'     the actual catch per population expected in the current year.
#' @param ... the ellipsis used in case any of the functions hcrfun, sampleCE,
#'     sampleFIS, sampleNas, and getdata require extra arguments not included
#'     in the default named collection
#'
#' @return a replacement for zoneDP containing the dynamics for all populations,
#'     for all replicates, and for all projection years
#' @export
#'
#' @examples
#' print("wait on suitable internal data sets")
doprojections <- function(ctrl,zoneDP,zoneCP,otherdata,glb,hcrfun,hsargs,
                             sampleCE,sampleFIS,sampleNaS,getdata,calcpopC,...) {
#  ctrl=ctrl;zoneDP=zoneDP;zoneCP=zoneCP;otherdata=condC$histCE;glb=glb; getdata=tasdata
#  hcrfun=mcdahcr; hsargs=hsargs; sampleCE=tasCPUE;sampleFIS=tasFIS; sampleNaS=tasNaS
#  calcpopC=calcexpectpopC
  reps <- ctrl$reps
  projyrs <- ctrl$projection
  sigmar <- ctrl$withsigR
  sigmab <- ctrl$withsigB
  sigce <- ctrl$withsigCE
  oldyrs <- as.numeric(rownames(otherdata))
  sauindex <- glb$sauindex
  saunames <- glb$saunames
  Nclass <- glb$Nclass
  movem <- glb$move
  r0 <- getvar(zoneCP,"R0") #sapply(zoneC,"[[","R0")
  b0 <- getvar(zoneCP,"B0") #sapply(zoneC,"[[","B0")
  exb0 <- getvar(zoneCP,"ExB0")
  for (iter in 1:reps) {
    for (year in 2:projyrs) { # iter=1; year=2
      hcrdata <- getdata(sampleCE,sampleFIS,sampleNaS,zoneDP,otherdata,
                         year=year,iter=iter)
      hcrout <- hcrfun(hcrdata,hsargs,saunames=glb$saunames)
      popC <- calcpopC(hcrout,exb=zoneDP$exploitB[year-1,,iter],
                             sauindex,sigmab=sigmab)
      outy <- oneyearsauC(zoneCC=zoneCP,inN=zoneDP$Nt[,year-1,,iter],
                          popC=popC,year=year,Ncl=Nclass,sauindex=sauindex,
                          movem=movem,sigmar=sigmar,
                          r0=r0,b0=b0,exb0=exb0)
      dyn <- outy$dyn
      saudyn <- poptosauCE(dyn["catch",],dyn["cpue",],sauindex)
      zoneDP$exploitB[year,,iter] <- dyn["exploitb",]
      zoneDP$midyexpB[year,,iter] <- dyn["midyexpB",]
      zoneDP$matureB[year,,iter] <- dyn["matureb",]
      zoneDP$catch[year,,iter] <- dyn["catch",]
      zoneDP$acatch[year,,iter] <- hcrout$acatch
      zoneDP$catsau[year,,iter] <- saudyn$saucatch
      zoneDP$harvestR[year,,iter] <- dyn["catch",]/dyn["exploitb",]
      zoneDP$cpue[year,,iter] <- dyn["cpue",]
      zoneDP$cesau[year,,iter] <- saudyn$saucpue
      zoneDP$recruit[year,,iter] <- dyn["recruits",]
      zoneDP$deplsB[year,,iter] <- dyn["deplsB",]
      zoneDP$depleB[year,,iter] <- dyn["depleB",]
      zoneDP$Nt[,year,,iter] <- outy$NaL
      zoneDP$catchN[,year,,iter] <- outy$catchN
      zoneDP$NumNe[,year,,iter] <- outy$NumNe
      zoneDP$TAC[year,iter] <- hcrout$TAC
    } # year loop
  }   # iter loop
  return(zoneDP)
} # end of doprojections

#' @title makezoneDP generates the container for the projection dynamics
#'
#' @description makezoneDP generates an object designed to hold the outputs
#'     from each replicate within a set of projections. This is identical to
#'     zoneD except it contains a repeat for each iteration and now includes
#'     a cesau, which is the population-catch-weighted cpue for each SAU.
#'
#' @param projyr the number of years of the projection
#' @param iter the number of replicates
#' @param glb the global object
#'
#' @return a large list containing an object ready for the projection dynamics
#' @export
#'
#' @examples
#' print("Could add variation to the harvest rates so that when ")
#' print("prepareprojection was run the range of initial H values would increase ")
makezoneDP <- function(projyr,iter,glb) {
  numpop <- glb$numpop
  nSAU <- glb$nSAU
  saunames <- glb$saunames
  N <- glb$Nclass
  SAU <- glb$SAUpop
  namedims <- list(1:projyr,1:numpop,1:iter)
  namendims <- list(glb$midpts,1:projyr,1:numpop,1:iter)
  MatB <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  ExplB <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  midyexpB <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  catch <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  acatch <- array(0,dim=c(projyr,nSAU,iter),
                  dimnames=list(1:projyr,saunames,1:iter)) #aspirational catches
  harvest <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  cpue <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  cesau <- array(0,dim=c(projyr,nSAU,iter),
                 dimnames=list(1:projyr,saunames,1:iter))
  catsau <- array(0,dim=c(projyr,nSAU,iter),
                  dimnames=list(1:projyr,saunames,1:iter))
  recruit <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  deplSpB <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  deplExB <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  catchN <- array(data=0,dim=c(N,projyr,numpop,iter),dimnames=namendims)
  Nt <- array(data=0,dim=c(N,projyr,numpop,iter),dimnames=namendims)
  NumNe <- array(data=0,dim=c(N,projyr,numpop,iter),dimnames=namendims)
  TAC <- array(0,dim=c(projyr,iter),dimnames=list(1:projyr,1:iter))
  zoneDP <- list(SAU=SAU,matureB=MatB,exploitB=ExplB,midyexpB=midyexpB,catch=catch,
                 acatch=acatch,harvestR=harvest,cpue=cpue,cesau=cesau,
                 catsau=catsau,recruit=recruit,deplsB=deplSpB,depleB=deplExB,
                 catchN=catchN,Nt=Nt,NumNe=NumNe,TAC=TAC)
  return(zoneDP)
} # end of makezoneDP

#' @title makezoneDR generates the container for the projection dynamics
#'
#' @description makezoneDR generates an object designed to hold the outputs
#'     from each replicate within a set of projections. This is identical to
#'     zoneD except it contains a repeat for each iteration and now includes
#'     a cesau and a acatch, which is the population-catch-weighted cpue for
#'     each SAU, and the aspirational catch. Variation should be added by
#'     addrecvar.
#'
#' @param projyr the number of years of the projection
#' @param iter the number of replicates
#' @param glb the global object
#' @param inzoneD the zoneD object after any preliminary depletion
#'
#' @return a large list containing an object ready for the projection dynamics
#' @export
#'
#' @examples
#' print("Could add variation to the harvest rates so that when the ")
#' print("prepareprojection was run the range of initial H values would increase ")
makezoneDR <- function(projyr,iter,glb,inzoneD) {
  # projyr=projyrs; iter=reps; glob=glb; inzoneD=zoneDD
  numpop <- glb$numpop
  nSAU <- glb$nSAU
  N <- glb$Nclass
  sauindex <- glb$sauindex
  endyr <- nrow(inzoneD$catch)
  namedims <- list(1:projyr,1:numpop,1:iter)
  namendims <- list(glb$midpts,1:projyr,1:numpop,1:iter)
  SAU <- inzoneD$SAU #as.numeric(sapply(zoneC,"[[","SAU"))
  MatB <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  MatB[1,,] <- inzoneD$matureB[endyr,]
  ExplB <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  ExplB[1,,] <- inzoneD$exploitB[endyr,]
  Catch <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  Catch[1,,] <- inzoneD$catch[endyr,]
  aCatch <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims) #aspirational
  Harvest <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  Harvest[1,,] <- inzoneD$harvestR[endyr,]
  cpue <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  cesau <- array(0,dim=c(projyr,nSAU,iter),
                 dimnames=list(1:projyr,1:nSAU,1:iter))
  catsau <- array(0,dim=c(projyr,nSAU,iter),
                             dimnames=list(1:projyr,1:nSAU,1:iter))
  catwt <- tapply(inzoneD$catch[endyr,],sauindex,sum,na.rm=TRUE)
  catsau[1,,] <- catwt
  wts <- inzoneD$catch[endyr,] / catwt[sauindex]
  for (i in 1:nSAU)
     cesau[1,i,] <- sum(inzoneD$cpue[endyr,(sauindex==i)] * wts[(sauindex==i)])
  Recruit <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  Recruit[1,,] <- inzoneD$recruit[endyr,]
  deplSpB <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  deplSpB[1,,] <- inzoneD$deplsB[endyr,]
  deplExB <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  deplExB[1,,] <- inzoneD$depleB[endyr,]
  CatchN <- array(data=0,dim=c(N,projyr,numpop,iter),dimnames=namendims)
  CatchN[,1,,] <- inzoneD$catchN[,endyr,]
  Nt <- array(data=0,dim=c(N,projyr,numpop,iter),dimnames=namendims)
  Nt[,1,,] <- inzoneD$Nt[,endyr,]
  zoneDP <- list(SAU=SAU,matureB=MatB,exploitB=ExplB,catch=Catch,acatch=aCatch,
                 catsau=catsau,harvestR=Harvest,cpue=cpue,cesau=cesau,
                 recruit=Recruit,deplsB=deplSpB,depleB=deplExB,catchN=CatchN,
                 Nt=Nt)
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
#'     projC. It can handle changes in LML through the years of projection.
#'
#' @param zoneC the constant part of the zone
#' @param glob the globals for the simulation
#' @param projC the project definition object from zone1 from readctrlfile
#'
#' @return an updated projC with projSel and projSelWt completed
#' @export
#'
#' @examples
#' print("wait on new example data")
modprojC <- function(zoneC,glob,projC) { # zoneC=zone$zoneC; glob=glb; projC=zone1$projC
  numpop <- glob$numpop
  midpts <- glob$midpts
  projyrs <- projC$projyrs
  popdefs <- getlistvar(zoneC,"popdef")
  projSel <- array(0,dim=c(glob$Nclass,numpop,projyrs),
                   dimnames=list(midpts,glob$SAUpop,1:projyrs))
  projSelWt <- array(0,dim=c(glob$Nclass,numpop,projyrs),
                   dimnames=list(midpts,glob$SAUpop,1:projyrs))
  pLML <- projC$projLML
  diffLML <- unique(pLML)
  nLML <- length(diffLML)
  for (iLML in 1:nLML) { # iLML = 1
    pickyr <- which(pLML %in% diffLML[iLML])
    nyr <- length(pickyr)
    for (yr in 1:nyr) {
      for (pop in 1:numpop) { #  yr=1; pop=1
        selL50 <- popdefs["SelP1",pop]
        selL95 <- popdefs["SelP2",pop]
        projSel[,pop,pickyr[yr]] <- logistic((pLML[pickyr[yr]] + selL50),
                                             selL95,midpts)
        projSelWt[,pop,pickyr[yr]] <- projSel[,pop,pickyr[yr]] * zoneC[[pop]]$WtL
      }
    }
  }
  projC$Sel <- projSel
  projC$SelWt <- projSelWt
  return(projC=projC)
}

#' @title modzoneCSel changes the selectivity characteristics in zoneC
#'
#' @description modzoneCSel  changes the selectivity characteristics in zoneC
#'     which is necessary when making a projection under a different LML.
#'     This function is designed to allow multiple LML during any
#'     projections.
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
    zoneC[[pop]]$Select <- sel[,pop,]
    zoneC[[pop]]$SelWt <- selwt[,pop,]
  }
  return(zoneC)
} # end of modzoneCSel


#' @title prepareprojection high level function that sets up a projection
#'
#' @description prepareprojection is a high level function that restructures
#'     zoneC by including the projected selectivity (in case the LML is set to
#'     change during the projections), and also selectivity x weight-at-size
#'     (for computational speed). It then generates a new zoneD with room for
#'     all replicates as well as the aspirational catches from each HS. It does
#'     this by converting the arrays of year x pop, to year x pop x replicate.
#'     It then uses the conditioned data in zoneDep to predict the first
#'     aspirational catches for the projections and conducts the initial
#'     replicate, thus starting the application of the HS.
#'     Finally, it adds recruitment variation
#'     to each of the replicates and keeps the last year of each iteration of
#'     the addrecvar function as the start of each replicate projection, with
#'     zeros in the catch, cpue, and cesau arrays
#'
#' @param projC the projection object from readctrlfile
#' @param condC historical conditioning data
#' @param zoneC the constant part of the zone
#' @param glb the global variables
#' @param zoneDD the zone after initial depletion through conditioning on the
#'     fishery
#' @param ctrl the ctrl object for the scenario run
#' @param varyrs how many years at the end to add recruitment variation
#' @param multTAC the TAC multiplication matrix from the HCR
#' @param calcpopC a function that takes the output from hcrfun and gernerates
#'     the actual catch per population expected in the current year.
#' @param lastsigR recruitment variation for when it is applied for varyrs
#' @return a list of the dynamic zone object as a list of arrays of projyrs x
#'     populations x replicates, plus the revised projC and revised zoneC
#' @export
#'
#' @examples
#' print("wait on data files")
prepareprojection <- function(projC,condC,zoneC,glb,zoneDD,ctrl,varyrs,
                                 multTAC,calcpopC,lastsigR = 0.3) {
  # projC=projC;condC=condC;zoneC=zoneC; glb=glb; zoneDD=zoneDD; ctrl=ctrl;
  # varyrs=6;lastsigR=0.1; multTAc=multTAC
  if (ctrl$randseedP > 0) set.seed(ctrl$randseedP)
  projyrs <- projC$projyrs
  projC <- modprojC(zoneC,glb,projC) # include selectivity into projC
  zoneCR <- modzoneCSel(zoneC,projC$Sel,projC$SelWt,glb,projyrs)
  zoneDR <- makezoneDP(projyrs,ctrl$reps,glb) #,zoneDep) # zoneDReplicates
  endyr <- nrow(zoneDD$matureB)
  endcatch <- tapply(zoneDD$catch[endyr,],glb$sauindex,sum,na.rm=TRUE)
  zoneDR$TAC[1,] <- sum(endcatch)
  arv <- addrecvar(zoneDD=zoneDD,zoneDP=zoneDR,zoneC=zoneC,glob=glb,
                   condC=condC,ctrl=ctrl,varyrs=varyrs,multTAC=multTAC,
                   calcpopC,lastsigR=lastsigR)
  return(list(zoneDP=arv$zoneDP,projC=projC,zoneCP=zoneCR,zoneDDR=arv$zoneDDR))
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

