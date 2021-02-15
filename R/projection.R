
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
#' @param aspcatch the aspirational catches for the first year of projection
#' @param nyrs defaults to 25, which appears to be enough to produce the
#'     required variation
#' @param reset should the first ten rows of the each output be set to the last
#'     10 rows of the initial variation run. default = TRUE. If FALSE the
#'     reps trajectories when exposed to recruitment variation are passed out
#'
#' @return a revised zoneDR object with the first row filled with the last row
#'     from addrecvar
#' @export
#'
#' @examples
#' print("wait on suitable data")
addrecvar <- function(zoneC,zoneDR,inHt,glb,ctrl,aspcatch,nyrs=25,reset=TRUE) {
  # zoneC=zoneC;zoneDR=zoneDR;inHt=zoneDR$harvestR;glb=glb;ctrl=ctrl;aspcatch=cmcda$acatch
  # nyrs=25; reset=TRUE
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
    varzoneDR$matureB[1,,] <- zoneDR$matureB[nyrs,,]
    varzoneDR$exploitB[1,,] <- zoneDR$exploitB[nyrs,,]
    varzoneDR$catch[1,,] <- zoneDR$catch[nyrs,,]
    varzoneDR$acatch[1,,] <- aspcatch
    varzoneDR$harvestR[1,,] <- zoneDR$harvestR[nyrs,,]
    varzoneDR$cpue[1,,] <- zoneDR$cpue[nyrs,,]
    varzoneDR$recruit[1,,] <- zoneDR$recruit[nyrs,,]
    varzoneDR$deplsB[1,,] <- zoneDR$deplsB[nyrs,,]
    varzoneDR$depleB[1,,] <- zoneDR$depleB[nyrs,,]
    varzoneDR$catchN[,1,,] <- zoneDR$catchN[,nyrs,,]
    varzoneDR$Nt[,1,,] <- zoneDR$Nt[,nyrs,,]
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

#' @title doprojection applies the Tasmanian MCDA to a simulated zone
#'
#' @description doprojection applies the Tasmanian MCDA to a simulated zone. It
#'     is still in need of optimization as 1000 iterations would currently take
#'     about 20 minutes.
#'
#' @param zoneCP the modified zoneC object ready for projections
#' @param zoneDP the modified zoneD object ready for projections
#' @param glob the global variables object
#' @param ctrl the ctrl object
#' @param projyrs the number of years of projection
#' @param applyHS the function describing the harvest strategy to be used
#' @param hsargs a list of any arguments used by applyHS
#' @param projpms the matrices of performance measures from calibrateHS
#' @param inityrs the number of years kept from the period of conditioning,
#' @param ... the ellipsis should contain any additional arguments
#'
#' @return the filled in zoneDP list
#' @export
#'
#' @examples
#' print("wait on data files")
#' # zoneCP=zoneCP;zoneDP=zoneDR;glob=glb;ctrl=ctrl;projyrs=projC$projyrs;inityrs=projC$inityrs
#' # wid=4;targqnt=0.55;pmwts=c(0.65, 0.25,0.1);hcr = c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2)
doprojection <- function(zoneCP,zoneDP,glob,ctrl,projyrs,applyHS,hsargs,
                         projpms,inityrs,...) {
  # get  important constants
  sigmaR <- ctrl$withsigR # needed to add recruitment variation
  npop <- glob$numpop
  nsau <- glob$nSAU
  Ncl <- glob$Nclass
  nyrs <- projyrs
  movem <- glob$move
  reps <- ctrl$reps
  sauindex <- glob$sauindex
  matb <- numeric(npop)                   # use the same initial TAC for all reps
  origTAC <- mean(colSums(zoneDP$catch[1,,])) # mean sum of catches in last year

  saucatch <- array(0,dim=c(nyrs,nsau,reps))
  saucpue <- saucatch

  # for (iter in 1:reps) { # generate SAU total catches and catch-weighted cpue iter=1
  #   for (yr in 1:inityrs) { # iter=1; yr=4
  #     saucatch[yr,,iter] <- tapply(zoneDP$catch[yr,,iter],sauindex,sum,na.rm=TRUE)
  #     wts <- zoneDP$catch[yr,,iter]/(saucatch[yr,sauindex,iter])
  #     saucpue[yr,,iter] <- tapply((zoneDP$cpue[yr,,iter] * wts),sauindex,sum,na.rm=TRUE)
  #     if (yr >= hsargs$wid) {
  #       grad1[yr,,iter] <- apply(saucpue[1:yr,,iter],2,getgradone,yr=yr)
  #       grad4[yr,,iter] <- apply(saucpue[1:yr,,iter],2,getgradwid,yr=yr,wid=hsargs$wid)
  #     }
  #   }
  #   for (sau in 1:nsau)
  #     targetce[sau] <- targscore(saucpue[1:inityrs,sau,iter],qnt=hsargs$targqnt)
  #   targsc[1:inityrs,,iter] <- targetce
  # }
  # now do replicates, updating saucatch and saucpue each year
  if(ctrl$randseedP > 0) set.seed(ctrl$randseedP) # set random seed if desired
  for (iter in 1:reps) {
    TAC <- origTAC  # use the same original TAC for each replicate
    for (year in (inityrs+1):nyrs) { # iter=1; year=11
      #  catpop <- colSums(zoneDP$catch[1:(year - 1),,iter])
      inexpB <- zoneDP$exploitB[(year - 1),,iter]
      sauexpB <- tapply(inexpB,sauindex,sum,na.rm=TRUE)
      catbysau <- TAC * sauexpB/sum(sauexpB)  # no error initially
      multh <- apply(saucpue[1:(year-1),,1],2,applyHS,yr=(year-1)) # apply mcdahcr
      TAC <- sum(catbysau * multh)
      divererr <- sauexpB * exp(rnorm(nsau,mean=0,sd=ctrl$withsigB))
      catbysau <- TAC * (divererr/sum(divererr)) # currently no error on TAC
      catbypop <- catbysau[sauindex] * (inexpB/sauexpB[sauindex]) # no error on pops
      for (popn in 1:npop) { # year=11; iter=1; pop=1
        out <- oneyearcat(inpopC=zoneCP[[popn]],inNt=zoneDP$Nt[,year-1,popn,iter],
                          Nclass=Ncl,incat=catbypop[popn],yr=year)
        zoneDP$exploitB[year,popn,iter] <- out$ExploitB
        zoneDP$matureB[year,popn,iter] <- out$MatureB
        zoneDP$catch[year,popn,iter] <- out$Catch
        zoneDP$harvestR[year,popn,iter] <- out$Harvest
        zoneDP$cpue[year,popn,iter] <- out$ce
        zoneDP$Nt[,year,popn,iter] <- out$Nt
        zoneDP$catchN[,year,popn,iter] <- out$CatchN
        matb[popn] <- out$MatureB
      } # pop
      steep <- getvect(zoneCP,"steeph")
      r0 <- sapply(zoneCP,"[[","R0")
      b0 <- sapply(zoneCP,"[[","B0")
      recs <- oneyearrec(steep,r0,b0,matb,sigR=ctrl$withsigR)
      newrecs <- movem %*% recs
      zoneDP$recruit[year,,iter] <- newrecs
      zoneDP$Nt[1,year,,iter] <- newrecs
      zoneDP$deplsB[year,,iter] <- zoneDP$matureB[year,,iter]/b0
      zoneDP$depleB[year,,iter] <- zoneDP$exploitB[year,,iter]/sapply(zoneCP,"[[","ExB0")
      saucatch[year,,iter] <- tapply(zoneDP$catch[year,,iter],sauindex,sum,na.rm=TRUE)
      wts <- zoneDP$catch[year,,iter]/(saucatch[year,sauindex,iter])
      saucpue[year,,iter] <- tapply((zoneDP$cpue[year,,iter] * wts),sauindex,sum,na.rm=TRUE)
    }   # year loop        zoneDR$matureB[,,1]
  }     # rep loop
  zoneDP$cesau <- saucpue
  return(zoneDP=zoneDP)
} # end of doprojection


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
  N <- glb$Nclass
  SAU <- glb$SAUpop
  namedims <- list(1:projyr,1:numpop,1:iter)
  namendims <- list(glb$midpts,1:projyr,1:numpop,1:iter)
  MatB <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  ExplB <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  catch <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  acatch <- array(0,dim=c(projyr,nSAU,iter),
                  dimnames=list(1:projyr,1:nSAU,1:iter)) #aspirational catches
  harvest <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  cpue <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  cesau <- array(0,dim=c(projyr,nSAU,iter),
                 dimnames=list(1:projyr,1:nSAU,1:iter))
  catsau <- array(0,dim=c(projyr,nSAU,iter),
                  dimnames=list(1:projyr,1:nSAU,1:iter))
  recruit <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  deplSpB <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  deplExB <- array(0,dim=c(projyr,numpop,iter),dimnames=namedims)
  catchN <- array(data=0,dim=c(N,projyr,numpop,iter),dimnames=namendims)
  Nt <- array(data=0,dim=c(N,projyr,numpop,iter),dimnames=namendims)
  zoneDP <- list(SAU=SAU,matureB=MatB,exploitB=ExplB,catch=catch,acatch=acatch,
                 harvestR=harvest,cpue=cpue,cesau=cesau,catsau=catsau,
                 recruit=recruit,deplsB=deplSpB,depleB=deplExB,catchN=catchN,
                 Nt=Nt)
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
#'     projC
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
  popdefs <- getlistvar(zoneC,"popdef")
  sau <- getvar(zoneC,"SAU")
  midpts <- glob$midpts
  projSel <- matrix(0,nrow=glob$Nclass,ncol=numpop,dimnames=list(midpts,sau))
  projSelWt <- matrix(0,nrow=glob$Nclass,ncol=numpop,dimnames=list(midpts,sau))
  pLML <- projC$projLML[1] # needs development to allow variation in projLML
  for (pop in 1:numpop) {
    selL50 <- popdefs["SelP1",pop]
    selL95 <- popdefs["SelP2",pop]
    projSel[,pop] <- logistic((pLML + selL50),selL95,midpts)
    projSelWt[,pop] <- projSel[,pop] * zoneC[[pop]]$WtL
  }
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

#' @title poptosau combines population cpue into sau as catch weighted sums
#'
#' @description poptosau combines the cpue from separate populations into their
#'     respective sau using a catch-weighted strategy. The sauindex is used to
#'     identify which populations to apply the sau total catches to.
#'
#' @param catvect the vector of catches x population for a given year
#' @param cpuevect the vector of cpue x population for a given year
#' @param sauindex the sau indices of each population
#'
#' @return a list of saucpue and saucatch
#' @export
#'
#' @examples
#' print("wait on appropriate built-in data files")
poptosau <- function(catvect,cpuevect,sauindex) {
  outvect <- tapply(catvect,sauindex,sum,na.rm=TRUE)
  wts <- catvect/outvect[sauindex]
  saucpue <- tapply((cpuevect * wts),sauindex,sum,na.rm=TRUE)
  return(list(saucpue=saucpue,saucatch=outvect))
} # end of poptosau


#' @title prepareprojection high level function that sets up a projection
#'
#' @description prepareprojection is a high level function that restructures
#'     zoneC by including the projected selectivity (in case the LML changes),
#'     and also selectivity x weight-at-size (for computational speed). It
#'     then generates a new zoneD with room for all replicates as well as the
#'     aspirational catches from each HS. Finally, it adds recruitment variation
#'     to each of the replicates and keeps the last year of each iteration of
#'     the addrecvar function as the start of each replicate projection, with
#'     zeros in the catch, cpue, and cesau arrays
#'
#' @param zone1 the original zone1 object from readctrlfile
#' @param zoneC the constant part of the zone
#' @param glb the global variables
#' @param zoneDep the zone after initial depletion
#' @param aspcatch the aspirational catch for the first year of every projection
#'
#' @return a list of the dynamic zone object as a list of arrays of projyrs x
#'     populations x replicates, plus the revised projC and revised zoneC
#' @export
#'
#' @examples
#' print("wait on data files")
prepareprojection <- function(zone1,zoneC,glb,zoneDep,aspcatch) {
  # zone1=zone$zone1;zoneC=zone$zoneC; glb=glb; zoneDep=zoneDD; aspcatch=cmcda$acatch
  ctrl <- zone1$ctrl
  projC <- modprojC(zoneC,glb,zone1$projC) # include selectivity into projC
  projyrs <- projC$projyrs
  zoneC <- modzoneCSel(zoneC,projC$Sel,projC$SelWt,glb,projyrs)
  # need to add aspirational catches to zoneDRp
  zoneDR <- makezoneDP(projyrs,ctrl$reps,glb) #,zoneDep) # zoneDReplicates
 # zoneDRp <- addrecvar(zoneC,zoneDR,zoneDR$harvestR,glb,ctrl,aspcatch)
  return(list(zoneDP=zoneDR,projC=projC,zoneCP=zoneC))
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

