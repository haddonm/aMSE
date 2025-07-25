
#' @title asSAU generates an output list of SAU values from the in put zoneDP
#'
#' @description asSAU takes the population based results in zoneDP and converts
#'     them into a list based around SAU. saucpue is already calculated during
#'     the application of the harvest strategy. Appears redundant
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
#'     appears to be redundant
#'
#' @param sauzone the output from the apply harvest strategy function for each
#'     jurisdiction, a large list of arrays or results
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
#'     population's weight-at-size relationship will be slightly different. Now
#'     appears to be redundant.
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

# removed just after calcpopC


#' @title doprojections conducts the replicate model runs for Tasmania
#'
#' @description doprojections conducts the replicate model runs for
#'     Tasmania using the mcdahcr and the input hsargs. Note the use of a list
#'     that takes up both the population and catch numbers-at-size and puts that
#'     into getdata. By referencing the iteration in the arrays being passed to
#'     getdata the object sizes are being reduced and hence, hopefully, the
#'     runtime will be reduced accordingly.
#'
#' @param ctrl the ctrl object from readctrlfile
#' @param zoneDP the object used to contain the dynamics from the replicate
#'     model runs
#' @param zoneCP the object used to contain the constants for each population
#'     used in model dynamics
#' @param glb the object containing the global constants for the given run
#' @param hcrfun the name of the harvest control rule that is used to
#'     calculate the multiplier for the previous aspirational catches (possibly
#'     for each SAU but possibly the TAC for the whole zone) so as to
#'     estimate the aspirational catches or TAC or the following year
#' @param hsargs the constants used to define the workings of the hcr
#' @param sampleCE a function that generates the CPUE statistics
#' @param sampleFIS a function that generates the FIS statistics
#' @param sampleNaS a function that generates the Numbers-at-size samples
#' @param getdata a function that gathers all the data required by the hcrfun
#'     and combines it into an hcrdata object ready for the hcrfun. It is
#'     expected to call sampleCE, sampleFIS, and sampleNAS, even if they only
#'     return NULL.
#' @param calcpopC a function that takes the output from hcrfun and generates
#'     the actual catch per population expected in the current year.
#' @param makehcrout is a function from HS.R that produces an object that is
#'     updated in each iteration by the hcrfun. If no such object is required
#'     then have a function that returns NULL.
#' @param fleetdyn a function that calculates the distribution of catch across
#'     the sau and populations. Currently not needed by Tas but needed by SA
#' @param verbose should the iterations be counted on the console?
#' @param fissettings an object containing settings used when calculating
#'     indices for the FIS within oneyearcat inside oneyearsauC
#' @param yearFIS a vector holding the years of the FIS survey
#' @param fisindexdata the fisindex data by SAU
#' @param ... the ellipsis used in case any of the functions hcrfun, sampleCE,
#'     sampleFIS, sampleNas, and getdata require extra arguments not included
#'     in the default named collection
#'
#' @seealso{
#'  \link{oneyearsauC}, \link[makehtml]{make_html}, \link{do_MSE}
#' }
#'
#' @return a list containing the full dynamics across all years zoneDP, and the
#'     final output from the HCR as hcrout and outhcr
#' @export
#'
#' @examples
#' print("wait on suitable internal data sets")
doprojections <- function(ctrl,zoneDP,zoneCP,glb,hcrfun,hsargs,
                          sampleCE,sampleFIS,sampleNaS,getdata,calcpopC,
                          makehcrout,fleetdyn,verbose=FALSE,fissettings=NULL,
                          yearFIS=NULL,fisindexdata=NULL,...) {
  # ctrl=ctrl; zoneDP=zoneDP; zoneCP=zoneCP; glb=glb; hcrfun=hcrfun; hsargs=hsargs
  # sampleCE=tasCPUE; sampleFIS=tasFIS; sampleNaS=tasNaS;  getdata=constdata#tasdata
  # calcpopC=calcexpectpopC; verbose=TRUE; fleetdyn=NULL;makehcrout=makeoutconst#makeouthcr
  # fissettings=NULL;fisindexdata=NULL;yearFIS=NULL
  reps <- ctrl$reps
  sigmar <- ctrl$withsigR
  sigmab <- ctrl$withsigB
  sigce <- ctrl$withsigCE
  hyrs <- glb$hyrs
  startyr <- hyrs + 1
  pyrs <- glb$pyrs
  endyr <- hyrs + pyrs
  sauindex <- glb$sauindex
  saunames <- glb$saunames
  nsau <- glb$nSAU
  yrnames <- c(glb$hyrnames,glb$pyrnames)
  Nclass <- glb$Nclass
  movem <- glb$move
  r0 <- getvar(zoneCP,"R0") #R0 by population
  b0 <- getvar(zoneCP,"B0") #sapply(zoneC,"[[","B0")
  exb0 <- getvar(zoneCP,"ExB0")
  envyr <- NULL
  survNt <- NULL
  proprec <- NULL
  if (!is.null(glb$envimpact)) {
    envimpact <- glb$envimpact
    envyr <- envimpact[["eyr"]]
    neyr <- length(envyr)
    npop <- glb$numpop
    survNt <- matrix(0,nrow=neyr,ncol=npop,dimnames=list(envyr,1:npop))
    proprec <- survNt
    for (i in 1:neyr) {
      survNt[i,] <- envimpact[["propNt"]][i,sauindex]
      proprec[i,] <- envimpact[["proprec"]][i,sauindex]
    }
  }
  outhcr <- makehcrout(glb,hsargs)
  for (year in startyr:endyr) { # iter=1; year=59
    if (verbose) cat(yrnames[year]," ")
      for (iter in 1:reps) { # iter=1
      hcrdata <- getdata(sampleCE,sampleFIS,sampleNaS,
                         sauCPUE=as.matrix(zoneDP$cesau[,,iter]),
                         sauacatch=as.matrix(zoneDP$acatch[,,iter]),  # add year in here
                         saucatch=as.matrix(zoneDP$catsau[,,iter]),
                         sauNAS=list(Nt=zoneDP$Nt[,,,iter],
                         catchN=zoneDP$catchN[,,,iter],
                         NumNe=zoneDP$NumNe[,,,iter]),year=year,
                         startCE=hsargs$startCE,decrement=hsargs$decrement)
      hcrout <- hcrfun(hcrdata,hsargs,glb=glb,projyear=year,outhcr=outhcr,
                       iter=iter)
      outhcr <- hcrout$outhcr  # needed so it can be updated each iteration
      calcpopCout <- calcpopC(hcrout,exb=zoneDP$exploitB[year-1,,iter],
                       sauCPUE= zoneDP$cesau[year-1,,iter],
                       catsau = zoneDP$catsau[year-1,,iter],
                       fleetacatch=fleetdyn,hsargs=hsargs,glb=glb,
                       sigmab=sigmab,year=year)
      outy <- oneyearsauC(zoneCC=zoneCP,inN=as.matrix(zoneDP$Nt[,year-1,,iter]),
                          popC=calcpopCout$popC,year=year,Ncl=Nclass,
                          sauindex=sauindex,movem=movem,sigmar=sigmar,
                          sigce=sigce,r0=r0,b0=b0,exb0=exb0,envyr=envyr,
                          envsurv=survNt,envrec=proprec,deltarec=glb$deltarec,
                          fissettings=fissettings,fisindexdata=fisindexdata,
                          useF=glb$useF)
      dyn <- outy$dyn
      saudyn <- poptosauCE(dyn["catch",],dyn["cpue",],sauindex)
      zoneDP$exploitB[year,,iter] <- dyn["exploitb",]
      zoneDP$midyexpB[year,,iter] <- dyn["midyexpB",]
      zoneDP$matureB[year,,iter] <- dyn["matureb",]
      zoneDP$catch[year,,iter] <- dyn["catch",]
      zoneDP$acatch[year,,iter] <- calcpopCout$acatch
      zoneDP$catsau[year,,iter] <- saudyn$saucatch
      zoneDP$harvestR[year,,iter] <- dyn["catch",]/dyn["midyexpB",]
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
  cat("\n")
  return(list(zoneDP=zoneDP,outhcr=outhcr,hcrout=hcrout))
} # end of doprojections

#' @title makezoneDP generates the container for the projection dynamics
#'
#' @description makezoneDP generates an object designed to hold the outputs
#'     from each replicate within a set of projections. This is identical to
#'     zoneD except it contains a repeat for each iteration and now includes
#'     a cesau, which is the population-catch-weighted cpue for each SAU.
#'
#' @param allyr number of years of historical catches plus projection years
#' @param iter the number of replicates
#' @param glb the object containing the global variables
#'
#' @return a large list containing an object ready for the projection dynamics
#' @export
#'
#' @examples
#' print("Could add variation to the harvest rates so that when ")
#' print("prepareprojection was run the range of initial H values would increase ")
makezoneDP <- function(allyr,iter,glb) { #  projyr=nyrs;iter=reps;glb=glob
  numpop <- glb$numpop
  nSAU <- glb$nSAU
  saunames <- glb$saunames
  N <- glb$Nclass
  SAU <- glb$SAUnum
  yrnames <- c(glb$hyrnames,glb$pyrnames)
  namedims <- list(yrnames,1:numpop,1:iter)
  namendims <- list(glb$midpts,yrnames,1:numpop,1:iter)
  MatB <- array(0,dim=c(allyr,numpop,iter),dimnames=namedims)
  ExplB <- array(0,dim=c(allyr,numpop,iter),dimnames=namedims)
  midyexpB <- array(0,dim=c(allyr,numpop,iter),dimnames=namedims)
  catch <- array(0,dim=c(allyr,numpop,iter),dimnames=namedims)
  acatch <- array(0,dim=c(allyr,nSAU,iter),
                  dimnames=list(yrnames,saunames,1:iter)) #aspirational catches
  harvest <- array(0,dim=c(allyr,numpop,iter),dimnames=namedims)
  cpue <- array(0,dim=c(allyr,numpop,iter),dimnames=namedims)
  cesau <- array(0,dim=c(allyr,nSAU,iter),
                 dimnames=list(yrnames,saunames,1:iter))
  catsau <- array(0,dim=c(allyr,nSAU,iter),
                  dimnames=list(yrnames,saunames,1:iter))
  recruit <- array(0,dim=c(allyr,numpop,iter),dimnames=namedims)
  deplSpB <- array(0,dim=c(allyr,numpop,iter),dimnames=namedims)
  deplExB <- array(0,dim=c(allyr,numpop,iter),dimnames=namedims)
  catchN <- array(data=0,dim=c(N,allyr,numpop,iter),dimnames=namendims)
  Nt <- array(data=0,dim=c(N,allyr,numpop,iter),dimnames=namendims)
  NumNe <- array(data=0,dim=c(N,allyr,numpop,iter),dimnames=namendims)
  TAC <- array(0,dim=c(allyr,iter),dimnames=list(yrnames,1:iter))
  zoneDP <- list(SAU=SAU,matureB=MatB,exploitB=ExplB,midyexpB=midyexpB,catch=catch,
                 acatch=acatch,harvestR=harvest,cpue=cpue,cesau=cesau,
                 catsau=catsau,recruit=recruit,deplsB=deplSpB,depleB=deplExB,
                 TAC=TAC,Nt=Nt,NumNe=NumNe,catchN=catchN)
  return(zoneDP=zoneDP)
} # end of makezoneDP

#' @title modprojC produces three objects used to condition the zone
#'
#' @description modprojC produces an object used to condition the zone when
#'     projecting it following any conditioning. Once the initial conditions
#'     for the projection have been attained through an initial depletion of an
#'     unfished equilibrium, or more particularly by conditioning on historical
#'     data then the projC object will be used to condition the projections. In
#'     essence, this puts the selectivity and selectivity x weight-at-size into
#'     projC. It can handle changes in LML through the years of projection, and
#'     it can now also handle different LML in different sau as well as a slot
#'     fishery selectivity.
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
modprojC <- function(zoneC,glob,projC) { # zoneC=zoneC; glob=glb; projC=projC
  numpop <- glob$numpop
  midpts <- glob$midpts
  sauLML <- glob$sauLML
  projyrs <- projC$projyrs
  popdefs <- getlistvar(zoneC,"popdef")
  projSel <- array(0,dim=c(glob$Nclass,numpop,projyrs),
                   dimnames=list(midpts,glob$SAUnum,1:projyrs))
  projSelWt <- array(0,dim=c(glob$Nclass,numpop,projyrs),
                   dimnames=list(midpts,glob$SAUnum,1:projyrs))
  if (sauLML) {
    nsau <- glob$nSAU
    pLML <- projC$projLML[,(1:nsau)]
    mLML <- projC$projLML[,(nsau+1):(2*nsau)]
    nLML <- nrow(pLML)
  } else {
    pLML <- as.matrix(projC$projLML[,"LML"])
    mLML <- as.matrix(projC$projLML[,"Max"])
    nLML <- length(pLML)
  }
  if (nLML != projyrs) stop(cat("Only ",nLML,"projLML instead of ",projyrs,"\n"))
  for (pop in 1:numpop) { #  yr=1; pop=1
    selL50 <- popdefs["SelP1",pop]
    selL95 <- popdefs["SelP2",pop]
    if (sauLML) {
      pLMLcol <- grep(popdefs["SAU",pop],colnames(pLML)) # which column to use
    } else {
      pLMLcol <- 1
    }
    for (yr in 1:projyrs) {
      projSel[,pop,yr] <- logistic((pLML[yr,pLMLcol] + selL50),selL95,
                                     midpts,maxLML=mLML[yr,pLMLcol])
      projSelWt[,pop,yr] <- projSel[,pop,yr] * zoneC[[pop]]$WtL
    }
  }
  projC$Sel <- projSel
  projC$SelWt <- projSelWt
  return(projC=projC)
}

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

#' @title addrecvar adds recruitment variation to end of conditioning step
#'
#' @description addrecvar replicates the historical depletion step in zoneDD
#'     reps times and copies the first Nyrs - varyrs of values into each
#'     replicate. Then it steps through each replicate adding recruitment
#'     variation to the dynamics for the last varyrs years. This is so that
#'     when the projections begin each replicate will already have recruitment
#'     variability fully developed and there will be no time lag on the
#'     recruitment and subsequent dynamics. Take especial note that prior to the
#'     projections the proposed harvest strategy needs to be applied to the
#'     historical fishery data (in Tasmania this is currently catches and CPUE).
#'     This application occurs in lines 60 - 65 in the addrecvar function. This
#'     WILL NEED ATTENTION with other jurisdictions.
#'
#' @param zoneDD the dynamic zone after historical fishery data has been used
#'     to finalize the conditioning
#' #param zoneDDR the empty dynamic zone object ready for expansion and variation
#' @param zoneC the zone's constants object for each population
#' @param glob the object containing the global constants
#' @param condC the object containing the historical fishery data
#' @param ctrl the control object
#' @param varyrs the number of years at the end of the historical period to
#'     which recruitment variation is to be added
#' @param calcpopC a function that takes the output from hcrfun and generates
#'     the actual catch per population expected in the current year.
#' @param fleetdyn a function that calculates the distribution of catch across
#'     the sau and populations. Currently not needed by Tas but needed by SA
#' @param hsargs the constants used to define the workings of the hcr
#' @param sigR the initial recruitment variation default=1e-08
#' @param sigB the initial biomass cpue variation default = 1e-08
#' @param lastsigR the recruitment variation to be added to the final varyrs
#'
#' @return an initialized dynamics zone object for the projections with the
#'     first year populated
#' @export
#'
#' @examples
#' print("wait on suitable data-sets")
addrecvar <- function(zoneDD,zoneC,glob,condC,ctrl,varyrs,calcpopC,
                      fleetdyn,hsargs,sigR=1e-08,sigB=1e-08,lastsigR=0.3) {
  #  zoneDD=zoneDD;zoneC=zoneC;glob=glb; calcpopC=calcpopC
  #  condC=condC;ctrl=ctrl;varyrs=7;lastsigR=lastsigR;
  #  sigR=1e-08; sigB=1e-08; fleetdyn=NULL
  sauindex <- glob$sauindex
  histC <- condC$histCatch
  yrs <- condC$histyr[,"year"]
  hyrs <- glob$hyrs
  finalyr <- hyrs - varyrs
  pyrs <- glob$pyrs
  totyr <- hyrs + pyrs
  reps <- ctrl$reps
  zoneDDR <- makezoneDP(totyr,reps,glob)
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
  zoneDDR$deplsB[1:finalyr,,] <- zoneDD$deplsB[1:finalyr,]
  zoneDDR$depleB[1:finalyr,,] <- zoneDD$depleB[1:finalyr,]
  zoneDDR$catchN[,1:finalyr,,] <- zoneDD$catchN[,1:finalyr,]
  zoneDDR$Nt[,1:finalyr,,] <- zoneDD$Nt[,1:finalyr,]
  zoneDDR$NumNe[,1:finalyr,,] <- zoneDD$NumNe[,1:finalyr,]
  for (i in 1:finalyr) { # i = 1
    zoneDDR$acatch[i,,] <- tapply(zoneDD$catch[i,],sauindex,sum,na.rm=TRUE)
    saudyn <- poptosauCE(zoneDD$catch[i,],zoneDD$cpue[i,],sauindex)
    zoneDDR$cesau[i,,] <- saudyn$saucpue
    zoneDDR$catsau[i,,] <- saudyn$saucatch
    zoneDDR$TAC[i,] <- sum(saudyn$saucatch,na.rm=TRUE)
  }
  for (iter in 1:reps) { # iter = 1
    for (year in (finalyr+1):hyrs) {  # year = finalyr + 1
      catchsau <- histC[year,]        # Use actual catches
      hcrout <- list(acatch=catchsau) # attention needed in other jurisdictions
      histCE <- condC$histCE
      ceyrs <- condC$yearCE
      catchsau <- histC[year,]
      if (yrs[year] %in% ceyrs) {
        pickyr <- which(yrs[year] == ceyrs)
        cpuesau <- histCE[pickyr,]
      } else {
        cpuesau <- numeric(glob$nSAU)
      }
      calcpopCout <- calcpopC(hcrout,exb=zoneDDR$exploitB[year-1,,iter],
                       sauCPUE=cpuesau,catsau=catchsau,
                       fleetacatch=fleetdyn,hsargs=hsargs,
                       glb=glob,sigmab=sigB,year=year)
      inN <- as.matrix(zoneDDR$Nt[,year-1,,iter])
      out <- oneyearsauC(zoneCC=zoneC,inN=inN,popC=calcpopCout$popC,year=year,
                         Ncl=glob$Nclass,sauindex=sauindex,movem=glob$move,
                         sigmar=lastsigR,sigce=1e-08,r0=r0,b0=b0,exb0=exb0,
                         rdev=condC$recdevs,envyr=NULL,envsurv=NULL,
                         fissettings=NULL,fisindexdata=NULL,
                         useF=glob$useF)
      dyn <- out$dyn
      saudyn <- poptosauCE(dyn["catch",],dyn["cpue",],sauindex)
      zoneDDR$exploitB[year,,iter] <- dyn["exploitb",]
      zoneDDR$midyexpB[year,,iter] <- dyn["midyexpB",]
      zoneDDR$matureB[year,,iter] <- dyn["matureb",]
      zoneDDR$catch[year,,iter] <- dyn["catch",]
      zoneDDR$acatch[year,,iter] <- calcpopCout$acatch
      zoneDDR$harvestR[year,,iter] <- dyn["catch",]/out$dyn["exploitb",]
      zoneDDR$cpue[year,,iter] <- dyn["cpue",]
      zoneDDR$cesau[year,,iter] <- saudyn$saucpue
      zoneDDR$catsau[year,,iter] <- saudyn$saucatch
      zoneDDR$recruit[year,,iter] <- dyn["recruits",]
      zoneDDR$deplsB[year,,iter] <- dyn["deplsB",]
      zoneDDR$depleB[year,,iter] <- dyn["depleB",]
      zoneDDR$TAC[year,iter] <- sum(hcrout$acatch,na.rm=TRUE)
      zoneDDR$Nt[,year,,iter] <- out$NaL
      zoneDDR$catchN[,year,,iter] <- out$catchN
      zoneDDR$NumNe[,year,,iter] <- out$NumNe
    }
  }
  return(zoneDP=zoneDDR)
} # end of addrecvar


#' @title modzoneCSel changes the selectivity characteristics in zoneC
#'
#' @description modzoneCSel  changes the selectivity characteristics in zoneC
#'     which is necessary when making a projection under a different LML.
#'     This function assumes a constant LML during any projections.
#'
#' @param zoneC the constant zone object from setupzone
#' @param sel the new selectivity as a matrix of selectivity by population
#' @param selwt new selectivity x WtL as a matrix of SelWt x Population
#' @param glb the globals object
#'
#' @return the zoneC object modified ready for projection
#' @export
#'
#' @examples
#' print("wait on suitable data")
modzoneCSel <- function(zoneC,sel,selwt,glb) {
  numpop <- glb$numpop
  Nclass <- glb$Nclass
  for (pop in 1:numpop){ # pop = 1
    zoneC[[pop]]$Select <- cbind(zoneC[[pop]]$Select,sel[,pop,])
    zoneC[[pop]]$SelWt <- cbind(zoneC[[pop]]$SelWt,selwt[,pop,])
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
#'     replicate, thus starting the application of the HS. Finally, it adds
#'     recruitment variation to each of the replicates and keeps the last year
#'     of each iteration of the addrecvar function as the start of each
#'     replicate projection, with zeros in the catch, cpue, and cesau arrays.
#'     If the repeatable projections are wanted either set _randseedP_ = 0,
#'     whereupon the pseudorandom numbers will proceed from those started using
#'     the _randomseed_ set in the control file. Alternatively set a new seed
#'     into _randseedP_ (perhaps use _getseed()_), or set it to NA, which means
#'     it will merely run _set.seed()_, which will use the system time as a
#'     starting seed (hence each scenario run will differ).
#'
#' @param projC the projection object from readctrlfile
#' @param condC historical conditioning data
#' @param zoneC the constant part of the zone
#' @param glb the global variables
#' @param zoneDD the zone after initial depletion through conditioning on the
#'     fishery
#' @param ctrl the ctrl object for the scenario run
#' @param varyrs how many years at the end to add recruitment variation
#' @param calcpopC a function that takes the output from hcrfun and gernerates
#'     the actual catch per population expected in the current year.
#' @param lastsigR recruitment variation for when it is applied for varyrs
#'
#' @return a list of the dynamic zone object as a list of arrays of projyrs x
#'     populations x replicates, plus the revised projC and revised zoneC
#' @export
#'
#' @examples
#' print("wait on data files")
#' #  projC <- modprojC(zoneC,glb,projC) # modify selectivity and SelWt in projC
prepareprojection <- function(projC=projC,condC=condC,zoneC=zoneC,glb=glb,
                              zoneDD=zoneDD,ctrl=ctrl,varyrs=varyrs,
                              calcpopC=calcpopC,lastsigR = 0.3) {
  # projC=projC;condC=condC;zoneC=zoneC;glb=glb;zoneDD=zoneDD;ctrl=ctrl
  # calcpopC=calcpopC; varyrs=varyrs;lastsigR = ctrl$withsigR
  if (ctrl$randseedP > 0) set.seed(ctrl$randseedP)
  if (is.na(ctrl$randseedP)) set.seed()
#  projyrs <- glb$pyrs
  zoneCR <- zoneC #modzoneCSel(zoneC,projC$Sel,projC$SelWt,glb) #now zoneC
  zoneDP <- addrecvar(zoneDD=zoneDD,zoneC=zoneC,glob=glb,
                         condC=condC,ctrl=ctrl,varyrs=varyrs,
                         calcpopC,lastsigR=lastsigR)
  return(list(zoneDP=zoneDP,projC=projC,zoneCP=zoneCR))
} # end of prepareprojection


