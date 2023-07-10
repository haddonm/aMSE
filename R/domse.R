# Tas context -----------------
# postfixdir <- "BCconstref"
# rundir <- rundir
# controlfile=controlfile
# hsargs=hsargs
# hcrfun=constantrefhcr  #mcdahcr  # consthcr   #constantrefhcr
# sampleCE=tasCPUE
# sampleFIS=tasFIS
# sampleNaS=tasNaS
# getdata=tasdata
# calcpopC=calcexpectpopC
# varyrs=7
# startyr=38
# verbose=TRUE
# ndiagprojs=4
# makehcrout=makeoutconst   #makeouthcr
# HSPMs=getcpueHS
# fleetdyn=NULL
# cutcatchN=56
# matureL = c(70,200)
# wtatL = c(80,200)
# mincount=100
# includeNAS=FALSE
# depensate=0
# pmwtSwitch = 0
# stablewts = c(0.4, 0.5, 0.1)
# hcrname="constantrefhcr"   #"consthcr"
# hcrscoreoutputs = extractandplotscores
# kobeRP = c(0.4,0.2,0.15)
# interimout=""

# SA context -------------
# rundir=rundir
# controlfile=controlfile
# hsargs=hsargs
# hcrfun=SA_HS
# sampleCE=saCPUE
# sampleFIS=saFIS
# sampleNaS=saNaS
# getdata=sadata
# calcpopC=sacalcexpectpopC
# makeouthcr=samakeouthcr
# makehcrout=makeouthcr
# hcrscoreoutputs = saextractandplotscores
# HSPMs=sagetHSscores
# fleetdyn = safleetdyn
# varyrs=7
# startyr=38
# verbose=TRUE
# ndiagprojs=3
# cutcatchN=56
# matureL=c(40,170)
# wtatL=c(50,200)
# mincount=120
# includeNAS = TRUE
# depensate=0
# kobeRP=c(0.4,0.2,0.15)
# hcrfun=SA_HS
# sampleCE=saCPUE
# sampleFIS=saFIS
# sampleNaS=saNaS
# getdata=sadata
#


#' @title do_MSE an encapsulating function to hold the details of a single run
#'
#' @description do_MSE is a form of meta-function that contains the code and
#'     function calls that constitute a single scenario run for aMSE. This is
#'     the reason the argument list is as long as it is. The code has been
#'     encapsulated like this to simplify development, maintenance, and
#'     distribution. Already with the number of variant conditioning examples
#'     produced as tests of the code-base, errors have arisen through a failure
#'     to propagate all changes to the code now in this function to all the
#'     separate MSE run files used to manage each scenario. A potential flaw
#'     lies with the need to apply the jurisdictionHS functions to any
#'     historical fishery data used to condition the operating model, ready for
#'     the first year of the projections. Currently, the code used is unique to
#'     the Tasmanian case and this obviously still requires generalization.
#'     This may require a new function to be produced in each jurisdictionHS.R
#'     file, but that is still under consideration. Part of the development of
#'     this meta-function will be to further articulate the generation of
#'     results, and the inclusion of error traps on the argument list entries.
#'     To aid in the conditioning the option of not calculating the productivity
#'     has been added. By default it will occur and when projecting it will also
#'     always occur. HSstats is a list that is always saved and contains
#'     performance statistics for the HS used (currently contains sum5 and sum10
#'     the cumulative catches to 5 and 10 years into the projections).
#'
#' @param rundir the full path to the directory in which all files relating to a
#'     particular run are to be held.
#' @param controlfile the filename of the control file present in rundir
#'     containing information regarding the particular run.
#' @param hsargs the constants used to define the workings of the hcr
#' @param hcrfun the name of the harvest control rule that is used to
#'     calculate the aspirational catches by SAU or TAC by zone for the
#'     following year
#' @param sampleCE a function from jurisdictionHS.R that generates the CPUE
#'     statistics from the projected data
#' @param sampleFIS a function from jurisdictionHS.R that generates the FIS
#'     statistics
#' @param sampleNaS a function from jurisdictionHS.R that generates the
#'     Numbers-at-size samples
#' @param getdata a function that gathers all the data required by the hcrfun
#'     and combines it into an hcrdata object ready for the hcrfun
#' @param calcpopC a function that takes the output from hcrfun (either the
#'     aspirational catch x SAU or TAC x zone) and generates the actual catch
#'     per population in each SAU expected in the current year.
#' @param makeouthcr Another JurisdictionHS.R function that generates (or not)
#'     an object that is continually updated by the hcrfun. To avoid very
#'     inefficient code this object (at least for Tasmania) is cycled through
#'     the iterations.
#' @param hcrscoreoutputs A function defined in an extra aMSE package or
#'     source file, that defines how the hcr scores should be dealt with. In
#'     Tasmania they recreated from the projections, plotted and their values
#'     snet out into the scores object.
#' @param HSPMs Another function defined in the HS package that reconstructs
#'     the cpue HSPM scores, and other measures ready for plotting
#' @param fleetdyn a function that calculates the distribution of catch across
#'     the sau and populations. Currently not needed by Tas but needed by SA
#' @param interimout should results be saved after projections have finished?
#'     the default="", which means nothing is saved if no risk of crashing in
#'     subsequent analyses. But if trying new stuff which may waste lots of time
#'     if the projections need to be repeated, the set interimout = c(outdir,
#'     postfixdir), which should be defined in the run_aMSE you are using,
#'     and that will save the main objects needed for later analysis in the
#'     directory used to save final results. The ordering of
#' @param varyrs how many years at the end of the conditioning on the fishery,
#'     data into zoneDD, to which to add recruitment variation, default = 7,
#'     which suits the Tasmanian west coast. Used in prepareprojection
#' @param startyr the index for the first year of the conditioning dynamics to
#'     include in the plotted cpue trajectories, to give startyr:indexoflastyear.
#'     used in sauplots.
#' @param verbose Should progress comments be printed to console, default=FALSE
#' @param ndiagprojs the number of replicate trajectories to plot in the
#'     diagnostics tab to illustrate ndiagprojs trajectories to ensure that
#'     such projections appear realistic; default=3
#' @param cutcatchN to reduce the size of the final array of numbers-at-size
#'     in the catch one can remove all the empty cells below a given size
#'     class. In the default there are 105 2mm size classes and setting
#'     cutcatchN = 56 removes size classes 1:55 2 - 110 mm; only from catchN
#' @param matureL is a vector of 2, default = c(70,200), that define the x-axis
#'     bounds on the biology maturity-at-Length plots
#' @param wtatL is a vector of 2, default = c(80,200), that define the x-axis
#'     bounds on the biology weight-at-Length plots
#' @param mincount given size-composition data what minimum sample size will
#'     be deemed acceptable for inclusion in the plots and conditioning.
#'     default=100.
#' @param includeNAS should the NAS projections be included in the output. If
#'     they are then the final output size is greatly increased. default=FALSE
#' @param depensate should depensation occur. If > 0, eg depensate=0.2, then
#'     depensation will occur when depletion is < whatever value one puts in for
#'     depensate. Implemented inside the function 'oneyearrec'
#' @param kobeRP reference points for the eHS pseudo kobe plot. A vector of the
#'     targdepl, default=0.4, limdepl, default=0.2, and the limH, default=0.15
#'
#' @seealso{
#'  \link{makeequilzone}, \link{dohistoricC}, \link{prepareprojection},
#'  \link{doprojections}, \link{getprojyrC}, \link{oneyearrec}
#' }
#'
#' @return a large list containing tottime, projtime, starttime, glb, ctrl,
#'     zoneDD, zoneDP, projC, condC, sauout, and outzone
#' @export
#'
#' @examples
#' print("wait on suitable data sets in data")
do_MSE <- function(rundir,controlfile,hsargs,hcrfun,sampleCE,sampleFIS,
                   sampleNaS,getdata,calcpopC,makeouthcr,hcrscoreoutputs,
                   HSPMs,fleetdyn=NULL,interimout="",
                   varyrs=7,startyr=42,verbose=FALSE,ndiagprojs=3,
                   cutcatchN=56,matureL=c(70,200),wtatL=c(80,200),mincount=100,
                   includeNAS=FALSE,depensate=0,kobeRP=c(0.4,0.2,0.15)) {
  # generate equilibrium zone -----------------------------------------------
  starttime <- (Sys.time())
  zone <- makeequilzone(rundir,controlfile,verbose=verbose)
  equiltime <- (Sys.time()); if (verbose) print(equiltime - starttime)
  # declare main objects ----------------------------------------------------
  glb <- zone$glb
  ctrl <- zone$ctrl
  zoneC <- zone$zoneC
  zoneD <- zone$zoneD
  zone1 <- zone$zone1
  saudat <- zone$saudat
  constants <- zone$constants
  production <- zone$product
  projC <- zone$zone1$projC
  condC <- zone$zone1$condC
  # save some equil results -------------------------------------------------
  biology_plots(rundir, glb, zoneC, matL=matureL,Lwt=wtatL)
  sauprod <- plotproductivity(rundir,production,glb)
  numbersatsize(rundir, glb, zoneD)
  #Condition on Fishery -----------------------------------------------------
  if (any(condC$initdepl < 1)) {
    initdepl <- condC$initdepl
    if (verbose) cat("Conducting initial depletions  ",initdepl,"\n")
    zoneD <- depleteSAU(zoneC,zoneD,glb,initdepl=initdepl,production)
  }
  if (verbose) cat("Conditioning on the Fishery data  \n")
  zoneDD <- dohistoricC(zoneD,zoneC,glob=glb,condC,calcpopC=calcpopC,
                        fleetdyn=fleetdyn,hsargs=hsargs,sigR=1e-08,sigB=1e-08)
  hyrs <- glb$hyrs
  propD <- as.data.frame(t(getzoneprops(zoneC,zoneDD,glb,year=hyrs)))
  columns <- c("B0","MSY","MSYDepl","bLML","propprot","SpBDepl","catch",
               "harvestR")
  plotpopprops(x=propD,rundir=rundir,glb=glb,varnames=columns,startyr=hyrs,
               bins=21,console=FALSE)
  propD[,"SAU"] <- c(glb$sauname[glb$sauindex],NA)
  addtable(propD,"propertyDD.csv",rundir,category="zoneDD",caption=
             "Properties of zoneD after conditioning on historical catches.")
  addtable(round(t(zoneDD$harvestR[(hyrs-9):hyrs,]),4),"final_harvestR.csv",
           rundir,category="zoneDD",
           caption="Last ten years of population vs harvest rate.")
  addtable(round(saudat,5),"saudat.csv",rundir,category="zoneDD",
           caption="SAU constant definitions")
  popdefs <- as.data.frame(t(getlistvar(zoneC,"popdef")))
  popdefs[,c(1:6,8:18)] <- round(popdefs[,c(1:6,8:18)],3) #so SAU can be txt
  popdefs[,"SAU"] <- glb$sauname[glb$sauindex]  #so SAU can be txt
  addtable(popdefs,"popdefs.csv",rundir,category="zoneDD",
           caption="Population vs Operating model parameter definitions")
  condout <- plotconditioning(zoneDD,glb,zoneC,condC$histCE,
                              histCatch=condC$histCatch,rundir,
                              condC$recdevs,console=FALSE)
  saurecdevs(condC$recdevs,glb,rundir,filen="saurecdevs.png")
  # plot predicted size-comp of catch vs observed size-comps
  catchN <- zoneDD$catchN
  sauCt <- popNAStosau(catchN,glb)
  compdat <- condC$compdat$lfs
  if (!is.null(compdat)) {
    for (plotsau in 1:glb$nSAU) {
      lfs <- preparesizecomp(compdat[,,plotsau],mincount=mincount)
      yrsize <- as.numeric(colnames(lfs))
      histyr <- condC$histyr
      pickyr <- match(yrsize,histyr[,"year"])
      LML <- histyr[pickyr,"histLML"]
      plotsizecomp(rundir=rundir,incomp=lfs,SAU=glb$saunames[plotsau],lml=LML,
                   catchN=sauCt[,,plotsau],start=NA,proportion=TRUE,
                   console=FALSE)
    }
  }
  # do projections ------------------------------------------------------------
  if (verbose) cat("Preparing for the projections \n")
  outpp <- prepareprojection(projC=projC,condC=condC,zoneC=zoneC,glb=glb,
                             calcpopC=calcpopC,zoneDD=zoneDD,
                             ctrl=ctrl,varyrs=varyrs,lastsigR = ctrl$withsigR)
  zoneDP <- outpp$zoneDP
  projC <- outpp$projC
  zoneCP <- outpp$zoneCP
 # Rprof()  #  ctrl$reps=25
  if (verbose) cat("Doing the projections \n")
  outproj <- doprojections(ctrl=ctrl,zoneDP=zoneDP,zoneCP=zoneCP,glb=glb,
                           hcrfun=hcrfun,hsargs=hsargs,sampleCE=sampleCE,
                           sampleFIS=sampleFIS,sampleNaS=sampleNaS,
                           getdata=getdata,calcpopC=calcpopC,
                           makehcrout=makeouthcr,fleetdyn=fleetdyn,verbose=TRUE)

#  Rprof(NULL)
  if (verbose) cat("Now generating final plots and tables \n")
  zoneDP=outproj$zoneDP
  hcrout <- outproj$hcrout; #str(hcrout)
  NAS <- list(Nt=zoneDP$Nt,catchN=zoneDP$catchN)
  zoneDP <- zoneDP[-c(17,16,15)]  # This removes the Nt etc from zoneDP
  if (length(interimout) == 2) {
    outfile <- filenametopath(interimout[1],paste0("temp",interimout[2],".Rdata"))
    if (includeNAS) {
      postprojout <- list(zoneCP=zoneCP,hcrout=hcrout,zoneDP=zoneDP,glb=glb,
                         ctrl=ctrl,condC=condC,constants=constants,NAS=NAS)
    } else {
      postprojout <- list(zoneCP=zoneCP,hcrout=hcrout,zoneDP=zoneDP,glb=glb,
                         ctrl=ctrl,condC=condC,constants=constants)
    }
    save(postprojout,file=outfile)
  }
 # histCE <- condC$histCE
  B0 <- getvar(zoneC,"B0")
  ExB0 <- getvar(zoneC,"ExB0")
  sauout <- sauplots(zoneDP,NAS,glb,rundir,B0,ExB0,
                     startyr=startyr,addCI=TRUE,histCE=condC$histCE)
  diagnosticsproj(sauout,glb,rundir,nrep=ndiagprojs)
  outzone <- poptozone(zoneDP,NAS,glb,
                       B0=sum(getvar(zoneC,"B0")),
                       ExB0=sum(getvar(zoneC,"ExB0")))
  zonesummary <- plotZone(outzone,rundir,glb,startyr=startyr,
                          CIprobs=c(0.05,0.5,0.95),addfile=TRUE)
  fishery_plots(rundir=rundir,glb=glb,select=zoneCP[[1]]$Select,
                histyr=condC$histyr,projLML=projC$projLML)
  historicalplots(rundir=rundir,condC=condC,glb=glb)
  NAS$catchN <- NAS$catchN[(cutcatchN:glb$Nclass),,,]
  projtime <- Sys.time()
  tottime <- round((projtime - starttime),3)
  # calculate HS performance statistics
  sum5 <- getprojyrC(catsau=zoneDP$catsau,glb=glb,period=5)
  sum10 <- getprojyrC(catsau=zoneDP$catsau,glb=glb)
  HSstats <- list(sum10=sum10,sum5=sum5)
  save(HSstats,file=paste0(rundir,"/HSstats.RData"))
  save(glb,file=paste0(rundir,"/glb.RData"))
  save_hsargs(rundir,hsargs)   # prints hsargs to HSPerfs tab
  if (verbose) cat("HSstats.RData, glb.RData, and hsargs.txt saved to rundir \n")
  plothsstats(rundir,HSstats,glb,average=FALSE)
  plothsstats(rundir,HSstats,glb,average=TRUE)
  addtable(hcrout$refpts,"hcrout_refpts.csv",rundir,category="HSperf",
           caption="HCR reference points")
  scores <- hcrscoreoutputs(rundir=rundir,HSPMs=HSPMs,hcrout=hcrout,
                            cpue=sauout$cpue,catches=sauout$catch,
                            glb=glb,yearCE=condC$yearCE,hsargs=hsargs)
  # generate sau phase plots
  nSAU <- glb$nSAU
  kobedata <- vector(mode="list",length=nSAU)
  names(kobedata) <- glb$saunames
  for (plotsau in 1:glb$nSAU) {
    kobedata[[plotsau]] <- HSphaseplot(dyn=sauout,glb=glb,sau=plotsau,
                                       rundir=rundir,startyr=condC$yearCE[1],
                                       console=FALSE,targdepl=kobeRP[1],
                                       limdepl=kobeRP[2],limH=kobeRP[3])
  }

  if (!includeNAS) NAS=NULL
  out <- list(tottime=tottime,projtime=projtime,starttime=starttime,glb=glb,
              ctrl=ctrl,zoneCP=zoneCP,zoneD=zoneD,zoneDD=zoneDD,zoneDP=zoneDP,
              NAS=NAS,projC=projC,condC=condC,sauout=sauout,outzone=outzone,
              hcrout=hcrout,production=production,condout=condout,
              HSstats=HSstats,saudat=saudat,constants=constants,hsargs=hsargs,
              sauprod=sauprod,scores=scores,zonesummary=zonesummary,
              kobedata=kobedata)
  return(out)
} # end of do_MSE








