#Tas context -----------------
# postfixdir <- "BCtestFIS" #
# rundir <- rundir
# controlfile=controlfile
# hsargs=hsargs
# hcrfun= constantrefhcr  #mcdahcr   #  # consthcr   #constantrefhcr
# sampleCE=tasCPUE
# sampleFIS=tasFIS
# sampleNaS=tasNaS
# getdata=tasdata    #  tasdata    constdata
# calcpopC=calcexpectpopC
# varyrs=7
# startyr=38
# verbose=TRUE
# ndiagprojs=4
# makeouthcr=makeouthcr    #makeouthcr   #makeoutconst - only used with consthcr#
# fleetdyn=NULL
# plotmultflags=plotmultandflags
# scoreplot = plotfinalscores
# cutcatchN=56
# matureL = c(70,200)
# wtatL = c(80,200)
# mincount=100
# includeNAS=FALSE
# depensate=0
# pmwtSwitch = 4
# stablewts = c(0.65, 0.25, 0.1)
# hcrname="constantrefhcr"  #"constantrefhcr"  #"mcdahcr"     #    #"consthcr"
# kobeRP = c(0.4,0.2,0.15)
# interimout=rundir
# nasInterval=5
# minsizecomp=c(100,135)
# uplimH=0.35
# incH=0.005
# deleteyrs=0
# calcpredfis=estpredfis
# # # # #

# SA context -------------
# rundir=rundir
# controlfile=controlfile
# hsargs=hsargs
# hcrfun=SA_HS # March 2024 file
# sampleCE=saCPUE
# sampleFIS=saFIS
# sampleNaS=saNaS
# getdata=sadata
# calcpopC=sacalcexpectpopC
# makeouthcr=samakeouthcr
# fleetdyn = fleetdynsa
# scoreplot=saplotfinalscores2
# plotmultflags=saplotmultandflags
# interimout= ""
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
# nasInterval=5
# minsizecomp=c(100,135)
# uplimH = 0.9
# incH = 0.01
# fisindex=NULL
# deleteyrs=0
#

# VIC HS hsargs
# rundir=rundir
# controlfile = controlfile
# hsargs=hsargs
# hcrfun=Vic_HS
# sampleCE=vicCPUE
# sampleFIS=vicFIS
# sampleNaS=vicNaS
# getdata=vicdata
# calcpopC=viccalcexpectpopC
# makeouthcr=vicmakeouthcr
# fleetdyn = NULL
# scoreplot=vicplotfinalscores2
# plotmultflags=vicplotmultandflags
# interimout= "" #outdir
# varyrs=7
# startyr=38
# verbose=TRUE
# ndiagprojs=3
# cutcatchN=56
# matureL=c(40,170)
# wtatL=c(50,200)
# mincount=120
# includeNAS = FALSE
# depensate=0
# kobeRP=c(0.4,0.2,0.15)
# nasInterval=5
# minsizecomp=c(100,135)
# uplimH=0.6
# incH=0.01




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
#'     an object that is continually updated by the hcrfun. The intent is that
#'     it will contain all the HS score values..
#' @param fleetdyn a function that calculates the distribution of catch across
#'     the sau and populations. Currently not needed by Tas but needed by SA
#' @param scoreplot the function that plots the hcr scores and final scores. In
#'     Tas = plotfinalscores, in SA = SAplotfinalscores. If more inputs than
#'     outhcr, zoneDP, sau, and filen=filename are needed then define your own
#'     function using the ellipsis at the end of the argument list.
#' @param plotmultflags the function that plots the TAC/acatch multipliers, and
#'     meta rule flags. In Tas = plotmultandflags. If more inputs than
#'     outhcr, zoneDP, sau, and filen=filename are needed then define your own
#'     function using the ellipsis at the end of the argument list.
#' @param interimout should results be saved after projections have finished?
#'     the default="", which means nothing is saved if no risk of crashing in
#'     subsequent analyses. But if trying new stuff which may waste lots of time
#'     if the projections need to be repeated, then set
#'     interimout = outdir, which should be a character string identifying the
#'     fullpath to a directory, which is usually defined in the run_aMSE you are
#'     using, and that will save the main objects needed for later analysis into
#'     the directory used to save final results. Importantly, if your simulation
#'     runs succeed there is no need to save the interimout results. Using
#'     interimout automatically saves the numbers-at-size so, if you want the
#'     numbers-at-size and your simulations work, then turn off interimout, ie
#'     set it = "", and set includeNAS=TRUE.
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
#' @param nasInterval default = 5, which selects the year just before
#'     projections and then the size-composition data for every such interval
#'     in years along the projection period. The years selected always include
#'     the last year of projection and the year just before projections started.
#'     Used by predsizecomp function inside do_MSE.
#' @param minsizecomp a vector of two values, the first being the minimum size-
#'     class to be used when plotting the predicted size-composition of the
#'     stock, and the second being the minimum size-class when plotting the
#'     numbers-at-size in the catch. The default=c(100,135)
#' @param uplimH defines the upper limit of harvest rate used when estimating
#'     the productivity (also important when initial depletion is not 1.0). The
#'     default = 0.4
#' @param incH defines the interval between H steps when estimating productivity
#'     default = 0.005
#' @param deleteyrs default = 0, meaning delete no years from the sizecomp data.
#'     if there are years that are to be removed then this should be a matrix
#'     of years to delete vs sau, ie years as rows and sau as columns. All
#'     sau need to be included. to complete the matrix 0 values can be used.
#' @param selectyr which year's LML should be used to estimate the LML. This
#'     should be a year index (eg in Tas 1-58). If set = 0, then the LML for
#'     the glb$hyrs + glb$pyrs, the last year of projections data will be used.
#'     default = 0
#' @param calcFIS default = NULL. jurisdictional function required to calculate
#'     the annual fisindex, only used when FIS data available
#' @param prepareforfis default = NULL. jurisdictional function required to
#'     prepare for using FIS data, only used when FIS data available
#' @param makefisproj default = NULL. jurisdictional function required to
#'     prepare the constants needed for using FIS data in theprojections, only
#'     used when FIS data available.
#'
#' @seealso{
#'  \link{makeequilzone}, \link{dohistoricC}, \link{prepareprojection},
#'  \link{doprojections}, \link{getprojyrC}, \link{oneyearrec}
#' }
#'
#' @return a vary large list containing tottime, projtime, starttime, glb, ctrl,
#'     zoneDD, zoneDP, projC, condC, sauout, and outzone
#' @export
#'
#' @examples
#' print("wait on suitable data sets in data")
do_MSE <- function(rundir,controlfile,hsargs,hcrfun,sampleCE,sampleFIS,
                   sampleNaS,getdata,calcpopC,makeouthcr,fleetdyn=NULL,
                   scoreplot,plotmultflags,interimout="",varyrs=7,startyr=42,
                   verbose=FALSE,ndiagprojs=3,cutcatchN=56,matureL=c(70,200),
                   wtatL=c(80,200),mincount=100,includeNAS=FALSE,
                   depensate=0,kobeRP=c(0.4,0.2,0.15),nasInterval=5,
                   minsizecomp=c(100,135),uplimH=0.4,incH=0.005,
                   deleteyrs=0,selectyr=0,calcFIS=NULL,prepareforfis=NULL,
                   makefisproj=NULL) {
# GENERATE EQUILIBRIUM ZONE -----------------------------------------------
  starttime <- (Sys.time())
  setuphtml(rundir)
  if (interimout  == rundir) {
    outfile <- pathtopath(rundir,paste0("zone.Rdata"))
    if (file.exists(outfile)) {
        load(file=outfile)
      } else {
        zone <- makeequilzone(rundir,controlfile,doproduct=TRUE,uplimH=uplimH,
                              incH=incH,verbose=verbose)
        save(zone,file=outfile)
        if (verbose) cat("Interim output file saved \n")
    }
  } else {
    zone <- makeequilzone(rundir,controlfile,doproduct=TRUE,uplimH=uplimH,
                          incH=incH,verbose=verbose)
  }
  if (verbose) {
    incrtime1 <- Sys.time(); timeinc <- incrtime1 - starttime
    cat("makeequilzone ",timeinc,attr(timeinc,"units"),"\n")
  }
  # declare main objects ----------------------------------------------------
  glb <- zone$glb
  # if TIMEVARY = 1 modify the content of deltarec
  if (!is.null(glb$deltarec)) glb <- updatedeltarec(glb,varyrs)
  ctrl <- zone$ctrl
  zoneC <- zone$zoneC
  zoneD <- zone$zoneD
  zone1 <- zone$zone1
  saudat <- zone$saudat
  constants <- zone$constants
  production <- zone$product
  projC <- zone$zone1$projC
  condC <- zone$zone1$condC
  hsargs$forfis <- NULL
  if (glb$useFIS) { # develop sau based FIS metrics
    hsargs$forfis <- prepareforfis(condC$fisindexdata,glb,hsargs$fispar)
  }
  # biology, recruits, and production tabs -------------------------------------------------
  outbiol <- biology_plots(rundir, glb, zoneC, matL=matureL,Lwt=wtatL)
  filen <- "population_defined_properties.csv"
  addtable(condC$poprec,filen=filen,rundir=rundir,category="popprops",
           caption=paste0("Specific population properties defined rather than ",
                          "randomly allocated away from a mean."))
  sauprod <- plotproductivity(rundir,production,glb,hsargs)
 # hsargs$saumsy <- sauprod[3,]
  # numbers-at-size tab------------------------------------------------------
  numbersatsize(rundir, glb, zoneD)
# CONDITION THE OM---------------------------------
  #zoneDD tab -----------------------------------------------------
  if (any(condC$initdepl < 1)) {
    initdepl <- condC$initdepl
    if (verbose) cat("Conducting initial depletions  ",initdepl,"\n")
    zoneD <- depleteSAU(zoneC,zoneD,glb,initdepl=initdepl,production)
  }
  if (verbose) cat("Conditioning on the Fishery data  \n")
  zoneDD <- dohistoricC(zoneD,zoneC,glob=glb,condC,calcpopC=calcpopC,
                        fleetdyn=fleetdyn,hsargs=hsargs,sigR=1e-08,sigB=1e-08)
# cbind(zoneDD$fisindex[,5],zoneDD$predfis[,5])
  if (verbose) {
    incrtime2 <- Sys.time(); timeinc <- incrtime2 - incrtime1
    cat("dohistoricC ",timeinc,attr(timeinc,"units"),"\n")
  }
    #pop properties
  popprops <- as.data.frame(sapply(zoneC,"[[","popdef"))
  columns <- NULL
  for (sau in 1:glb$nSAU)
    columns <- c(columns,paste0(glb$saunames[sau],"-",1:glb$SAUpop[sau]))
  colnames(popprops) <- columns
  msy <- sapply(zoneC,"[[","MSY")
  bLML <- sapply(zoneC,"[[","bLML")
  B0 <- sapply(zoneC,"[[","B0")
  pops <- rbind(1:sum(glb$SAUpop),popprops,msy,B0,bLML)
  rownames(pops) <- c("popnum",rownames(popprops),"msy","B0","bLML")
  pops <- t(pops)
  label <- paste0(ctrl$runlabel," population properties. This is a bigtable",
                  " so, if you scroll across or down then use the topleft",
                  "  arrow to return to the home page.")
  addtable(round(pops,4),filen="popprops.csv",rundir=rundir,category="poptable",
           caption=label,big=TRUE)
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
  # condition tab----------------------------------------
  condout <- plotconditioning(zoneDD,glb,zoneC,condC$histCE,
                              histCatch=condC$histCatch,rundir,
                              condC$recdevs,console=FALSE)
  if (verbose) {
    incrtime1 <- Sys.time(); timeinc <- incrtime1 - incrtime2
    cat("Conditioning plots completed ",timeinc,attr(timeinc,"units") ,"\n")
  }
  saurecdevs(condC$recdevs,glb,rundir,filen="saurecdevs.png")
  # FIS tab------------------------------
  if (glb$useFIS) {
    fisout <- plotindices(zoneDD=zoneDD,glb=glb,condC=condC,hsargs=hsargs,
                          rundir=rundir,console=FALSE)
    plotfisfit(condC=condC,predfis=zoneDD$predfis,forfis=hsargs$forfis,
               rundir=rundir,console=FALSE)
    if (verbose) cat("FIS plots completed \n")
  }
  # predictedcatchN tab-----------------------------------------
  catchN <- zoneDD$catchN
  sauCt <- popNAStosau(catchN,glb)
  compdat <- condC$compdat$lfs
  if (!is.null(compdat)) {
    for (plotsau in 1:glb$nSAU) { # plotsau=1
      if (length(ncol(deleteyrs)) == 0) delyrs=0 else delyrs=deleteyrs[,plotsau]
      lfs <- preparesizecomp(compdat[,,plotsau],mincount=mincount,
                             deleteyears=delyrs)
      yrsize <- as.numeric(colnames(lfs))
      histyr <- condC$histyr
      pickyr <- match(yrsize,histyr[,"year"])
      LML <- condC$histyr[pickyr,]
      plotsizecomp(rundir=rundir,incomp=lfs,SAU=glb$saunames[plotsau],lml=LML,
                   catchN=sauCt[,,plotsau],start=NA,proportion=TRUE,
                   console=FALSE)
    }
  }
  # OrigComp tab--------------------------------------------
  nsau <- glb$nSAU  #    sampsize <- round(colSums(compdata,na.rm=TRUE),1)
  saunames <- glb$saunames
  for (sau in 1:nsau) {  # sau=1
    labely <- paste0("Size-Composition of Catches for ",saunames[sau])
    ans <- plotcompdata(compdata=compdat[,,sau],analysis=saunames[sau],
                        ylabel=labely,console=FALSE,outdir=rundir,
                        barcol="red",bordercol=1,horizline=140)
      addplot(ans$filename,rundir=rundir,category="OrigComp",ans$caption)
   }
  # popgrowth tab ------------------------------implied size-at-age
  popgrowth(rundir=rundir,zoneC=zoneC,glb=glb,console=FALSE,maxage=30,
            startsize= 2.0)
  # DO-PROJECTIONS ------------------------------------------------------------
  if (verbose) cat("Preparing for the projections \n")
  outpp <- prepareprojection(projC=projC,condC=condC,zoneC=zoneC,glb=glb,
                             zoneDD=zoneDD,ctrl=ctrl,varyrs=varyrs,
                             calcpopC=calcpopC,hsargs=hsargs,
                             lastsigR = ctrl$withsigR)
  zoneDP <- outpp$zoneDP

  ## new from CM
  temp <- zoneDP$acatch
  temp[] <- 0
  zoneDP$flagstate <- zoneDP$closedyrs <- zoneDP$recovyrs <- temp
  rm(temp)

  projC <- outpp$projC
  zoneCP <- outpp$zoneCP
  if (verbose) {
    incrtime2 <- Sys.time(); timeinc <- incrtime2 - incrtime1
    cat("Projection preparation completed ",timeinc,attr(timeinc,"units"),"\n")
  }
 # Rprof()  #  ctrl$reps=25
  if (verbose) cat("Doing projections \n")
  outproj <- doprojections(ctrl=ctrl,zoneDP=zoneDP,zoneCP=zoneCP,glb=glb,
                           hcrfun=hcrfun,hsargs=hsargs,sampleCE=sampleCE,
                           sampleFIS=sampleFIS,sampleNaS=sampleNaS,
                           getdata=getdata,calcpopC=calcpopC,
                           makehcrout=makeouthcr,fleetdyn=fleetdyn,verbose=TRUE,
                           makefisproj=makefisproj)
  if (verbose) {
    incrtime1 <- Sys.time(); timeinc <- incrtime1 - incrtime2
    cat("All projections finished ",timeinc,attr(timeinc,"units") ,"\n")
    cat("Now generating final plots and tables \n")
  }
 # GENERATE SCENARIO OUTPUTS------------------------
  zoneDP <- outproj$zoneDP
  hcrout <- outproj$hcrout
  outhcr <- outproj$outhcr
  NAS <- list(Nt=zoneDP$Nt,catchN=zoneDP$catchN,sauNumNe=zoneDP$sauNumNe)
  pickNAS <-
    sort(which(names(zoneDP) %in% c("Nt","catchN","NumNe","sauNumNe")),decreasing=TRUE)
  zoneDP <- zoneDP[-pickNAS]  # This removes the Nt etc from zoneDP
  if ((nchar(interimout) > 0) & (interimout != rundir)) {
    confirmdir(interimout,ask=verbose)
    outfile <- pathtopath(interimout[1],paste0("temp_projections.Rdata"))
    if (includeNAS) {
      postprojout <- list(zoneCP=zoneCP,hcrout=hcrout,zoneDP=zoneDP,glb=glb,
                         ctrl=ctrl,condC=condC,constants=constants,NAS=NAS)
    } else {
      postprojout <- list(zoneCP=zoneCP,hcrout=hcrout,zoneDP=zoneDP,glb=glb,
                         ctrl=ctrl,condC=condC,constants=constants)
    }
    save(postprojout,file=outfile)
    if (verbose) cat("Interim output file saved \n")
  }
  B0 <- getvar(zoneC,"B0")
  ExB0 <- getvar(zoneC,"ExB0")
  # add to FIS tab?------------------------------------
  if (glb$useFIS) {
    plotpredfis(predfis=zoneDP$predfis,glb=glb,forfis=hsargs$forfis,
                rundir=rundir,console=FALSE)
  }
  # projSAU tab----------------------------------------------------------
  if (verbose) cat("Starting the sau related plots \n")
  sauout <- sauplots(zoneDP,NAS,glb,rundir,B0,ExB0,startyr=startyr,
                     addCI=TRUE,histCE=condC$histCE,hlines=sauprod$sauprod)
  sauNAS <- list(Nt=sauout$Nt,catchN=sauout$catchN)
  sauout <- sauout[-c(12,11)]  # This removes the Nt and catchN from sauout
  if (verbose) {
    incrtime2 <- Sys.time(); timeinc <- incrtime2 - incrtime1
    cat("Finished all sau plots ",timeinc,attr(timeinc,"units"),"\n")
  }
  # catchN tab------------------------------------------------------
  if (verbose) cat("Starting size-composition plots \n")
  nsau <- glb$nSAU  # first the stock size
  for (sau in 1:nsau) # sau=1
    predsizecomp(sau=sau, NSC=sauNAS$Nt, glb=glb, minSL=minsizecomp[1],
                 interval=nasInterval,prop=TRUE,console=FALSE,rundir=rundir)
  for (sau in 1:nsau)
    prob <- predsizecomp(sau=sau, NSC=sauNAS$catchN, glb=glb,
                         minSL=minsizecomp[2],interval=nasInterval,prop=TRUE,
                         console=FALSE,rundir=rundir)
  if ((verbose) & (prob != "OK")) cat(prob,"\n")
  saucomp <- getsaunas(NAS$catchN,glb) # mean numbers-at-sizex yr across reps
  saunames <- glb$saunames
  for (sau in 1:nsau) {
    labely <- paste0("Mean Shell Length mm across reps, CatchN ",saunames[sau])
    prob <- plotcompdata(saucomp[cutcatchN:glb$Nclass,,sau],
                         analysis=paste0("CatchN ",saunames[sau]),
                         ylabel=labely,console=FALSE,outdir=rundir)
    addplot(filen=prob$filename,rundir=rundir,category="CatchN",
            caption=prob$caption)
  }
  if (verbose) {
    incrtime1 <- Sys.time(); timeinc <- incrtime1 - incrtime2
    cat("Finished size-composition plots",timeinc,attr(timeinc,"units"),"\n")
  }
  # DiagProj tab------------------------------------
  diagnosticsproj(sauout,glb,rundir,nrep=ndiagprojs)
  # zonescale tab ----------------------------------
  outzone <- poptozone(zoneDP,NAS,glb,
                       B0=sum(getvar(zoneC,"B0")),
                       ExB0=sum(getvar(zoneC,"ExB0")))
  zonesummary <- plotZone(outzone,rundir,glb,startyr=startyr,
                          CIprobs=c(0.05,0.5,0.95),addfile=TRUE)

  # Fishery tab-----------------------------------------------
  if (verbose) cat("Plotting fishery information \n")
  if (glb$sauLML) {
    for (sau in 1:glb$nSAU) { # sau = 1
      pickS <- which(glb$sauindex == sau)
      outfish <- fishery_plots(rundir=rundir,glb=glb,
                               select=zoneCP[[pickS[1]]]$Select,
                               histyr=condC$histyr,projLML=projC$projLML,
                               SAU=sau)
    }
  } else {
    outfish <- fishery_plots(rundir=rundir,glb=glb,select=zoneCP[[1]]$Select,
                             histyr=condC$histyr,projLML=projC$projLML)
  }
  LMLs <- outfish$outyrs$LMLs
  tabcap <- paste0(" The LML and maximum LML in each year where a change occurs.",
                   " Once a change in LML occurs it continues until the next ",
                   "change, or until the final year.")
  addtable(LMLs,"LMLs_through_time.csv",rundir,category="Fishery",caption=tabcap)
  historicalplots(rundir=rundir,condC=condC,glb=glb)
  NAS$catchN <- NAS$catchN[(cutcatchN:glb$Nclass),,,]
  # HSperf tab--------------------------------
  projtime <- Sys.time()
  tottime <- round((projtime - starttime),3)
  sum5 <- getprojyrC(catsau=zoneDP$catsau,glb=glb,period=5)
  sum10 <- getprojyrC(catsau=zoneDP$catsau,glb=glb)
  HSstats <- list(sum10=sum10,sum5=sum5)
  #save(HSstats,file=paste0(rundir,"/HSstats.RData"))
  #save(glb,file=paste0(rundir,"/glb.RData"))
  save(hsargs,file=pathtopath(rundir,"hsargs.RData"))
  if (verbose) cat("hsargs saved as RData file to rundir \n")
  scoremed <- NULL
  nSAU <- glb$nSAU
  if (hsargs$hcrname != "consthcr") {
    if (verbose) cat("plotting HS performance statistics \n")
    plothsstats(rundir,HSstats,glb,average=FALSE)
    plothsstats(rundir,HSstats,glb,average=TRUE)
    addtable(hcrout$refpts,"hcrout_refpts.csv",rundir,category="HSperf",
             caption="HCR reference points")
  # scores tab----------------------------------------------------
    for (sau in 1:nSAU) { # sau=7
      filen <- paste0("HS_Score_plots_",glb$saunames[sau],".png")
      filename <- pathtopath(rundir,filen)
      caption <- paste0("Harvest strategy scores and outputs for ",
                        glb$saunames[sau],". Blue line is median target CE.")
      scoremed <- scoreplot(outhcr,zoneDP,sau,filen=filename)
      addplot(filen,rundir=rundir,category="scores",caption)
      filen <- paste0("HS_Mult_MetaFlag_plots_",glb$saunames[sau],".png")
      filename <- pathtopath(rundir,filen)
      caption <- paste0("catchmult and meta flags for ",glb$saunames[sau])
      plotmultflags(outhcr=outhcr,sauans=sauout,sau=sau,filen=filename)
      addplot(filen,rundir=rundir,category="scores",caption)
      #medscores <- cbind(scoremed)
    }
  }
  # poplevelplots tab ---------------------------------
  # generate population plots
  if (verbose) cat("plotting Population level dynamics \n")
  popmedcatch <- vector(mode="list",length=nSAU)
  names(popmedcatch) <- glb$saunames
  popmedcpue <- vector(mode="list",length=nSAU)
  names(popmedcpue) <- glb$saunames
  popmeddepleB <- vector(mode="list",length=nSAU)
  names(popmeddepleB) <- glb$saunames
  for (sau in 1:nSAU) { # sau = 2
    saumed <- poplevelplot(rundir=rundir,sau=sau,popvar=zoneDP$catch,glb=glb,
                 label="Catch",console=FALSE)
    popmedcatch[[sau]] <- saumed
    saumed <- poplevelplot(rundir=rundir,sau=sau,popvar=zoneDP$cpue,glb=glb,
                           label="CPUE",console=FALSE)
    popmedcpue[[sau]] <- saumed
    saumed <- poplevelplot(rundir=rundir,sau=sau,popvar=zoneDP$depleB,glb=glb,
                           label="Depletion_ExpB",console=FALSE)
    popmeddepleB[[sau]] <- saumed
  }
  # phaseplot tab--------------------------------------------------
  nsau <- glb$nSAU
  columns <- c("meddepl","medH","medTargCE","medGrad4")
  kobedata <- matrix(0,nrow=nsau,ncol=4,dimnames=list(glb$saunames,columns))
  for (sau in 1:glb$nSAU) {
    kobedata[sau,] <- twophaseplots(dyn=sauout,glb=glb,outhcr=outhcr,sau=sau,
                                    kobeRP=kobeRP,rundir=rundir,
                                    startyr=condC$yearCE[1],console=FALSE,fnt=7)
  }
  addtable(kobedata,"kobe_plot_endpoints.csv",rundir,category="phaseplot",
           caption=paste0("The median depletion and Harvest rate, and median ",
                          "Target CE and Grad4 values."))
  if (verbose) {
    incrtime2 <- Sys.time(); timeinc <- incrtime2 - incrtime1
    cat("All plots and tables completed ",timeinc,attr(timeinc,"units"),"\n")
  }
  # out list generation-------------------------------------------------------
  if (!includeNAS) NAS=NULL
  out <- list(tottime=tottime,runtime=projtime,starttime=starttime,glb=glb,
              ctrl=ctrl,zoneCP=zoneCP,zoneD=zoneD,zoneDD=zoneDD,zoneDP=zoneDP,
              NAS=NAS,projC=projC,condC=condC,sauout=sauout,outzone=outzone,
              production=production,condout=condout,saucomp=saucomp,
              HSstats=HSstats,saudat=saudat,constants=constants,hsargs=hsargs,
              sauprod=sauprod,zonesummary=zonesummary,
              kobedata=kobedata,outhcr=outhcr,scoremed=scoremed,
              popmedcatch=popmedcatch,popmedcpue=popmedcpue,
              popmeddepleB=popmeddepleB,pops=pops,fisproj=outproj$fisproj)
  return(out)
} # end of do_MSE








