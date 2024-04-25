# Tas context -----------------
# postfixdir <- "BCmeta"
# rundir <- rundir
# controlfile=controlfile
# hsargs=hsargs
# hcrfun=mcdahcr  #mcdahcr   #  # consthcr   #constantrefhcr
# sampleCE=tasCPUE
# sampleFIS=tasFIS
# sampleNaS=tasNaS
# getdata=tasdata
# calcpopC=calcexpectpopC
# varyrs=7
# startyr=38
# verbose=TRUE
# ndiagprojs=4
# makeouthcr=makeouthcr   #makeoutconst - only used with consthcr# fleetdyn=NULL
# plotmultflags=plotmultandflags
# scoreplot = plotfinalscores
# cutcatchN=56
# matureL = c(70,200)
# wtatL = c(80,200)
# mincount=100
# includeNAS=FALSE
# depensate=0
# pmwtSwitch = 0
# stablewts = c(0.4, 0.5, 0.1)
# hcrname="mcdahcr"  #"constantrefhcr"  #"mcdahcr"     #    #"consthcr"
# kobeRP = c(0.4,0.2,0.15)
# interimout=""
# nasInterval=5
# minsizecomp=c(100,135)
# uplimH=0.4
# incH=0.005
#

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
# makeouthcr=makeouthcr
# hcrscoreoutputs = saextractandplotscores
# fleetdyn = safleetdyn
# scoreplot = plotfinalscores
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
#' @param minsizecomp a vector of two values, the first being the minimum size-
#'     class to be used when plotting the predicted size-composition of the
#'     stock, and the second being the minimum size-class when plotting the
#'     numbers-at-size in the catch. The default=c(100,135)
#' @param uplimH defines the upper limit of harvest rate used when estimating
#'     the productivity (also important when initial depletion is not 1.0). The
#'     default = 0.4
#' @param incH defines the interval between H steps when estimating productivity
#'     default = 0.005
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
                   minsizecomp=c(100,135),uplimH=0.4,incH=0.005) {
  # generate equilibrium zone -----------------------------------------------
  starttime <- (Sys.time())
  zone <- makeequilzone(rundir,controlfile,doproduct=TRUE,uplimH=uplimH,
                        incH=incH,verbose=verbose)
  if (verbose) {
    incrtime1 <- Sys.time(); timeinc <- incrtime1 - starttime
    cat("makeequilzone ",timeinc,attr(timeinc,"units"),"\n")
  }
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
  # biology, recruits, and production tabs -------------------------------------------------
  biology_plots(rundir, glb, zoneC, matL=matureL,Lwt=wtatL)
  sauprod <- plotproductivity(rundir,production,glb,hsargs)
 # hsargs$saumsy <- sauprod[3,]
  # numbers-at-size tab------------------------------------------------------
  numbersatsize(rundir, glb, zoneD)
  #zoneDD tab -----------------------------------------------------
  if (any(condC$initdepl < 1)) {
    initdepl <- condC$initdepl
    if (verbose) cat("Conducting initial depletions  ",initdepl,"\n")
    zoneD <- depleteSAU(zoneC,zoneD,glb,initdepl=initdepl,production)
  }
  if (verbose) cat("Conditioning on the Fishery data  \n")
  zoneDD <- dohistoricC(zoneD,zoneC,glob=glb,condC,calcpopC=calcpopC,
                        fleetdyn=fleetdyn,hsargs=hsargs,sigR=1e-08,sigB=1e-08)

  if (verbose) {
    incrtime2 <- Sys.time(); timeinc <- incrtime2 - incrtime1
    cat("dohistoricC ",timeinc,attr(timeinc,"units"),"\n")
  }
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
  # condition tab----------------------------------------
  condout <- plotconditioning(zoneDD,glb,zoneC,condC$histCE,
                              histCatch=condC$histCatch,rundir,
                              condC$recdevs,console=FALSE)
  if (verbose) {
    incrtime1 <- Sys.time(); timeinc <- incrtime1 - incrtime2
    cat("Conditioning plots completed ",timeinc,attr(timeinc,"units") ,"\n")
  }
  saurecdevs(condC$recdevs,glb,rundir,filen="saurecdevs.png")
  # predictedcatchN tab-----------------------------------------
  catchN <- zoneDD$catchN
  sauCt <- popNAStosau(catchN,glb)
  compdat <- condC$compdat$lfs
  if (!is.null(compdat)) {
    for (plotsau in 1:glb$nSAU) { # plotsau=1
      lfs <- preparesizecomp(compdat[,,plotsau],mincount=mincount)
      yrsize <- as.numeric(colnames(lfs))
      histyr <- condC$histyr
      pickyr <- match(yrsize,histyr[,"year"])
      LML <- histyr[pickyr,"histLML"]
      plotsizecomp(rundir=rundir,incomp=lfs,SAU=glb$saunames[plotsau],lml=LML,
                   catchN=sauCt[,,plotsau],start=NA,proportion=TRUE,
                   console=FALSE)
    }
    # OrigComp tab--------------------------------------------
    saucompdata(allcomp=compdat,glb=glb,horizline=5,console=FALSE,rundir=rundir,
                ylabel="Size-Composition of Catches",tabname="OrigComp")
  }
  # do projections ------------------------------------------------------------
  if (verbose) cat("Preparing for the projections \n")
  outpp <- prepareprojection(projC=projC,condC=condC,zoneC=zoneC,glb=glb,
                             calcpopC=calcpopC,zoneDD=zoneDD,
                             ctrl=ctrl,varyrs=varyrs,lastsigR = ctrl$withsigR)
  zoneDP <- outpp$zoneDP
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
                           makehcrout=makeouthcr,fleetdyn=fleetdyn,verbose=TRUE)
  if (verbose) {
    incrtime1 <- Sys.time(); timeinc <- incrtime1 - incrtime2
    cat("All projections finished ",timeinc,attr(timeinc,"units") ,"\n")
    cat("Now generating final plots and tables \n")
  }
  zoneDP=outproj$zoneDP
  hcrout <- outproj$hcrout
  outhcr <- outproj$outhcr
  NAS <- list(Nt=zoneDP$Nt,catchN=zoneDP$catchN)
  zoneDP <- zoneDP[-c(17,16,15)]  # This removes the Nt etc from zoneDP
  if (nchar(interimout) > 0) {
    confirmdir(interimout,ask=verbose)
    outfile <- pathtopath(interimout[1],paste0("temp",interimout[2],".Rdata"))
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
  # projSAU tab----------------------------------------------------------
  if (verbose) cat("Starting the sau related plots \n")
  sauout <- sauplots(zoneDP,NAS,glb,rundir,B0,ExB0,startyr=startyr,
                     addCI=TRUE,histCE=condC$histCE,hlines=sauprod)
  sauNAS <- list(Nt=sauout$Nt,catchN=sauout$catchN)
  sauout <- sauout[-c(12,11)]  # This removes the Nt and catchN from sauout
  if (verbose) {
    incrtime2 <- Sys.time(); timeinc <- incrtime2 - incrtime1
    cat("Finished all sau plots ",timeinc,attr(timeinc,"units"),"\n")
  }
  # Add to numbers-at-size tab------------------------------------------------------
  if (verbose) cat("Starting size-composition plots \n")
  for (sau in 1:glb$nSAU) # sau=1
    predsizecomp(sau=sau, NSC=sauNAS$Nt, glb=glb, minSL=minsizecomp[1],
                 interval=nasInterval,prop=TRUE,console=FALSE,rundir=rundir)
  for (sau in 1:glb$nSAU)
    prob <- predsizecomp(sau=sau, NSC=sauNAS$catchN, glb=glb, minSL=minsizecomp[2],
                         interval=nasInterval,prop=TRUE,console=FALSE,rundir=rundir)
  if ((verbose) & (prob != "OK")) cat(prob,"\n")
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
  outfish <- fishery_plots(rundir=rundir,glb=glb,select=zoneCP[[1]]$Select,
                           histyr=condC$histyr,projLML=projC$projLML)
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
  save_hsargs(rundir,hsargs)   # prints hsargs to HSPerfs tab
  if (verbose) cat("hsargs.txt saved to rundir \n")
  if (verbose) cat("plotting HS performance statistics \n")
  plothsstats(rundir,HSstats,glb,average=FALSE)
  plothsstats(rundir,HSstats,glb,average=TRUE)
  addtable(hcrout$refpts,"hcrout_refpts.csv",rundir,category="HSperf",
           caption="HCR reference points")
  # scores tab----------------------------------------------------
  nSAU <- glb$nSAU
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
  # poplevelplots tab ---------------------------------
  # generate population plots
  if (verbose) cat("plotting Population level dynamics \n")
  popmedcatch <- vector(mode="list",length=nSAU)
  names(popmedcatch) <- glb$saunames
  popmedcpue <- vector(mode="list",length=nSAU)
  names(popmedcpue) <- glb$saunames
  popmeddepleB <- vector(mode="list",length=nSAU)
  names(popmeddepleB) <- glb$saunames
  for (sau in 1:nSAU) {
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
              production=production,condout=condout,
              HSstats=HSstats,saudat=saudat,constants=constants,hsargs=hsargs,
              sauprod=sauprod,zonesummary=zonesummary,
              kobedata=kobedata,outhcr=outhcr,scoremed=scoremed,
              popmedcatch=popmedcatch,popmedcpue=popmedcpue,
              popmeddepleB=popmeddepleB)
  return(out)
} # end of do_MSE








