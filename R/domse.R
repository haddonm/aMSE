#
# rundir=rundir
# controlfile="controlsau.csv"
# datadir=datadir
# hsargs=hsargs
# hcrfun=mcdahcr
# sampleCE=tasCPUE
# sampleFIS=tasFIS
# sampleNaS=tasNaS
# getdata=tasdata
# calcpopC=calcexpectpopC
# varyrs=7
# startyr=42
# cleanslate=FALSE
# verbose=TRUE
# doproject=doproject
# ndiagprojs=4

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
#'     always occur.
#'
#' @param rundir the full path to the directory in which all files relating to a
#'     particular run are to be held. If datadir != rundir then the main data
#'     files are kept in datadir. The juridictionHS.R can be held in either.
#' @param controlfile the filename of the control file present in rundir
#'     containing information regarding the particular run.
#' @param datadir the directory in which the SAU data file is to be found. This
#'     will usually be the rundir for the scenario run but if the data file is
#'     to be shared among a set of scenarios it will be more efficient to define
#'     a separate datadir
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
#' @param varyrs how many years at the end of the conditioning on the fishery,
#'     data into zoneDD, to which to add recruitment variation, default = 7,
#'     which suits the Tasmanian west coast. Used in prepareprojection
#' @param startyr the index for the first year of the conditioning dynamics to
#'     include in the plotted cpue trajectories, to give startyr:indexoflastyear.
#'     used in sauplots.
#' @param cleanslate should the directory be emptied of all files first? All
#'     html, png, RData, and css files in the directory will be deleted.
#'     default = FALSE. This is obviously a very powerful and potentially
#'     dangerous argument, hence it needs to be set =TRUE explicitly. It does
#'     not delete any .csv files so if the rundir is used to store the data for
#'     the run then 'cleanslate' will not affect the data or any .R files.
#' @param verbose Should progress comments be printed to console, default=FALSE
#' @param doproject should the projections be run. default = TRUE. When
#'     conditioning the operating model this should be set to FALSE to focus
#'     the program on the production of the completely conditioned zone and its
#'     outputs.
#' @param ndiagprojs the number of replicate trajectories to plot in the
#'     diagnostics tab to illustrate ndiagprojs trajectories to ensure that
#'     such projections appear realistic; default=3
#' @param openfile should the HTML pages be opened automatically? default=TRUE
#' @param savesauout should the sau dynamics object be saved as an sauoutD.RData
#'     file? 100 replicates of 56 populations for 58 years of conditioning and
#'     30 years of projection = about 5.7 Mb. default=FALSE.
#'
#' @seealso{
#'  \link{makeequilzone}, \link{dohistoricC}, \link{prepareprojection},
#'  \link{doprojections}
#' }
#'
#' @return a large list containing tottime, projtime, starttime, glb, ctrl,
#'     zoneDD, zoneDP, projC, condC, sauout, and outzone
#' @export
#'
#' @examples
#' print("wait on suitable data sets in data")
do_MSE <- function(rundir,controlfile,datadir,hsargs,hcrfun,sampleCE,sampleFIS,
                   sampleNaS,getdata,calcpopC,varyrs=7,startyr=42,
                   cleanslate=FALSE,verbose=FALSE,doproject=TRUE,ndiagprojs=3,
                   openfile=TRUE,savesauout=FALSE) {
  # generate equilibrium zone ----------------------------------------------------
  starttime <- (Sys.time())
  zone <- makeequilzone(rundir,controlfile,datadir,cleanslate=cleanslate,
                        verbose=verbose)
  equiltime <- (Sys.time()); if (verbose) print(equiltime - starttime)
  # declare main objects -------------------------------------------------------
  glb <- zone$glb
  ctrl <- zone$ctrl
  zone1 <- zone$zone1
  projC <- zone$zone1$projC
  condC <- zone$zone1$condC
  zoneC <- zone$zoneC
  zoneD <- zone$zoneD
  production <- zone$product
  # save some equil results ----------------------------------------------------
  biology_plots(rundir, glb, zoneC, matL=c(70,200))
  plotproductivity(rundir,production,glb)
  numbersatsize(rundir, glb, zoneD)
  #Condition on Fishery --------------------------------------------------------
  if (verbose) cat("Conditioning on the Fishery data  \n")
  zoneDD <- dohistoricC(zoneD,zoneC,glob=glb,condC,calcpopC=calcpopC,
                        sigR=1e-08,sigB=1e-08)
  hyrs <- glb$hyrs
  propD <- t(getzoneprops(zoneC,zoneDD,glb,year=hyrs))
  addtable(round(propD,4),"propertyDD.csv",rundir,category="zoneDD",caption=
             "Properties of zoneD after conditioning on historical catches.")
  addtable(round(t(zoneDD$harvestR[(hyrs-9):hyrs,]),4),"final_harvestR.csv",
           rundir,category="zoneDD",
           caption="Last ten years of population vs harvest rate.")
  popdefs <- getlistvar(zone$zoneC,"popdef")
  addtable(round(t(popdefs),3),"popdefs.csv",rundir,category="zoneDD",
           caption="Population vs Operating model parameter definitions")
  condout <- plotconditioning(zoneDD,glb,zoneC,condC$histCE,rundir)
  # do projections ------------------------------------------------------------
  if (doproject) {
    if (verbose) cat("Doing the projections \n")
    outpp <- prepareprojection(projC=projC,condC=condC,zoneC=zoneC,glb=glb,
                               calcpopC=calcpopC,zoneDD=zoneDD,
                               ctrl=ctrl,varyrs=varyrs,lastsigR = ctrl$withsigR)
    zoneDP <- outpp$zoneDP
    projC <- outpp$projC
    zoneCP <- outpp$zoneCP
   # Rprof()
    outproj <- doprojections(ctrl=ctrl,zoneDP=zoneDP,zoneCP=zoneCP,glb=glb,
                             hcrfun=hcrfun,hsargs=hsargs,sampleCE=sampleCE,
                             sampleFIS=sampleFIS,sampleNaS=sampleNaS,
                             getdata=getdata,calcpopC=calcpopC,verbose=TRUE)
  #  Rprof(NULL)
    if (verbose) cat("Now generating final plots and tables \n")
    zoneDP=outproj$zoneDP
    hcrout <- outproj$hcrout; #str(hcrout)
    NAS <- list(Nt=zoneDP$Nt,NumNe=zoneDP$NumNe,catchN=zoneDP$catchN)
    zoneDP <- zoneDP[-c(17,16,15)]
   # zoneDP <- zoneDP[-16]
   #  zoneDP <- zoneDP[-15]
    histCE <- condC$histCE
    B0 <- getvar(zoneC,"B0")
    ExB0 <- getvar(zoneC,"ExB0")
    sauout <- sauplots(zoneDP,NAS,glb,rundir,B0,ExB0,
                       startyr=startyr,addCI=TRUE,histCE=histCE)
    diagnosticsproj(sauout$zonePsau,glb,rundir,nrep=ndiagprojs)
    outzone <- poptozone(zoneDP,NAS,glb,
                         B0=sum(getvar(zoneC,"B0")),
                         ExB0=sum(getvar(zoneC,"ExB0")))
    plotZone(outzone,rundir,glb,startyr=startyr,CIprobs=c(0.05,0.5,0.95),
             addfile=TRUE)
    fishery_plots(rundir=rundir,glb=glb,select=zoneCP[[1]]$Select,
                  histyr=condC$histyr,projLML=projC$projLML)
  } else { # in case doproject = FALSE
    zoneDP <- NULL
    zoneCP <- NULL
    NAS <- NULL
    sauout <- NULL
    outzone <- NULL
    outpp=NULL
    hcrout=NULL
  }
  projtime <- Sys.time()
  tottime <- round((projtime - starttime),3)
  replist <- list(starttime=as.character(starttime),
                  endtime=as.character(projtime))
  if (savesauout) {
    sauoutD <- sauout[-c(12,11,10)]
    save(sauoutD,file=paste0(rundir,"/sauoutD.RData"))
  }
  projy <- ifelse(doproject,glb$pyrs,0)
  runnotes <- paste0(ctrl$runlabel,":  RunTime = ",tottime,
                     "  replicates = ",glb$reps,",   years projected = ",projy,
                     "  Populations = ",glb$numpop," and SAU = ",glb$nSAU,
                     "  Randomseed for conditioning = ",ctrl$randseed)

  make_html(replist = replist,  rundir = rundir,  datadir=datadir,
            controlfile=controlfile, datafile=ctrl$datafile, width = 500,
            openfile = openfile,  runnotes = runnotes,   verbose = FALSE,
            packagename = "aMSE",  htmlname = postdir)

  out <- list(tottime=tottime,projtime=projtime,starttime=starttime,glb=glb,
              ctrl=ctrl,zoneCP=zoneCP,zoneDD=zoneDD,zoneDP=zoneDP,NAS=NAS,
              projC=projC,condC=condC,sauout=sauout,outzone=outzone,
              hcrout=hcrout,production=production,condout=condout)
  return(out)
} # end of do_MSE
