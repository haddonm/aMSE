
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

#' @title do_condition is a function to assist with conditioning the model
#'
#' @description do_condition is a utility function that can be used to speed the
#'     conditioning of the operating model on available fishery data. During the
#'     generation of the simulated zone the operation that takes the longest is
#'     the estimation of the productivity of each population. The productivity
#'     is estimated at equilibrium and getting there is what tales the time.
#'     When conditioning the model the objective is to match the predicted
#'     fishery response to the historical catches to the observed responses.
#'     'do_condition' generates the simulated zone much more rapidly, which
#'     simplifies any searches for optimal values of AvRec (average unfished
#'     recruitment), and especially for suitable values of recruitment deviates.
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
#' @param calcpopC a function that takes the output from hcrfun (either the
#'     aspirational catch x SAU or TAC x zone) and generates the actual catch
#'     per population in each SAU expected in the current year.
#' @param cleanslate should the directory be emptied of all files first? All
#'     html, png, RData, and css files in the directory will be deleted.
#'     default = FALSE. This is obviously a very powerful and potentially
#'     dangerous argument, hence it needs to be set =TRUE explicitly. It does
#'     not delete any .csv files so if the rundir is used to store the data for
#'     the run then 'cleanslate' will not affect the data or any .R files.
#' @param verbose Should progress comments be printed to console, default=FALSE
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
do_condition <- function(rundir,controlfile,datadir,calcpopC,cleanslate=FALSE,
                         verbose=FALSE) {
  # generate equilibrium zone ----------------------------------------------------
  starttime <- Sys.time()
  zone <- makeequilzone(rundir,controlfile,datadir,doproduct=FALSE,
                        cleanslate=cleanslate,verbose=verbose)
  equiltime <- (Sys.time()); if (verbose) print(equiltime - starttime)
  # declare main objects -------------------------------------------------------
  glb <- zone$glb
  ctrl <- zone$ctrl
  zone1 <- zone$zone1
  projC <- zone$zone1$projC
  condC <- zone$zone1$condC
  zoneC <- zone$zoneC
  zoneD <- zone$zoneD
  # save some equil results ----------------------------------------------------
 # biology_plots(rundir, glb, zoneC)
 #numbersatsize(rundir, glb, zoneD)
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
  condtime <- Sys.time()
  tottime <- round((condtime - starttime),3)
  times <- list(tottime=tottime,condtime=condtime,starttime=starttime)
  out <- list(times=times,glb=glb,ctrl=ctrl,zoneC=zoneC,zoneDD=zoneDD,
              condC=condC,condout=condout)
  return(out)
} # end of do_MSE


