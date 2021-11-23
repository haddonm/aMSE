





#' @title saureccpuessq applies the AvRec and returns the ssq from the cpue
#'
#' @description saureccpuessq applies the AvRec vector while adjusting just
#'     one for the picksau, and then compares the predicted CPUE with the
#'     observed and returns a sum-of-squared differences. This is used by
#'     optim while using the Brent option to find an optimum for a single
#'     parameter.
#'
#' @param param is the AvRec for the selected SAU
#' @param rundir the rundir for the scenario being considered
#' @param datadir the directory holding the saudata file, usually = rundir
#' @param controlfile the name of the control file, no path
#' @param datafile the name of the data file, no path
#' @param linenum a vector of two linenumbers, the first for the AvRec vector,
#'     the second for the MaxCEpars vector
#' @param calcpopC the HS function that calculates the aspirational catches
#'     for each SAU
#' @param extrarec the vector of SAU AvRecs minus the selected SAU value
#' @param extrace the vector of SAU MaxCEpars minus the selected SAU value
#' @param picksau which SAU is to be worked on as in 1:nsau
#' @param nsau the total number of SAU
#'
#' @return the SSQ for the selected SAU in a comparison of predicted and
#'     observed CPUE
#' @export
#'
#' @examples
#' print("Wait on suitable internal data files")
saureccpuessq <- function(param,rundir,datadir,controlfile,datafile,linenum,
                        calcpopC,extrarec,extrace,picksau,nsau) {
#  param=param;rundir=rundir;datadir=rundir;controlfile=controlfile; datafile=datafile
#  linenum=linenum;calcpopC=calcpopC;extrarec=extrarec;extrace=extrace;picksau=sau;nsau=8

  nsaum1 <- nsau - 1
  if (picksau == 1) {
    replaceRec <- paste0("AvRec,",as.character(param[1]),",",
                         paste0(as.character(extrarec[1:nsaum1]),collapse=","),
                         collapse=",")
    replaceCE <- paste0("MaxCEpars,",as.character(param[2]),",",
                         paste0(as.character(extrace[1:nsaum1]),collapse=","),
                         collapse=",")
  }
  if ((picksau > 1) & (picksau < nsau)) {
    replaceRec <- paste0("AvRec,",
                         paste0(as.character(extrarec[1:(picksau-1)]),collapse=","),
                         ",",as.character(param[1]),",",
                         paste0(as.character(extrarec[picksau:nsaum1]),collapse=","))
    replaceCE <- paste0("MaxCEpars,",
                         paste0(as.character(extrace[1:(picksau-1)]),collapse=","),
                         ",",as.character(param[2]),",",
                         paste0(as.character(extrace[picksau:nsaum1]),collapse=","))
  }
  if (picksau == nsau) {
    replaceRec <- paste0("AvRec,",
                         paste0(as.character(extrarec[1:nsaum1]),collapse=","),
                         ",",as.character(param[1]),collapse=",")
    replaceCE <- paste0("MaxCEpars,",
                         paste0(as.character(extrace[1:nsaum1]),collapse=","),
                         ",",as.character(param[2]),collapse=",")
  }
  changeline(datadir,datafile,linenum[1],replaceRec)
  changeline(datadir,datafile,linenum[2],replaceCE)
  zone <- makeequilzone(rundir,controlfile,datadir,doproduct=FALSE,
                        verbose=FALSE)
  # declare main objects
  glb <- zone$glb
  condC <- zone$zone1$condC
  zoneC <- zone$zoneC
  zoneD <- zone$zoneD
  #Condition on Fishery
  zoneDD <- dohistoricC(zoneD,zoneC,glob=glb,condC,calcpopC=calcpopC,
                        sigR=1e-08,sigB=1e-08)
  hyrs <- glb$hyrs
  sauindex <- glb$sauindex
  popB0 <- getlistvar(zoneC,"B0")
  B0 <- tapply(popB0,sauindex,sum)
  popExB0 <- getlistvar(zoneC,"ExB0")
  ExB0 <- tapply(popExB0,sauindex,sum)
  sauZone <- getsauzone(zoneDD,glb,B0=B0,ExB0=ExB0)
  ssq <- compareCPUE(condC$histCE,sauZone$cpue,glb,rundir,filen="SAU_compareCPUE.png")
  return(ssq[picksau])
} # end of saureccpuessq





