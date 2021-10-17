
#' @title changecolumn alters a selected column of recdevs in the control file
#'
#' @description changecolumn is designed to be used when conditioning the OM
#'     when optimizing ad hoc recruitment deviates by calculating the SSQ
#'     between the predicted and observed CPUE for a selected range of years.
#'     If the length of new values is not the same as the selected linerange
#'     the function will stop with a warning. When a line is changed then all
#'     -1 values are changed to +1
#'
#' @param rundir the rundir for the scenario
#' @param filename the character name of the control file
#' @param linerange the linerange within the control file containing the rows
#'     representing the selected years of recruitment deviates. Obtained by
#'     using 'getrecdevcolumn
#' @param column the sau + 1 to account for the initial year in each row of
#'     deviates
#' @param newvect the new values to replace those present. If using optim these
#'     would be exp(ans$par), if 'nlm' exp(ans$estimate)
#' @param verbose should console reports of before and after be made?
#'     default = FALSE
#'
#' @return nothing although values within the control file will be changed and,
#'     if verbose=TRUE, it will write to the console.
#' @export
#'
#' @examples
#' print("wait on suitable internal data sets")
changecolumn <- function(rundir, filename, linerange, column,newvect,
                         verbose=FALSE) {
  num <- length(newvect)
  if (num != length(linerange))
    stop("length of line-range differs from number of parameters")
  filen <- filenametopath(rundir,filename)
  dat <- readLines(filen)
  if (verbose) {
    print(dat[linerange])
    cat("\n\n")
  }
  for (linen in 1:num) { # linen=linerange[1]
    origtext <- dat[linerange[linen]]
    txt <- unlist(strsplit(origtext,","))
    txt[column] <- as.character(newvect[linen])
    pickt <- which(txt == "-1")
    if (length(pickt) > 0) txt[pickt] <- "1"
    outtxt <- paste0(txt,collapse=",")
    dat[linerange[linen]] <- outtxt
  }
  writeLines(dat,filen)
  if (verbose) print(dat[linerange])
} # end of changecolumn

#' @title changeline replaces a given line in a given file with new text
#'
#' @description changeline enables a text file to be changed line by line.
#'     One identifies a given line, perhaps by using 'findlinenumber', and
#'     this can then be replaced by an input character string. Obviously this
#'     is a fine way to mess up a data file so use with care.
#'
#' @param indir the directory path in which to find the text file. Usually,
#'     rundir or datadir
#' @param filename the full name of the text file in quotations.
#' @param linenumber either the line number within the text file to be changed,
#'     or, the character name of the variable to be changed, e.g. 'AvRec'
#' @param newline the character string with which to replace the line
#' @param verbose should confirmation be output to the console. default=FALSE
#'
#' @seealso{
#'  \link{findlinenumber}, \link{changevar}
#' }
#'
#' @return nothing but it does alter a line in a text file. Optionally it may
#'     confirm the action to the console
#' @export
#'
#' @examples
#' print("wait on an example")
changeline <- function(indir, filename, linenumber, newline,verbose=FALSE) {
  # indir=datadir; filename="saudataM15h5.csv"; linenumber="AvRec"; newline=chgevals
  filen <- filenametopath(indir,filename)
  dat <- readLines(filen)
  if (class(linenumber) == "numeric") {
    origtext <- dat[linenumber]
    dat[linenumber] <- newline
    writeLines(dat,filen)
  }
  if (class(linenumber) == "character") { # ie the name of the variable
    pickL <- grep(linenumber,dat)[1]
    if (length(pickL) == 0)
      stop(cat(paste0(linenumber," does not appear in ",filename),"\n"))
    origtext <- dat[pickL]
    dat[pickL] <- newline
    writeLines(dat,filen)
  }
  if (verbose) {
    cat(origtext," \n")
    cat("replaced with ",newline)
  }
} # end of changeline


#
# rundir=rundir
# controlfile="controltest_SA_trim.csv"
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
# verbose=TRUE
# doproduct=FALSE
# ndiagprojs=4
# makehcrout=makeouthcr

#' @title do_condition is a function to assist with conditioning the model
#'
#' @description do_condition is a utility function that can be used to speed the
#'     conditioning of the operating model on available fishery data. During the
#'     generation of the simulated zone the operation that takes the longest is
#'     the estimation of the productivity of each population. The productivity
#'     is estimated at equilibrium and getting there is what takes the time.
#'     When conditioning the model the objective is to match the predicted
#'     fishery response to the historical catches to the observed responses.
#'     'do_condition' generates the simulated zone without running the
#'     projections, which simplifies any searches for optimal values of AvRec
#'     (average unfished recruitment), and especially for suitable values of
#'     recruitment deviates. If any of the initdepl (initial depletion) values,
#'     as listed in the control file are less than 1.0 then before applying
#'     the historical catches it first depletes each population within each
#'     SAU to the initdepl value for each SAU. For this to work doproduct
#'     must be TRUE. This will slow down any searches for an optimum set of
#'     AvRec and recdevs, so a different strategy for doing that will be
#'     needed.
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
#' @param verbose Should progress comments be printed to console, default=FALSE
#' @param doproduct should production estimates be made. default=FALSE
#' @param dohistoric should the historical catches be applied. Default=TRUE
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
#' # rundir=rundir; controlfile=controlfile;datadir=datadir;calcpopC=calcexpectpopC
#' #verbose=TRUE; doproduct=TRUE; dohistoric=FALSE
do_condition <- function(rundir,controlfile,datadir,calcpopC,
                         verbose=FALSE,doproduct=FALSE,dohistoric=TRUE) {
  starttime <- Sys.time()
  zone <- makeequilzone(rundir,controlfile,datadir,doproduct=doproduct,
                        verbose=verbose)
  equiltime <- (Sys.time()); if (verbose) print(equiltime - starttime)
  # declare main objects
  glb <- zone$glb
  ctrl <- zone$ctrl
  zone1 <- zone$zone1
  projC <- zone$zone1$projC
  condC <- zone$zone1$condC
  zoneC <- zone$zoneC
  zoneD <- zone$zoneD
  production <- NULL
  if (doproduct) production <- zone$product
  #Condition on Fishery
  if (any(condC$initdepl < 1)) {
    if (!doproduct) stop("doproduct must be TRUE for Initial Depletion to work \n")
    initdepl <- condC$initdepl
    if (verbose) cat("Conducting initial depletions  ",initdepl,"\n")
    zoneD <- depleteSAU(zoneC,zoneD,glb,initdepl=initdepl,production)
  }
  if (dohistoric) {
    if (verbose) cat("Conditioning on the Fishery data  \n")
    zoneDD <- dohistoricC(zoneD,zoneC,glob=glb,condC,calcpopC=calcpopC,
                          sigR=1e-08,sigB=1e-08)
  } else {
    zoneDD <- zoneD
  }
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
  if (dohistoric) {
    condout <- plotconditioning(zoneDD,glb,zoneC,condC$histCE,rundir)
  } else {
    condout <- NULL
  }
  condtime <- Sys.time()
  tottime <- round((condtime - starttime),3)
  times <- list(tottime=tottime,condtime=condtime,starttime=starttime)
  out <- list(times=times,glb=glb,ctrl=ctrl,zoneC=zoneC,zoneD=zoneD,zoneDD=zoneDD,
              condC=condC,condout=condout,projC=projC,production=production)
  return(out)
} # end of do_condition

#' @title gettasdevssq calculates the SSQ between observed and predicted CPUE
#'
#' @description gettasdevssq is used when conditioning the aMSE operating model
#'     using ad hoc recruitment deviates after conditioning on AvRec. It focuses
#'     on a single SAU at a time. log transformed recdevs are used to avoid the
#'     possibility of negative recdevs (a no-no).
#'
#' @param param the log transformed sequence of recdevs for the selected years
#' @param rundir the rundir for the scenario
#' @param datadir the datadir for the scenario
#' @param ctrlfile the character name of the control file
#' @param calcpopC the function used to distribute catches among populations.
#'     This is from the external JurisdictionHS.R file that is 'sourced' in.
#' @param locyrs the rows within the matrix of recruitment deviates
#'     corresponding to the selected years.
#' @param sau which sau (in TAS 1 - 8) is to be worked on
#' @param verbose should console reports be made? default = FALSE
#'
#' @return a scalar value which is the SSQ for the selected years of CPUE
#' @export
#'
#' @examples
#' print("wait on suitable internal data sets")
gettasdevssq <- function(param,rundir,datadir,ctrlfile,calcpopC,locyrs,sau,
                         verbose=FALSE) {
  zone1 <- readctrlfile(rundir,infile=ctrlfile,datadir=datadir,verbose=verbose)
  # pickX <- which(param > 1.0)
  # if (length(pickX) > 0) param[pickX] <- 1.0
  zone1$condC$recdevs[locyrs,sau] <- exp(param)
  ctrl <- zone1$ctrl
  glb <- zone1$globals     # glb without the movement matrix
  bysau <- ctrl$bysau
  if (is.null(bysau)) bysau <- 0
  if (bysau) {
    constants <- readsaudatafile(datadir,ctrl$datafile)
  } else {
    constants <- readdatafile(glb$numpop,datadir,ctrl$datafile)
  }
  out <- setupzone(constants,zone1,doproduct=FALSE,verbose=verbose) # make operating model
  zoneC <- out$zoneC
  zoneD <- out$zoneD
  glb <- out$glb             # glb now has the movement matrix
  zone1$globals <- glb
  # did the larval dispersal level disturb the equilibrium?
  zoneD <- testequil(zoneC,zoneD,glb,verbose=verbose)
  zoneC <- resetexB0(zoneC,zoneD) # rescale exploitB to avexplB after dynamics
  setuphtml(rundir)
  zone <- list(zoneC=zoneC,zoneD=zoneD,glb=glb,constants=constants,
               product=NULL,ctrl=ctrl,zone1=zone1)
  projC <- zone1$projC
  condC <- zone1$condC
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
  return(ssq[sau])
} # end of gettasdevssq


#' @title getrecdevcolumn extracts a column of recdevs from a control file
#'
#' @description getrecdevcolumn is used when conditioning the operating model
#'     by modifying the recruitment deviates. getrecdevcolumn extracts a column
#'     of years selected from a selected controlfile. It does this to use those
#'     values as the starting parameter vector for optim.
#'
#' @param rundir the rundir for the scenario
#' @param filename the character name of the control file
#' @param yearrange the range of years of rec devs to be selected
#' @param sau which SAU or block is to be fitted?
#' @param verbose should details of the run be made. default=FALSE
#'
#' @return a list of the vewctor of recdevs, the vector of the yearrange, and
#'     a vector of the linenumbers within the control file that are to be
#'     changed
#' @export
#'
#' @examples
#' print("wait on suitable internal data-sets")
getrecdevcolumn <- function(rundir, filename, yearrange, sau,verbose=FALSE) {
  column <- sau+1
  nyr <- length(yearrange)
  recdevs <- numeric(nyr)
  filen <- filenametopath(rundir,filename)
  dat <- readLines(filen)
  nlines <- length(dat)
  begin <- grep("RECDEV",dat) + 2
  reclines <- length(begin:nlines)
  years <- matrix(0,nrow=reclines,ncol=2)
  years[,1] <- begin:nlines
  for (i in begin:nlines) {
    txt <- unlist(strsplit(dat[i],","))
    years[(i-begin+1),2] <- as.numeric(txt[1])
  }
  picky <- match(yearrange,years[,2])
  line1 <- years[picky[1],1]
  line2 <- years[tail(picky,1),1]
  for (i in line1:line2) { # i = line1
    txt <- unlist(strsplit(dat[i],","))
    x <- as.numeric(txt[column])
    if (x < 0) x <- 1
    recdevs[i-line1+1] <- x
  }
  return(list(recdevs=recdevs,yearrange=yearrange,linerange=line1:line2))
} # end of getrecdevcolumn




#' @title sauavrecssq applies the AvRec and returns the ssq from the cpue
#'
#' @description sauavrecssq applies the AvRec vector while adjusting just
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
#' @param linenum the linenumber containing the AvRec vector
#' @param calcpopC the HS function that calculates the aspirational catches
#'     for each SAU
#' @param extra the vector of SAU AvRecs minus the selected SAU value
#' @param picksau which SAU is to be worked on as in 1:nsau
#' @param nsau the total number of SAU
#'
#' @return the SSQ for the selected SAU in a comparison of predicted and
#'     observed CPUE
#' @export
#'
#' @examples
#' print("Wait on suitable internal data files")
sauavrecssq <- function(param,rundir,datadir,controlfile,datafile,linenum,
                        calcpopC,extra,picksau,nsau) {
  nsaum1 <- nsau - 1
  if (picksau == 1)
    replacetxt <- paste0("AvRec,",as.character(param),",",
                         paste0(as.character(extra[1:nsaum1]),collapse=","),
                         collapse=",")
  if ((picksau > 1) & (picksau < nsau))
    replacetxt <- paste0("AvRec,",
                         paste0(as.character(extra[1:(picksau-1)]),collapse=","),
                         ",",as.character(param),",",
                         paste0(as.character(extra[picksau:nsaum1]),collapse=","))
  if (picksau == nsau)
    replacetxt <- paste0("AvRec,",
                         paste0(as.character(extra[1:nsaum1]),collapse=","),
                         ",",as.character(param),collapse=",")
  changeline(datadir,datafile,linenum,replacetxt)
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
} # end of sauavrecssq
