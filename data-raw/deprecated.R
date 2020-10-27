
#' @title checkctrldat checks datadir contains the required csv files
#'
#' @description checkctrldat checks datadir contains the required csv
#'     files including the named control file, which then contains
#'     the names of the region data file, and the population data
#'     file. The run stops if any are not present or are misnamed.
#'
#' @param datadir the directory in which all all files relating to a
#'     particular run are to be held.
#' @param ctrlfile default="control.csv", the filename of the control
#'     file present in datadir containing information regarding the
#'     run.
#'
#' @return the control list for the run
#' @export
#'
#' @examples
#' resdir <- tempdir()
#' ctrlfiletemplate(resdir)
#' zonefiletemplate(resdir)
#' datafiletemplate(6,resdir,filename="zone1sau2pop6.csv")
#' ctrl <- checkctrldat(resdir)
#' ctrl
checkctrldat <- function(datadir,ctrlfile="control.csv") { # resdir=resdir; ctrlfile="control.csv"
  filenames <- dir(datadir)
  if (length(grep(ctrlfile,filenames)) != 1)
    stop(cat(ctrlfile," not found in resdir \n"))
  ctrol <- readctrlfile(datadir,infile=ctrlfile)
  if (length(grep(ctrol$zonefile,filenames)) != 1)
    stop("zone data file not found \n")
  if (length(grep(ctrol$datafile,filenames)) != 1)
    stop("population data file not found \n")
  cat("All required files appear to be present \n")
  return(ctrol)
} # end of checctrldat

#' @title ctrlfiletemplate generates a template input control file
#'
#' @description ctrlfiletemplate generates a standardized control file
#'     template. Generate this and then modify the contents to suit
#'     the system you are attempting to simulate. Defaults to 10
#'     replicates, and with
#'     effectively no variation in recruitment, mature biomass or cpue
#'     calculations.
#'
#' @param indir directory in which to place the control.csv file
#' @param filename the name for the generated ctrlfile, a character
#'     string that defaults to control.csv.
#'
#' @return invisibly the fill path and name of the control file. More
#'     importantly, it write a control file template to that directory.
#'
#' @export
#'
#' @examples
#'  yourdir <- tempdir()
#'  ctrlfiletemplate(yourdir,filename="testctrl.csv")   #
#'  control <- readctrlfile(yourdir,"testctrl.csv")
#'  str(control,max.level=1)
ctrlfiletemplate <- function(indir,filename="control.csv") {
  filename <- filenametopath(indir,filename)
  cat("Control file containing details of a particular run \n",
      file=filename,append=FALSE)
  cat("Modify the contents to suit your own situation \n",
      file=filename,append=TRUE)
  cat("\n",file=filename,append=TRUE)
  cat("START \n",file=filename,append=TRUE)
  cat("runlabel, testrun, label for particular run \n",
      file=filename,append=TRUE)
  cat("zonefile, zone1.csv, name of zone wide constants \n",
      file=filename,append=TRUE)
  cat("datafile, zone1sau2pop6.csv, name of popdefs file \n",
      file=filename,append=TRUE)
  cat("\n",file=filename,append=TRUE)
  cat("\n",file=filename,append=TRUE)
  cat("zoneCOAST \n",file=filename,append=TRUE)
  cat("batch,  0, run job interactively or as a batch (0 = FALSE) \n",
      file=filename, append=TRUE)
  cat("replicates,  10, number of replicates, usually 1000  \n",
      file=filename, append=TRUE)
  cat("withsigR,  1.0e-08, recruitment variability eg 0.3 \n",
      file=filename, append=TRUE)
  cat("withsigB,  1.0e-08, process error on mature biomass \n",
      file=filename, append=TRUE)
  cat("withsigCE, 1.0e-08, process error on cpue calculations  \n",
      file=filename, append=TRUE)
  return(invisible(filename))
} # end of ctrlfiletemplate

#' @title readctrlfile reads a csv file for controlling aMSE
#'
#' @description readctrlfile complements the readdatafile. The MSE requires
#'     a data file to condition the operating model but it also needs a
#'     control file to setup the details of the simulation test being
#'     conducted. See ctrlfileTemplate./
#'
#' @param indir directory in which to find the control file
#' @param infile filename of the control file, default="control.csv"
#'
#' @return a list object containing the control variables
#' @export
#'
#' @examples
#' \dontrun{
#'   direct="./../../rcode2/aMSE/data-raw"
#'   infile="control.csv"
#'   ctrl <- readctrlfile(indir=direct,infile=infile)
#'   str(ctrl)
#'  }
readctrlfile <- function(indir,infile="control.csv") {
  filename <- filenametopath(indir,infile)
  indat <- readLines(filename)   # reads the whole file as character strings
  begin <- grep("START",indat) + 1
  runlabel <- getStr(indat[begin],1)
  zonefile <- getStr(indat[begin+1],1)
  datafile <- getStr(indat[begin+2],1)
  batch <- getsingleNum("batch",indat)
  reps <- getsingleNum("replicates",indat)
  #  initdepl <- getsingleNum("initdepl",indat)
  withsigR <- getsingleNum("withsigR",indat)
  withsigB <- getsingleNum("withsigB",indat)
  withsigCE <- getsingleNum("withsigCE",indat)
  outctrl <- list(runlabel,zonefile,datafile,batch,
                  reps,withsigR,withsigB,withsigCE)
  names(outctrl) <- c("runlabel","zonefile","datafile","batch",
                      "reps","withsigR","withsigB","withsigCE")
  return(outctrl)
} # end of readctrlfile

#' @title readzonefile reads in the constants for the zone
#'
#' @description with the zone filename from the control file the
#'     readzonefile will read the data from a csv file arranged with
#'     a standard layout. Once again this uses the utility functions
#'     for reading and parsing lines of text. This is illustrated by
#'     the function makezonefile, which produces an example *.csv
#'     file that can then be customized to suit your own simulations.
#'     Each required section contains a series of constants which are
#'     read in individually,so their labels are equally important. This
#'     important function reads in the details of any conditioning and also
#'     of any projections.
#'
#' @param indir directory in which to find the zone file
#' @param infile character string with filename of the zone file
#'
#' @return 15 objects relating to constnats for the zone
#' \itemize{
#'   \item SAUnames the labels given to each SAU
#'   \item SAUpop the number of populations in each SAU, in sequence
#'   \item minc scaler defining the mid point of the smallest size class
#'   \item cw scaler defining the class width
#'   \item randomseed used to ensure that each simulation run starts
#'       in the same place. Could use getseed() to produce this.
#'   \item projyrs, he number of years of projection
#'   \item outyear a vector of three, with projyrs, fixyear, and firstyear
#'   \item projLML the LML expected in each projection year, a vector
#'   \item HS the name of the harvest strategy to be used in the projection
#'   \item condition the number of years of historical catch data
#'   \item histCatch the years if historical catch and historical LML
#'   \item histyr the actual years for which SAU catch data available
#'   \item histCE the historical cpue data for each SAU
#'   \item yearCE the actual years fo which cpue data available
#'   \item globals a list containing numpop, nSAU, midpts, Nclass, and
#'       Nyrs
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' datadir <- "./../../rcode2/aMSE/data-raw/"
#' ctrlfile <- "control.csv"
#' ctrl <- readctrlfile(datadir,ctrlfile)
#' reg1 <- readzonefile(datadir,ctrl$zonefile)
#' str(reg1)
#' }
readzonefile <- function(indir,infile) {  # infile=ctrl$zonefile; indir=resdir
  context <- "zone_file"
  filename <- filenametopath(indir,infile)
  indat <- readLines(filename)   # reads the whole file as character strings
  nSAU <-  getsingleNum("nSAU",indat) # number of spatial management units
  begin <- grep("SAUpop",indat)
  SAUpop <-  getConst(indat[begin],nSAU) # how many populations per SAU
  numpop <- sum(SAUpop)
  SAUnames <- getStr(indat[begin+1],nSAU)
  minc <-  getsingleNum("minc",indat) # minimum size class
  cw    <- getsingleNum("cw",indat) # class width
  Nclass <- getsingleNum("Nclass",indat) # number of classes
  midpts <- seq(minc,minc+((Nclass-1)*cw),2)
  larvdisp <- getsingleNum("larvdisp",indat)
  randomseed <- getsingleNum("randomseed",indat)
  initLML <- getsingleNum("initLML",indat)
  projyrs <- getsingleNum("PROJECT",indat)
  firstyear <- getsingleNum("firstyear",indat)
  fixyear <- getsingleNum("fixyear",indat)
  outyear <- c(projyrs,fixyear,firstyear)
  projLML <- NULL
  HS <- NULL
  if (projyrs > 0) {
    projLML <- numeric(projyrs)
    from <- grep("PROJLML",indat)
    for (i in 1:projyrs) {
      from <- from + 1
      projLML[i] <- getConst(indat[from],1)
    }
    start <- grep("HARVESTS",indat)
    HS <- getStr(indat[start],1)
  }
  condition <- getsingleNum("CONDITION",indat)
  if (condition > 0) {
    Nyrs <- condition  # don't forget to add an extra year for initiation
    begin <- grep("CondYears",indat)
    histCatch <- matrix(0,nrow=condition,ncol=nSAU)
    colnames(histCatch) <- SAUnames
    histyr <- matrix(0,nrow=condition,ncol=2)
    colnames(histyr) <- c("year","histLML")
    for (i in 1:condition) {
      begin <- begin + 1
      asnum <- as.numeric(unlist(strsplit(indat[begin],",")))
      histyr[i,] <- asnum[1:2]
      histCatch[i,] <- asnum[3:(nSAU+2)]
    }
    rownames(histCatch) <- histyr[,1]
    rownames(histyr) <- histyr[,1]
    yrce <- getsingleNum("CEYRS",indat)
    begin <- grep("CPUE",indat)
    histCE <- matrix(NA,nrow=yrce,ncol=nSAU)
    yearCE <- numeric(yrce) # of same length as nSAU
    colnames(histCE) <- SAUnames
    for (i in 1:yrce) {
      begin <- begin + 1
      cenum <- as.numeric(unlist(strsplit(indat[begin],",")))
      yearCE[i] <- cenum[1]
      histCE[i,] <- cenum[2:(nSAU+1)]
    }
    rownames(histCE) <- yearCE
    sizecomp <- getsingleNum("SIZECOMP",indat)
    if (sizecomp > 0) {
      lffiles <- NULL
      locsizecomp <- grep("SIZECOMP",indat)
      for (i in 1:sizecomp)
        lffiles <- c(lffiles,getStr(indat[locsizecomp+i],1))
      compdat <- vector("list",sizecomp)
      for (i in 1:sizecomp) {
        filename <- filenametopath(indir,lffiles[i])
        compdat[[i]] <- read.csv(file=filename,header=TRUE)
      }
    }
  } else {
    histCatch <- NULL
    histyr <- NULL
    histCE <- NULL
    yearCE <- NULL
    compdat=NULL
  }
  if ((projyrs == 0) & (condition == 0)) Nyrs=40
  globals <- list(numpop=numpop, nSAU=nSAU, midpts=midpts,
                  Nclass=Nclass, Nyrs=Nyrs, larvdisp=larvdisp)
  totans <- list(SAUnames,SAUpop,minc,cw,larvdisp,randomseed,projyrs,outyear,
                 initLML,projLML,HS,condition,histCatch,histyr,histCE,yearCE,
                 compdat,globals)
  names(totans) <- c("SAUnames","SAUpop","minc","cw","larvdisp","randomseed",
                     "projyrs","outyear","initLML","projLML","HS",
                     "condition","histCatch","histyr","histCE","yearCE",
                     "compdat","globals")
  return(totans)
}  # end of readzonefile

#' @title zonefiletemplate generates a template zone file
#'
#' @description zonefiletemplate generates a standardized zone file
#'     template containing constants tha toperate at a zoneal level.
#'     Generate this and then modify the contents to suit
#'     the system you are attempting to simulate.
#'
#' @param indir directory in which to place the zone file
#' @param filename the name for the generated zone file, a character
#'     string that defaults to zone1.csv. It is best to include the
#'     complete path
#'
#' @return nothing, but it creates a template zone file within indir
#'
#' @export
#'
#' @examples
#'  yourdir <- tempdir()
#'  zonefiletemplate(yourdir,filename="zone1.csv")   #
#'  zone1 <- readzonefile(yourdir,"zone1.csv")
#'  str(zone1,max.level=1)
zonefiletemplate <- function(indir,filename="zone1.csv") {
  filename <- filenametopath(indir,filename)
  cat("zone file containing zone wide constants for a run \n",
      file=filename,append=FALSE)
  cat("Modify the contents to suit your own situation \n",
      file=filename,append=TRUE)
  cat("\n",file=filename,append=TRUE)
  cat("zone \n",file=filename,append=TRUE)
  cat("nSAU, 2, number of spatial management units eg 2 \n",
      file=filename,append=TRUE)
  cat("SAUpop, 2, 4, number of populations per SAU in sequence \n",
      file=filename,append=TRUE)
  cat("SAUnames, block1, block2, labels for each SAU \n",file=filename,
      append=TRUE)
  cat("\n",file=filename,append=TRUE)
  cat("SIZE \n",file=filename,append=TRUE)
  cat("minc, 2, centre of minimum size class \n",file=filename,append=TRUE)
  cat("cw, 2, class width mm \n",file=filename,append=TRUE)
  cat("Nclass, 105, number of size classes \n",file=filename,append=TRUE)
  cat("\n",file=filename,append=TRUE)
  cat("RECRUIT \n",file=filename,append=TRUE)
  cat("larvdisp, 0.03, rate of larval dispersal eg 0.03=3precent \n",
      file=filename,append=TRUE)
  cat("\n",file=filename,append=TRUE)
  cat("RANDOM \n",file=filename,append=TRUE)
  cat("randomseed, 4024136, for repeatability of results if wished \n",
      file=filename,append=TRUE)
  cat("\n",file=filename,append=TRUE)
  cat("initLML, 132, the initial LML for generating the unfished zone \n",
      file=filename,append=TRUE)
  cat("\n",file=filename,append=TRUE)
  cat("PROJECT, 0, number of projection years for each simulation \n",
      file=filename,append=TRUE)
  cat("firstyear, 2020, first year of simulation if 1 hypothetical 2014 conditioned \n",
      file=filename,append=TRUE)
  cat("fixyear, 3, conditions fixed for the first three years \n",
      file=filename,append=TRUE)
  cat("\n",file=filename,append=TRUE)
  cat("\n",file=filename,append=TRUE)
  cat("HARVESTS, MQMF, the name of the harvest strategy to apply \n",
      file=filename,append=TRUE)
  cat("\n",file=filename,append=TRUE)
  cat("PROJLML, same number of projection years \n",file=filename,
      append=TRUE) # ensure there are Nyrs lines
  cat("2020, 132,  Legal Minimum Length (LML, MLL, MLS) e.g. 140 \n",
      file=filename,append=TRUE)
  for (i in 2:10) {
    yr <- 2020 + i - 1
    cat(as.character(yr),", 132 \n",file=filename,append=TRUE)
  }
  cat("\n",file=filename,append=TRUE)
  cat("CONDITION, 0, if > 1 then how many years in the histLML \n",
      file=filename,append=TRUE)
  cat("\n",file=filename,append=TRUE)
  cat("HistoricalLML, 0, if >1 then number reflects number of histLML \n",
      file=filename,append=TRUE)
} # end of zonefiletemplate
