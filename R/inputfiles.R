
# rutilsMH::listFunctions("C:/Users/User/Dropbox/rcode2/aMSE/R/inputfiles.R")

#' @title ctrlfiletemplate generates a template input control file
#'
#' @description ctrlfiletemplate generates a standardized control file
#'     template. Generate this and then modify the contents to suit
#'     the system you are attempting to simulate. Defaults to 10
#'     replicates, starting at 100% depletion,i.e. at B0 and with
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
   cat("regionfile, region1.csv, name of region wide constants \n",
       file=filename,append=TRUE)
   cat("datafile, reg1smu2pop6.csv, name of popdefs file \n",
       file=filename,append=TRUE)
   cat("hcrfile, HCRfile.csv, HCR details file name \n",file=filename,
       append=TRUE)
   cat("outdir, C:/Users/user/Dropbox/rcode2/aMSEUse, output directory \n",
       file=filename,append=TRUE)
   cat("\n",file=filename,append=TRUE)
   cat("\n",file=filename,append=TRUE)
   cat("REGIONCOAST \n",file=filename,append=TRUE)
   cat("batch,  0, run job interactively or as a batch (0 = FALSE) \n",
       file=filename, append=TRUE)
   cat("replicates,  10, number of replicates, usually 1000  \n",
       file=filename, append=TRUE)
   cat("initdepl,  1.0, the initial depletion for each population \n",
       file=filename, append=TRUE)
   cat("withsigR,  1.0e-08, recruitment variability eg 0.3 \n",
       file=filename, append=TRUE)
   cat("withsigB,  1.0e-08, process error on mature biomass \n",
       file=filename, append=TRUE)
   cat("withsigCE, 1.0e-08, process error on cpue calculations  \n",
       file=filename, append=TRUE)
   return(invisible(filename))
} # end of ctrlfiletemplate

#' @title diagnostics generates an array of diagnostic plots
#'
#' @description diagnostics can generates an array of diagnostic plots,
#'     however, currently it only plots the selectivity
#'
#' @param regC the constants for the region simulated
#' @param regD the dynamic aspects of the region simulated
#' @param glob the general globals
#' @param plot plot the output or not; default = TRUE
#'
#' @return invisibly returns the unique values of LML
#' @export
#'
#' @examples
#' \dontrun{
#'   print("stiull to devise an example.")
#' }
diagnostics <- function(regC,regD,glob,plot=TRUE) {   # inzone <- testzone
   useLML <- sapply(regC,"[[","LML")
   colLML <- apply(useLML,2,unique)
   valLML <- unique(colLML)
   nLML <- length(valLML)
   midpts <- glob$midpts
   Nclass <- glob$Nclass
   plotprep(width=6,height=4)
   plot(midpts,seq(0,1,length=Nclass),type="n",xlab="",
        ylab="",ylim=c(0,1.025),yaxs="i",panel.first=grid())
   for (i in 1:nLML) {
      pick <- which(colLML == valLML[i])
      lines(midpts,regC[[pick[1]]]$Select[,1],lwd=2,col=i)
   }
   title(ylab=list("Selectivity", cex=1.0, font=7),
         xlab=list("Shell Length (mm)", cex=1.0, font=7))
   return(invisible(valLML))
}  # end of diagnostics


#' @title datafiletemplate generates a template input datafile
#'
#' @description datafiletemplate generates a standard input datafile
#'     to use as a template, go in and edit it appropriately to suit
#'     your own needs. It contains the probability distributions that
#'     are sampled to provide the necessary biological constants for
#'     each population.
#'
#' @param numpop the totla number of populations in the region
#' @param indir the directory into which to place the data template
#' @param filename the name for the generated datafile, a character
#'     string, defaults to tmpdat.csv
#'
#' @return a standard definition data file, to be read by readdatafile
#'     whose name and path is returned invisibly
#' @export
#'
#' @examples
#' \dontrun{
#'  yourdir <- tempdir()
#'  datafiletemplate(numpop=6,yourdir,"tmpdat.csv")
#'  data(region1)
#'  glb <- region1$globals
#'  constants <- readdatafile(numpop,yourdir,"tmpdat.csv")
#'  str(constants,max.level=1)
#' }
datafiletemplate <- function(numpop,indir,filename="tmpdat.csv") {
   genconst <- function(invect) {
      nlab <- length(invect)
      invect <- as.character(invect)
      ans <- invect[1]
      if (nlab > 1) for (i in 2:nlab) ans <- paste(ans,invect[i],sep=",")
      return(ans)
   }
   filename <- filenametopath(indir,filename)
   cat("Population definitions containing the Probability Density Functions by parameter  popdefs \n",
       file=filename,append=FALSE)
   cat("Define the populations in the same sequence in which they occur along the coast \n",
       file=filename,append=TRUE)
   cat("\n",file=filename,append=TRUE)
   cat("PDFs,  32   \n",file=filename, append=TRUE)
   cat("popnum,1, 2, 3, 4, 5, 6, the index number for each population \n",
       file=filename, append=TRUE)
   cat("SMU, 1, 1, 2, 2, 2, 2, define the SMU index for each population \n",
       file=filename,append=TRUE)
   cat("MaxDL ,",genconst(rep(32,numpop)),", maximum growth increment \n",
       file=filename, append=TRUE)
   cat("sMaxDL ,",genconst(rep(2.5,numpop)),", variation of MaxDL \n",
       file=filename, append=TRUE)
   cat("L50 ,",genconst(rep(123.0,numpop)),", Length at 50% MaxDL \n",
       file=filename, append=TRUE)
   cat("sL50 ,",genconst(rep(5,numpop)),", variation of L50 \n",
       file=filename, append=TRUE)
   cat("L50inc ,",genconst(rep(45,numpop)),", L95 - L50 = delta =L50inc \n",
       file=filename, append=TRUE)
   cat("sL50inc ,",genconst(rep(1.25,numpop)),", variation of L50inc \n",
       file=filename, append=TRUE)
   cat("SigMax  ,",genconst(rep(4.6,numpop)),", max var around growth \n",
       file=filename, append=TRUE)
   cat("sSigMax ,",genconst(rep(0.1,numpop)),", var of SigMax \n",
       file=filename, append=TRUE)
   cat("LML ,",genconst(rep(140,numpop)),", current legal minimum length \n",
       file=filename, append=TRUE)
   cat("Wtb ,",genconst(rep(3.161963,numpop)),", weight-at-length exponent \n",
       file=filename, append=TRUE)
   cat("sWtb ,",genconst(rep(0.148,numpop)),", var of Wtb \n",
       file=filename, append=TRUE)
   cat("Wtbtoa ,",genconst(rep(962.8098,numpop)),
       ", intercept of power curve between Wtb and Wta \n",file=filename,
       append=TRUE)
   cat("sWtbtoa ,",genconst(rep(-14.3526,numpop)),
       ", exponent of power curve between Wtb and Wta \n",file=filename,
       append=TRUE)
   cat("Me ,",genconst(rep(0.2,numpop)),", Nat Mort, 0.2 maxage=23, 0.1 maxage=46 \n",
       file=filename, append=TRUE)
   cat("sMe ,",genconst(rep(0.003,numpop)),", var of Me\n",
       file=filename, append=TRUE)
   cat("AvRec ,",genconst(rep(13.0,numpop)),", log of unfished recruitment LnR0 \n",
       file=filename, append=TRUE)
   cat("sAvRec ,",genconst(rep(1,numpop)),", \n",file=filename, append=TRUE)
   cat("defsteep ,",genconst(rep(0.6,numpop)),", Beverton-Holt steepness \n",
       file=filename, append=TRUE)
   cat("sdefsteep ,",genconst(rep(0.02,numpop)),", \n",file=filename,
       append=TRUE)
   cat("L50C ,",genconst(rep(126.4222,numpop)),", length at 50% emergent \n",
       file=filename, append=TRUE)
   cat("sL50C ,",genconst(rep(0.5,numpop)),", \n",file=filename, append=TRUE)
   cat("deltaC ,",genconst(rep(10.0,numpop)),", length at 95% emergent \n",
       file=filename, append=TRUE)
   cat("sdeltaC ,",genconst(rep(0.1,numpop)),", \n",file=filename, append=TRUE)
   cat("MaxCEpar,",genconst(rep(0.37,numpop)),", max cpue t-hr \n",
       file=filename, append=TRUE)
   cat("MaxCEvar,",genconst(rep(0.02,numpop)),", \n",file=filename, append=TRUE)
   cat("selL50p ,",genconst(rep(0.0,numpop)),", L50 of selectivity \n",
       file=filename, append=TRUE)
   cat("selL95p ,",genconst(rep(1.5,numpop)),", L95 of selectivity \n",
       file=filename, append=TRUE)
   cat("SaMa,",genconst(rep(-16,numpop)),", matuirity logistic a par \n",
       file=filename, append=TRUE)
   cat("L50Mat,",genconst(rep(123.384,numpop)),", L50 for maturity b = -1/L50\n",
       file=filename, append=TRUE)
   cat("sL50Mat,",genconst(rep(4,numpop)),", \n",file=filename, append=TRUE)
   return(invisible(filename))
}  # end of datafileTemplate

# Utility functions used within parseFile, not exported


#' @title makefilename Generates a filename for output files
#'
#' @description makefilename Generates a filename for output files based upon
#'     he name of the HCR (hcrLabel), the runlabel, the number of reps and
#'     the initial depletion (initDepl), which should be sufficient to
#'     uniquely identify any different scenarios to be run
#'     Uses splitDate.
#'
#' @param hcrLabel the name of the harvest strategy beign used
#' @param runlabel the optional additional label for each run
#' @param reps the number of iterations of the MSE run of the dynamics
#' @param initDepl the initial depletion for the zone
#' @return  makefilename generates a filename for output files
#'
#' @export
#' @examples
#' \dontrun{
#'  tmp <- splitDate()
#'  print(tmp)
#'  print(names(tmp))
#'  print(as.numeric(tmp[1:3]))
#'  print("still need a full example for this function")
#' }
makefilename <- function(hcrLabel,runlabel,reps,initDepl) {
   dates <- splitDate()   # Devise filenames for this run from its details and date
   if (nchar(runlabel) == 0) {
      fileadd <- paste(hcrLabel,reps,initDepl,dates[5],sep="_")
   } else { fileadd <- paste(hcrLabel,reps,initDepl,dates[5],runlabel,sep="_")
   }
   return(fileadd)
}  # end of makefilename



#' @title makeLabel converts a vector of numbers or strings into a label
#'
#' @description makeLabel converts a vector of numbers or strings into a
#'     single label
#'
#' @param invect the input vector of numbers or strings
#' @param insep the seperator between each part of the label
#'
#' @return a text string
#' @export
#'
#' @examples
#' \dontrun{
#'  x <- 1:5
#'  makeLabel(x,insep="_")
#'  makeLabel(x,insep="-")
#' }
makeLabel <- function(invect,insep="_") {
   nlab <- length(invect)
   invect <- as.character(invect)
   ans <- invect[1]
   if (nlab > 1) for (i in 2:nlab) ans <- paste(ans,invect[i],sep=insep)
   return(ans)
}  # end of makeLabel

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
   regionfile <- getStr(indat[begin+1],1)
   datafile <- getStr(indat[begin+2],1)
   hcrfile <- getStr(indat[begin+3],1)
   outdir <- getStr(indat[begin+4],1)

   batch <- getsingleNum("batch",indat)
   reps <- getsingleNum("replicates",indat)
   initdepl <- getsingleNum("initdepl",indat)
   withsigR <- getsingleNum("withsigR",indat)
   withsigB <- getsingleNum("withsigB",indat)
   withsigCE <- getsingleNum("withsigCE",indat)
   outctrl <- list(runlabel,regionfile,datafile,hcrfile,outdir,
                   batch,reps,initdepl,withsigR,withsigB,withsigCE)
   names(outctrl) <- c("runlabel","regionfile","datafile","hcrfile",
                       "outdir","batch","reps","initdepl",
                       "withsigR","withsigB","withsigCE")
   return(outctrl)
} # end of readctrlfile

#' @title readdatafile reads in a matrix of data defining each population
#'
#' @description readdatafile expects a matrix of probability density
#'     function definitions that are used to define the populaitons
#'     used in the simualtion. These constitute the definition of
#'     popdefs.
#'
#' @param numpop the total number of populations across the region
#' @param indir directory in which to find the date file
#' @param infile character string with filename of the data file
#'
#' @return a matrix of values defining the PDFs used to define the
#'     properties of each population. The contents of popdefs
#' @export
#'
#' @examples
#' \dontrun{
#' data(region1)
#' glb <- region1$globals
#' glb
#' data(constants)
#' constants
#' ctrlfile <- "control.csv"
#' ctrl <- readctrlfile(glb$numpop,datadir,ctrlfile)
#' reg1 <- readregionfile(datadir,ctrl$regionfile)
#' popdefs <- readdatafile(reg1$globals,datadir,ctrl$datafile)
#' print(popdefs)
#' }
readdatafile <- function(numpop,indir,infile) {  # indir=datadir;infile="reg1smu2pop6.csv";numpop=6
   filename <- filenametopath(indir,infile)
   indat <- readLines(filename)   # reads the whole file as character strings
   begin <- grep("PDFs",indat)
   npar <- getConst(indat[begin],1)
   rows <- c("popnum","SMU","DLMax","sMaxDL","L50","sL50","L50inc","sL50inc","SigMax",
             "sSigMax","LML","Wtb","sWtb","Wtbtoa","sWtbtoa","Me","sMe",
             "AvRec","sAvRec","defsteep","sdefsteep","L50C","sL50C",
             "deltaC","sdeltaC","MaxCEpars","sMaxCEpars","selL50p",
             "selL95p","SaMa","L50Mat","sL50Mat")
   ans <- matrix(0,nrow=length(rows),ncol=numpop)
   begin <- begin + 1
   for (i in 1:npar) {
      ans[i,] <- getConst(indat[begin],numpop)
      begin <- begin + 1
   } # completed filling ans matrix
   rownames(ans) <- rows
   colnames(ans) <- ans["popnum",]
   return(ans)
} # end of readdatafile

#' @title readhcrfile literally reads a csv file for controling the AbMSE
#'
#' @description readhcrfile complements the readdatafile. The MSE requires
#'     a data file to condition the operating model but it also needs a
#'     control file to setup the details of the simulation test being
#'     conducted. See ctrlfileTemplate./
#'
#' @param infile the filename of the control file
#'
#' @return a list object containing the contro variables
#' @export
#'
#' @examples
#' \dontrun{
#'   print("Still to be developed.")
#' }
readhcrfile <- function(infile) {  # infile <- "C:/A_CSIRO/Rcode/AbMSERun/ctrl_west.csv"
   indat <- readLines(infile)   # reads the whole file as character strings
   begin <- grep("batch",indat)
   batch <-  getLogical(indat[begin],1) # minimum size class
   reps <- getsingleNum("replicates",indat)
   initDepl <- getsingleNum("initDepl",indat)
   assessInterval <- getsingleNum("assessInterval",indat)
   recthreshold <- getsingleNum("recthreshold",indat)
   begin <- grep("hcrLabel",indat)
   hcrLabel <- getStr(indat[begin],1)
   begin <- grep("ConstC",indat)
   ConstC <- getLogical(indat[begin],1)
   begin <- grep("mcdaHCR",indat)
   if (length(begin) > 1) begin <- begin[2]
   mcdaHCR <- getLogical(indat[begin],1)
   begin <- grep("ConstH",indat)
   ConstH <- getLogical(indat[begin],1)
   pickSched <- getsingleNum("pickSched",indat)
   begin <- grep("TACadj",indat)
   TACadj <- getConst(indat[begin],11,2)
   begin <- grep("TACadj2",indat)
   TACadj2 <- getConst(indat[begin],11,2)
   begin <- grep("mcdaWts",indat)
   mcdaWts <- getConst(indat[begin],3,2)
   runlabel <-  paste0("_",assessInterval,"_",pickSched,"_",mcdaWts[1],"_",mcdaWts[2])
   begin <- grep("postmcdaWts",indat)
   postmcdaWts <- getConst(indat[begin],3,2)
   begin <- grep("withVariation",indat)
   withVariation <- getLogical(indat[begin],1)
   cpuePeriod <- getsingleNum("cpuePeriod",indat)
   maxGrad4 <- getsingleNum("maxGrad4",indat)
   maxRate1 <- getsingleNum("maxRate1",indat)
   begin <- grep("CETarg",indat)
   CETarg <- getConst(indat[begin],4,2)
   begin <- grep("deltaCE",indat)
   deltaCE <- getConst(indat[begin],4,2)
   implementE <- getsingleNum("implementE",indat)
   begin <- grep("LRPTAC",indat)
   LRPTAC <- getLogical(indat[begin],1)
   TACLower <- getsingleNum("TACLower",indat)
   TACUpper <- getsingleNum("TACUpper",indat)
   refyr <- getsingleNum("refyr",indat)
   withsigR <- getsingleNum("withsigR",indat)
   withsigB <- getsingleNum("withsigB",indat)
   withsigCE <- getsingleNum("withsigCE",indat)
   outctrl <- list(batch,reps,initDepl,assessInterval,runlabel,
                   recthreshold,hcrLabel,mcdaHCR,ConstC,mcdaWts,postmcdaWts,
                   pickSched,TACadj,TACadj2,withVariation,cpuePeriod,maxGrad4,
                   maxRate1,CETarg,deltaCE,
                   implementE,LRPTAC,TACLower,TACUpper,refyr,ConstH,
                   withsigR,withsigB,withsigCE)
   names(outctrl) <- c("batch","reps","initDepl","assessInterval","runlabel",
                       "recthreshold","hcrLabel","mcdaHCR","ConstC","mcdaWts",
                       "postmcdaWts","picksched","TACadj","TACadj2",
                       "withVariation","cpuePeriod","maxGrad4","maxRate1",
                       "CETarg","deltaCE","implementE","LRPTAC","TACLower",
                       "TACUpper","refyr","ConstH","withsigR","withsigB",
                       "withsigCE")
   return(outctrl)
} # end of readhcrfile


#' @title readregionfile reads in the constants for the region
#'
#' @description with the region filename from the control file the
#'     readregionfile will read the data from a csv file arranged with
#'     a standard layout. Once again this uses the utility functions
#'     for reading and parsing lines of text. This is illustrated by
#'     the function makeregionfile, which produces an example *.csv
#'     file that can then be customized to suit your own simulations.
#'     Each required section contains a series of constants which are
#'     read in individually,so their labels are equally important.
#'
#' @param indir directory in which to find the region file
#' @param infile character string with filename of the region file
#'
#' @return eight objects, four vectors, 4 numbers, and a list of five
#'     objects
#' \itemize{
#'   \item SMUnames the labels given to each SMU
#'   \item SMUpop the number of populations in each SMU, in sequence
#'   \item minc scaler defining the mid point of the smallest size class
#'   \item cw scaler defining the class width
#'   \item randomseed used to ensure that each simulation run starts
#'       in the same place. Could use getseed() to produce this.
#'   \item outyear a vector of three, with Nyrs, fixyear, and firstyear
#'   \item projLML the LML expected in each projection year, a vector
#'   \item globals a list containing numpop, nSMU, midpts, Nclass, and
#'       Nyrs
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' datadir <- "./../../rcode2/aMSE/data-raw/"
#' ctrlfile <- "control.csv"
#' ctrl <- readctrlfile(datadir,ctrlfile)
#' reg1 <- readregionfile(datadir,ctrl$regionfile)
#' str(reg1)
#' }
readregionfile <- function(indir,infile) {  # infile="region1.csv"; indir=datadir
   context <- "region file"
   filename <- filenametopath(indir,infile)
   indat <- readLines(filename)   # reads the whole file as character strings
   nSMU <-  getsingleNum("nSMU",indat) # number of spatial management units
   begin <- grep("SMUpop",indat)
   SMUpop <-  getConst(indat[begin],nSMU) #
   numpop <- sum(SMUpop)
   SMUnames <- getStr(indat[begin+1],nSMU)
   minc <-  getsingleNum("minc",indat) # minimum size class
   cw    <- getsingleNum("cw",indat) # class width
   Nclass <- getsingleNum("Nclass",indat) # number of classes
   midpts <- seq(minc,minc+((Nclass-1)*cw),2)
   larvdisp <- getsingleNum("larvdisp",indat)
   randomseed <- getsingleNum("randomseed",indat)
   Nyrs <- getsingleNum("Nyrs",indat)
   firstyear <- getsingleNum("firstyear",indat)
   fixyear <- getsingleNum("fixyear",indat)
   outyear <- c(Nyrs,fixyear,firstyear)
   projLML <- numeric(Nyrs)
   begin <- grep("PROJLML",indat)
   for (i in 1:Nyrs) {
      from <- begin + 1
      projLML[i] <- getConst(indat[from],1)
   }
   begin <- grep("CONDITION",indat)
   condition <- getConst(indat[begin],1)
   histLML <- NA
   if (condition > 0) {
      begin <- grep("HistoricalLML",indat)
      ncond <- getConst(indat[begin],1)
      histLML <- numeric(ncond)
      for (i in 1:ncond) {
         begin <- begin+1
         histLML[i] <- getConst(indat[begin],1)
      }
   }
   globals <- list(numpop=numpop, nSMU=nSMU, midpts=midpts,
                   Nclass=Nclass, Nyrs=Nyrs, larvdisp=larvdisp)
   totans <- list(SMUnames,SMUpop,minc,cw,larvdisp,randomseed,outyear,
                  projLML,condition,histLML,globals)
   names(totans) <- c("SMUnames","SMUpop","minc","cw","larvdisp",
                      "randomseed","outyear","projLML","condition",
                      "histLML","globals")
   return(totans)
}  # end of readregionfile

#' @title regionfiletemplate generates a template region file
#'
#' @description regionfiletemplate generates a standardized region file
#'     template containing constants tha toperate at a regional level.
#'     Generate this and then modify the contents to suit
#'     the system you are attempting to simulate.
#'
#' @param indir directory in which to place the region file
#' @param filename the name for the generated region file, a character
#'     string that defaults to region1.csv. It is best to include the
#'     complete path
#'
#' @return nothing, but it creates a template region file within indir
#'
#' @export
#'
#' @examples
#'  yourdir <- tempdir()
#'  regionfiletemplate(yourdir,filename="region1.csv")   #
#'  region1 <- readregionfile(yourdir,"region1.csv")
#'  str(region1,max.level=1)
regionfiletemplate <- function(indir,filename="region1.csv") {
   filename <- filenametopath(indir,filename)
   cat("region file containing region wide constants for a run \n",
       file=filename,append=FALSE)
   cat("Modify the contents to suit your own situation \n",
       file=filename,append=TRUE)
   cat("\n",file=filename,append=TRUE)
   cat("REGION \n",file=filename,append=TRUE)
   cat("nSMU, 2, number of spatial management units eg 2 \n",
       file=filename,append=TRUE)
   cat("SMUpop, 2, 4, number of populations per SMU in sequence \n",
       file=filename,append=TRUE)
   cat("SMUnames, block1, block2, labels for each SMU \n",file=filename,
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
   cat("YEARS \n",file=filename,append=TRUE)
   cat("Nyrs, 30, number of projection years for each simulation \n",
       file=filename,append=TRUE)
   cat("firstyear, 2014, first year of simulation if 1 hypothetical 2014 conditioned \n",
       file=filename,append=TRUE)
   cat("fixyear, 3, conditions fixed for the first three years \n",
       file=filename,append=TRUE)
   cat("\n",file=filename,append=TRUE)
   cat("PROJLML, same number of projection years \n",file=filename,
       append=TRUE) # ensure there are Nyrs lines
   cat("projLML1, 132,  Legal Minimum Length (LML, MLL, MLS) e.g. 140 \n",
       file=filename,append=TRUE)
   for (i in 2:30) {
      label <- paste0("projLML",i)
      cat(label,"132 \n",file=filename,append=TRUE)
   }
   cat("\n",file=filename,append=TRUE)
   cat("CONDITION, 0, if > 1 then how many years in the histLML \n",
       file=filename,append=TRUE)
   cat("\n",file=filename,append=TRUE)
   cat("HistoricalLML, 0, if >1 then number reflects number of histLML \n",
       file=filename,append=TRUE)
} # end of regionfiletemplate




#' @title replaceVar replaces values of a variable in the input datafile
#'
#' @description replaceVar replaces the values of a variable in the input
#'     datafile. TAKE CARE, it overwrites the original file! However, it
#'     also saves the original file by adding an '_old' to the filename.
#'     Obviously if the function is used multiple times the original
#'     'original' file will be over-written on the second use.
#'     Alternatively, one could readin the datafile using readdateFile and
#'     then directly alter the condDat, e.g., in the case of a four block
#'     zone in which the AvRec vale is changed: condDat$constants["AvRec",]
#'     <- c(12.5,12.2,12.4,12.1)
#' @param infile the name, and path,  of the data file to be altered
#' @param invar the text name of the variable to be changed
#' @param newval the new value with which to replace the current values
#' @return The function over-writes the original file but saves the original
#'     by adding an '_old' to the end of the original filename.
#' @export
#' @examples
#' \dontrun{
#' filename <- datafiletemplate(numblock=1,filename="oneblock.csv")
#' replace(filename,"AvRec",15.75)
#' condDat <- readdatafile(filename)
#' print(round(condDat$constants,4))
#' }
replaceVar <- function(infile,invar,newval) {
  indat <- readLines(infile)   # reads the whole file as character strings
  filenamelen <- nchar(infile)
  oldfilename <- paste0(substr(infile,1,(filenamelen-4)),"_old.csv")
  write(indat,file=oldfilename) #,row.names=FALSE)
  begin <- grep("BLOCKNAMES",indat)
  numblk <- getConst(indat[begin],1) # number of blocks
  begin <- grep(invar,indat)
  if (length(begin) > 0) {
    begin <- begin[1]
  } else {
    stop("Error in 'replace' unknown variable used \n")
  }
  ans <- newval
  if (numblk > 1) for (i in 2:numblk) ans <- paste(ans,newval,sep=",")
  ans <- paste(invar,ans,",",sep=",")
  indat[begin] <- ans
  write(indat,file=infile) #,row.names=FALSE)
}  # end of replaceVar

# end of utility functions used in parseFile


#' @title resetSel alters the selectivity of all popualtions
#'
#' @description resetSel alters the selectivity of all populations before
#'     proceeding with further dynamics
#'
#' @param inzone the simulated zone
#' @param inLML the new LML
#' @param glob the global variables object
#'
#' @return a revised zone
#' @export
#'
#' @examples
#' \dontrun{
#'   print("develop an example.")
#'
#' }
resetSel <- function(inzone,inLML,glob) {  # inzone <- zone; inLML <- testLML
   outzone <- inzone
   Nyrs <- glob$Nyrs
   nblock <- glob$nblock
   Nclass <- glob$Nclass
   midpts <- glob$midpts
   numpop <- glob$numpop
   useLML <- matrix(rep(inLML,Nyrs),nrow=Nyrs,ncol=nblock,byrow=TRUE)
   zSelect <- matrix(0,nrow=Nclass,ncol=Nyrs,dimnames=list(midpts,1:Nyrs))
   for (pop in 1:numpop) {  #  pop <- 1
      popparam <- inzone[[pop]]$popdef
      blk <- popparam["block"]
      verLML <- unique(useLML[,blk])
      selL50 <- popparam["SelP1"]
      selL95 <- popparam["SelP2"]
      for (LML in verLML) {
         Sel <- logistic((LML+selL50),(LML+selL95),midpts)    # uses SelP1 and SelP2
         pick <- which(useLML[,blk] == LML)
         zSelect[,pick] <- rep(Sel,length(pick))
         outzone[[pop]]$LML[pick] <- LML
      }
      outzone[[pop]]$Select <- zSelect
      outzone[[pop]]$SelWt <- zSelect * outzone[[pop]]$WtL
   }
   return(outzone)
}

#' @title scaletoOne divides a vector by its means to scale it to mean = 1.0
#'
#' @description scaletoOne divides a vector by its mean which scales it
#'     to have a mean of 1.0
#'
#' @param invect te vector in need of re-scaling
#'
#' @return a rescaled vector
#' @export
#'
#' @examples
#' \dontrun{
#'   x <- 1:9
#'   y <- scaletoOne(x)
#'   cbind(x,y)
#' }
scaletoOne <- function(invect) {
   avCE <- mean(invect,na.rm=T)
   invect <- invect/avCE
   return(invect)
} # end of scaleto1

