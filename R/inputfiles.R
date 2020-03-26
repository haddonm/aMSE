
# codeutils::listFunctions("C:/A_Mal/Rcode/Abalone/AbMSE/R/inputfiles.R")

#' @title ctrlfileTemplate generates a template input control file
#'
#' @description ctrlfileTemplate generates a standardized input control
#'     file template. Generate this and then modify the contents to suit
#'     the system you are attempting to simulate. Defaults to 4 replicates,
#'     starting at 100% depletion,i.e. at B0, with an assessInterval = 1
#'     year, a very small probability of an unusual recruitment event,
#'     and using the mcdaHCR harvest strategy.
#'
#' @param filename the name for the generated ctrlfile, character string,
#'     defaults to tmpctl.csv
#'
#' @return a standard definition control file
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  ctrlfile <- paste0(yourdir,"/newctrl.csv")
#'  ctrlfileTemplate(filename=ctrlfile)   #
#'  control <- readcrtlFile(ctrlfile)
#'  str(control,max.level=1)
#' }
ctrlfileTemplate <- function(filename="tmpctl.csv") {
   cat("Modify the contents to suit your own situation \n",file=filename,
       append=FALSE)
   cat("batch,  FALSE","\n",file=filename, append=TRUE)
   cat("replicates,  4  \n",file=filename, append=TRUE)
   cat("initDepl,  1.0 \n",file=filename, append=TRUE)
   cat("assessInterval,    1  \n",file=filename, append=TRUE)
   cat("recthreshold, 0.0000001  \n",file=filename, append=TRUE)
   cat("hcrLabel,  mcdaHCR  \n",file=filename, append=TRUE)
   cat("mcdaHCR,  TRUE   \n",file=filename, append=TRUE)
   cat("ConstC,   FALSE  \n",file=filename, append=TRUE)
   cat("ConstH,   FALSE  \n",file=filename, append=TRUE)
   TACadj <- c(0.75,0.8,0.85,0.9,0.95,1.0,1.05,1.10,1.15,1.2,1.25)
   cat(paste0("TACadj,  ",makeLabel(TACadj,insep=",")," \n"),file=filename,
       append=TRUE)
   mcdaWeights <- c(0.4,0.5,0.1)
   cat(paste0("mcdaWts,  ",makeLabel(mcdaWeights,insep=",")," \n"),
       file=filename, append=TRUE)
   postmcdaWeights <- mcdaWeights  # wts after target
   cat(paste0("postmcdaWts, ",makeLabel(postmcdaWeights,insep=",")," \n"),
       file=filename,
       append=TRUE)
   cat(paste0("withVariation,    TRUE  \n"),file=filename, append=TRUE)
   # Gradient based HCR
   cat(paste0("cpuePeriod,  4   \n"),file=filename, append=TRUE)
   cat(paste0("maxGrad4,   0.2   \n"),file=filename, append=TRUE)
   cat(paste0("maxRate1,   0.4   \n"),file=filename, append=TRUE)
   # Target Based HCR
   targCE <- rep(80.0,4)
   cat(paste0("CETarg,  ",makeLabel(targCE,insep=",")," \n"),file=filename,
       append=TRUE)   # if using rep then same target in each block
   deltCE <- rep(45,4)
   cat(paste0("deltaCE,  ",makeLabel(deltCE,insep=",")," \n"),file=filename,
       append=TRUE)
   cat(paste0("implementE, 0.0  \n"),file=filename, append=TRUE) # 0 = no delay, 1 = 1 year delay, used in all HCR
   cat(paste0("LRPTAC,  TRUE   \n"),file=filename, append=TRUE)
   cat(paste0("TACLower, 500.0   \n"),file=filename, append=TRUE)# TAC can only decline to >= TACLower * origTAC
   cat(paste0("TACUpper, 1000.0  \n"),file=filename, append=TRUE)
   cat(paste0("refyr,  20     \n"),file=filename, append=TRUE)
   cat(paste0("withsigR,  0.5   \n"),file=filename, append=TRUE)
   cat(paste0("withsigB,  0.2   \n"),file=filename, append=TRUE)
   cat(paste0("withsigCE, 0.11  \n"),file=filename, append=TRUE)
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
#' @description datafiletemplate generates a standard input datafile to use
#'     as a template. It is possible to define the number of blocks and
#'     then, once the data file is created, go in an edit it appropriately
#'     to suit exactly your own needs.
#'
#' @param numblock literally the number of blocks in the simulated zone. The
#'     number of populations per block is pre-set but can easily be altered by
#'     editing the generated data file; a scalar, defaults to 3
#' @param Nyrs the number of years of projection from year 2 to Nyrs; year 1
#'     is the initial state, either depleted to a predetermined state or in an
#'     unfished state
#' @param filename the name for the generated datafile, character string,
#'     defaults to tmpdat.csv
#' @param FisheryData a boolean determining whether or not fishery data,
#'     meaning CPUE, historical catches, effort, and historical LML are to be
#'     included in the data file. These would be used in any stock reduction
#'     that becomes part of the conditioning; defaults to FALSE
#'
#' @return a standard definition data file ready to be read by readdatafile
#' @export
#'
#' @examples
#' \dontrun{
#'  filename <- dataTemplate()
#'  condDat <- readdatafile(filename)
#'  str(condDat,max.level=1)
#' }
datafiletemplate <- function(numblock=3,Nyrs=NA,filename="tmpdat.csv",
                             FisheryData=FALSE) {
   genconst <- function(invect) {
      nlab <- length(invect)
      invect <- as.character(invect)
      ans <- invect[1]
      if (nlab > 1) for (i in 2:nlab) ans <- paste(ans,invect[i],sep=",")
      return(ans)
   }
   label <- paste0("Label,",genconst(1:numblock))
   cat(label,"\n",file=filename, append=FALSE)
   cat("#SIZECLASS  \n",file=filename, append=TRUE)
   cat("minc,  2  \n",file=filename, append=TRUE)
   cat("cw,    2  \n",file=filename, append=TRUE)
   cat("Nclass, 105  \n",file=filename, append=TRUE)
   cat("\n",file=filename, append=TRUE)
   cat("BLOCKNAMES,",numblock,"\n",file=filename, append=TRUE)
   addon <- paste0("Block",1:numblock,insep=",")
   cat("Blockname,",addon,"\n",file=filename, append=TRUE)
   cat("numpop,",genconst(rep(2,numblock)),"\n",file=filename, append=TRUE)
   cat("PDFS,  32   \n",file=filename, append=TRUE)
   cat("MaxDL ,",genconst(rep(32,numblock)),"\n",file=filename, append=TRUE)
   cat("sMaxDL ,",genconst(rep(2.5,numblock)),"\n",file=filename, append=TRUE)
   cat("L50 ,",genconst(rep(123.0,numblock)),"\n",file=filename, append=TRUE)
   cat("sL50 ,",genconst(rep(5,numblock)),"\n",file=filename, append=TRUE)
   cat("L50inc ,",genconst(rep(45,numblock)),"\n",file=filename, append=TRUE)
   cat("sL50inc ,",genconst(rep(1.25,numblock)),"\n",file=filename, append=TRUE)
   cat("SigMax  ,",genconst(rep(4.6,numblock)),"\n",file=filename, append=TRUE)
   cat("sSigMax ,",genconst(rep(0.1,numblock)),"\n",file=filename, append=TRUE)
   cat("LML ,",genconst(rep(140,numblock)),"\n",file=filename, append=TRUE)
   cat("Wtb ,",genconst(rep(3.161963,numblock)),"\n",file=filename, append=TRUE)
   cat("sWtb ,",genconst(rep(0.148,numblock)),"\n",file=filename, append=TRUE)
   cat("Wtbtoa ,",genconst(rep(962.8098,numblock)),"\n",file=filename, append=TRUE)
   cat("sWtbtoa ,",genconst(rep(-14.3526,numblock)),"\n",file=filename, append=TRUE)
   cat("Me ,",genconst(rep(0.2,numblock)),"\n",file=filename, append=TRUE)
   cat("sMe ,",genconst(rep(0.003,numblock)),"\n",file=filename, append=TRUE)
   cat("Mc ,",genconst(rep(0.2,numblock)),"\n",file=filename, append=TRUE)
   cat("sMc ,",genconst(rep(0.003,numblock)),"\n",file=filename, append=TRUE)
   cat("AvRec ,",genconst(rep(13.0,numblock)),"\n",file=filename, append=TRUE)
   cat("sAvRec ,",genconst(rep(1,numblock)),"\n",file=filename, append=TRUE)
   cat("defsteep ,",genconst(rep(0.6,numblock)),"\n",file=filename, append=TRUE)
   cat("sdefsteep ,",genconst(rep(0.02,numblock)),"\n",file=filename, append=TRUE)
   cat("L50C ,",genconst(rep(126.4222,numblock)),"\n",file=filename, append=TRUE)
   cat("sL50C ,",genconst(rep(0.5,numblock)),"\n",file=filename, append=TRUE)
   cat("L95C ,",genconst(rep(136.3749,numblock)),"\n",file=filename, append=TRUE)
   cat("sL95C ,",genconst(rep(1.0,numblock)),"\n",file=filename, append=TRUE)
   cat("MaxCEpar,",genconst(rep(0.37,numblock)),"\n",file=filename, append=TRUE)
   cat("MaxCEpar,",genconst(rep(0.02,numblock)),"\n",file=filename, append=TRUE)
   cat("selL50p ,",genconst(rep(0.0,numblock)),"\n",file=filename, append=TRUE)
   cat("selL95p ,",genconst(rep(1.5,numblock)),"\n",file=filename, append=TRUE)
   cat("SaMa,",genconst(rep(-16,numblock)),"\n",file=filename, append=TRUE)
   cat("L50Mat,",genconst(rep(123.384,numblock)),"\n",file=filename, append=TRUE)
   cat("sL50Mat,",genconst(rep(4,numblock)),"\n",file=filename, append=TRUE)
   cat("\n",file=filename, append=TRUE)
   cat("CONDITIONING, 0  \n",file=filename, append=TRUE)
   cat("\n",file=filename, append=TRUE)
   cat("RANDOMSEED,  10  \n",file=filename, append=TRUE)
   cat("\n",file=filename, append=TRUE)
   cat("YEARS,   3  \n",file=filename, append=TRUE)
   if (is.na(Nyrs)) Nyrs <- 20
   firstYear <- 2014
   cat("Nyrs,",Nyrs,"\n",file=filename, append=TRUE)
   cat("firstYear,",firstYear,"\n",file=filename, append=TRUE)
   cat("fixYear, 3  \n",file=filename, append=TRUE)
   cat("\n",file=filename, append=TRUE)
   cat("PRODUCTIVITY,undefined,  \n",file=filename, append=TRUE)
   cat("\n",file=filename, append=TRUE)
   cat("PROJLML,",Nyrs,"\n",file=filename, append=TRUE)
   addon <- paste0("Block",1:numblock,insep=",")
   cat("Parameter,",addon, "\n",file=filename, append=TRUE)
   addon <- makeLabel(rep(132,numblock),insep=",")
   for (yr in firstYear:(firstYear+Nyrs-1)) cat(yr,",",addon,"\n",
                                                file=filename, append=TRUE)
   cat("data placed into ",filename,"\n\n")
   return(filename)
}  # end of dataTemplate

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
#' @description readctrlFile complements the readdatafile. The MSE requires
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
#'  direct="./../../rcode2/aMSE/data-raw"
#'  infile="control.csv"
#'  ctrl <- readctrlfile(indir=direct,infile=infile)
#'  str(ctrl)
readctrlfile <- function(indir,infile="control.csv") {
   filename <- filenametopath(indir,infile)
   indat <- readLines(filename)   # reads the whole file as character strings
   begin <- grep("START",indat) + 1
   runlabel <- getStr(indat[begin],1)
   regionfile <- getStr(indat[begin+1],1)
   datafile <- getStr(indat[begin+2],1)
   hcrfile <- getStr(indat[begin+3],1)
   outdir <-getStr(indat[begin+4],1)
   reps <- getsingleNum("replicates",indat)
   initdepl <- getsingleNum("initdepl",indat)
   assessinterval <- getsingleNum("assessinterval",indat)
   withsigR <- getsingleNum("withsigR",indat)
   withsigB <- getsingleNum("withsigB",indat)
   withsigCE <- getsingleNum("withsigCE",indat)
   outctrl <- list(runlabel,regionfile,datafile,hcrfile,
                   outdir,reps,initdepl,assessinterval,
                   withsigR,withsigB,
                   withsigCE)
   names(outctrl) <- c("runlabel","regionfile","datafile","hcrfile",
                       "outdir","reps","initdepl","assessinterval",
                       "withsigR","withsigB","withsigCE")
   return(outctrl)
}

#' @title readdatafile reads in a matrix of data defining each population
#'
#' @description readdatafile expects a matrix of probability density
#'     function definitions that are used to define the populaitons
#'     used in the simualtion. These constitute the definition of
#'     popdefs.
#'
#' @param indir directory in which to find the date file
#' @param infile character string with filename of the data file
#' @param glb the globals variable from the region file, so obviously
#'     the readregionfile needs to run before readdatafile
#'
#' @return a matrix of values defining the PDFs used to define the
#'     properties of each population. The contents of popdefs
#' @export
#'
#' @examples
#' data(region1)
#' glb <- region1$globals
#' glb
#' data(constants)
#' constants
#' ctrlfile <- "control.csv"
#' ctrl <- readctrlfile(datadir,ctrlfile)
#' reg1 <- readregionfile(datadir,ctrl$regionfile)
#' popdefs <- readdatafile(datadir,ctrl$datafile,reg1$globals)
#' print(popdefs)
readdatafile <- function(indir,infile,glb) {  # indir=datadir;infile=ctrl$datafile;glb=ctrl$globals
   numpop <- glb$numpop
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
   label <- c("numpop","nSMU","midpts","Nclass","Nyrs","larvdisp")
   globals <- vector("list",length(label))
   names(globals) <- label
   globals$numpop <- numpop
   globals$nSMU <- nSMU
   globals$midpts <- midpts
   globals$Nclass <- Nclass   # still have Nyrs to fill in
   globals$Nyrs <- Nyrs
   globals$larvdisp <- larvdisp
   totans <- list(SMUnames,SMUpop,minc,cw,larvdisp,randomseed,outyear,
                  projLML,condition,histLML,globals)
   names(totans) <- c("SMUnames","SMUpop","minc","cw","larvdisp",
                      "randomseed","outyear","projLML","condition",
                      "histLML","globals")
   return(totans)
}  # end of readregionfile




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


#' @title readdatafile reads in the constants and matrices of data required
#'
#' @description given a filename readdatafile will expect a given format.
#'     Each required section starts with the section names in capitals. The
#'     required sections are: SIZECLASS, BLOCKNAMES, PDFs, CONDITIONING,
#'     PROJLML, HISTLML, YEARS, RANDOMSEED, PRODUCTIVITY, CATCHES,
#'     PROJCATCH, CPUE, STANDARDIZED, EFFORT. An example data file is
#'     included in the package to illustrate the required formatting for the
#'     data
#'
#' @param infile a character string containing the filename of the data file
#' @return Twenty different objects, some matrices, others single numbers
#' \itemize{
#'   \item constants the array of variables as means and stdev used to
#'       define the variation within each of the populations used to define
#'       a zone
#'   \item projLML the LML expected in each block in each projection year;
#'       a matrix
#'   \item histLML LML imposed in each block in each historical year used
#'       in the conditioning
#'   \item randomseed literally the randomseed used to ensure that each
#'       simulation run starts in the same place
#'   \item catches catches in each block in each year during the
#'       conditioning
#'       period
#'   \item geomCE a matrix of the geometric mean CPUE for each block in each
#'       year of the conditioning period
#'   \item standCE a matrix of the standardized mean CPUE for each block in
#'       each year of the conditioning period
#'   \item effort a matrix of the effort imposed for each block in each year
#'       of the conditioning period
#'   \item nblock a scaler defining the number of blocks simulated
#'   \item blockNames a vector of character with the names of each block
#'   \item blkpop defines the number of populations within each block
#'   \item numpop scaler defining the total number of populations in zone
#'   \item minc scaler defining the mid point of the smallest size class
#'   \item cw scaler defining the class width
#'   \item Nclass scaler defiing the number of size classes
#'   \item midpts vector of midpoints calculated from minc, cw, and Nclass
#'   \item prodfile filename for zone production data; used in zoneStart
#'   \item projCatch a matrix of catches expected in each block in each
#'      projection year
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' #need to have an example data file before I can run an example
#' filename <- datafiletemplate(numblock=3,filename="block3.csv")
#' condDat <- readdatafile(filename)
#' str(condDat,max.level=1)
#' }
readdatafileold <- function(infile) {  # infile <- infile
   ## pickVar = "CATCHES"; dat = indat; nb <- numblk
   # indat <- filename
   indat <- readLines(infile)   # reads the whole file as character strings
   begin <- grep("SIZECLASS",indat)
   minc <-  getConst(indat[begin+1],1) # minimum size class
   cw    <- getConst(indat[begin+2],1) # class width
   Nclass <- getConst(indat[begin+3],1) # number of classes
   midpts <- seq(minc,minc+((Nclass-1)*cw),2)

   begin <- grep("BLOCKNAMES",indat)
   numblk <- getConst(indat[begin],1) # number of blocks
   blockNames <- getStr(indat[begin+1],numblk)
   popbyblk <- getConst(indat[begin+2],numblk)
   popnum <- sum(popbyblk)
   label <- c("numpop","nblock","midpts","Nclass","Nyrs")
   globals <- vector("list",length(label))
   names(globals) <- label
   globals$numpop <- popnum
   globals$nblock <- numblk
   globals$midpts <- midpts
   globals$Nclass <- Nclass   # still have Nyrs to fill in
   rows <- c("DLMax","sMaxDL","L50","sL50","L50inc","sL50inc","SigMax","sSigMax",
             "LML","Wtb","sWtb","Wtbtoa","sWtbtoa","Me","sMe","AvRec",
             "sAvRec","defsteep","sdefsteep","L50C","sL50C","deltaC","sdeltaC",
             "MaxCEpars","sMaxCEpars","selL50p","selL95p","SaMa","L50Mat","sL50Mat")
   ans <- matrix(0,nrow=length(rows),ncol=numblk,dimnames=list(rows,blockNames))
   begin <- grep("PDFS",indat)
   gr <- length(rows) # number of parameters
   pickdat <- indat[(begin+1):(begin+gr+1)]
   for (i in 1:gr) {
      ans[i,] <- getConst(pickdat[i],numblk)
   } # completed filling ans matrix
   begin <- grep("CONDITIONING",indat)
   tmp <- getConst(indat[begin],1) # number of maturity parameters
   if (tmp > 0) condition = TRUE else condition = FALSE

   fillMat <- function(pickVar,dat,nb) {
      begin <- grep(pickVar,dat)
      if (length(begin) > 0) {
         nVar <- getConst(dat[begin],1) # number of rows in the matrix
         if (nVar > 0) {
            outmat <- matrix(0,nrow=nVar,ncol=nb,dimnames=list(1:nVar,1:nb))
            pickdat <- dat[(begin+2):(begin+nVar+2)]
            for (i in 1:nVar) {
               outmat[i,] <- getConst(pickdat[i],nb)
               colnames(outmat) <- blockNames
            }
         } else {
            print(paste0("No data for ",pickVar))
            outmat <- NA
         }
      }  else {
         print(paste0("No data for ",pickVar))
         outmat <- NA
      }
      return(outmat)
   } # end of fillmat

   lmlProj <- fillMat("PROJLML",indat,numblk)
   lmlhist <- fillMat("HISTLML",indat,numblk)
   if (condition) Nyrs <- length(lmlhist[,1]) else Nyrs <- length(lmlProj[,1])
   globals$Nyrs <- Nyrs
   begin <- grep("YEARS",indat)
   #Nyrs <- getConst(indat[(begin+1)],1)     # now derived from the LML matrices
   firstYear <- getConst(indat[(begin+2)],1)
   fixYear <- getConst(indat[(begin+3)],1)
   outYear <- c(Nyrs,fixYear,firstYear)

   begin <- grep("RANDOMSEED",indat)
   randomseed <- getConst(indat[begin],1) # number of maturity parameters

   begin <- grep("PRODUCTIVITY",indat)
   prodfile <- getStr(indat[begin],1) # nname of prpoductivity file
   catches <- fillMat("CATCHES",indat,numblk)
   projCatch <- fillMat("PROJCATCH",indat,numblk)
   gcpue <- fillMat("CPUE",indat,numblk)  # geometric mean CPUE
   scpue <- fillMat("STANDARDIZED",indat,numblk)  # Standardized CPUE
   effort <- fillMat("EFFORT",indat,numblk)
   totans <- list(ans,lmlProj,randomseed,catches,gcpue,scpue,effort,numblk,
                  blockNames,popbyblk,popnum,minc,cw,Nclass,Nyrs,midpts,outYear,
                  prodfile,lmlhist,condition,projCatch,globals)
   names(totans) <- c("constants","projLML","randomseed","catches","geomCE",
                      "standCE","effort","nblock","blockNames","blkpop",
                      "numpop","minc","cw","Nclass","Nyrs","midpts","outYear",
                      "prodfile","histLML","Condition","projCatch","globals")
   return(totans)
}  # end of readdatafile

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

