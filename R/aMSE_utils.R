# x <- hutils::describefunctions("C:/Users/User/Dropbox/A_code/aMSE/R/","aMSE_utils.R")


#' @title putNA can add NAs to the start and end of a vector
#'
#' @description putNA fulfils a common requirement to expand a vector with NAs
#'     to assist a plot by ensuring the x and y vectors are the same length.
#'     It can accept both character and numeric vectors. This used to be called
#'     addNA, but there is a function in base R with that name that does
#'     something completely different.
#'
#' @param x the vector to which NAs are to be added
#' @param pre how many NA to be added the front
#' @param post how many to be added to the end
#'
#' @return a vector of length(x) + pre + post
#' @export
#'
#' @examples
#' vect <- rnorm(10,mean=5,sd=1)
#' putNA(vect,3,5)
putNA <- function(x,pre,post) { # x=med14; pre=0; post=5
  n <- length(x) + pre + post
  if (pre > 0) x <- c(rep(NA,pre),x)
  if (post > 0) x <- c(x,rep(NA,post))
  return(x)
} # end of addNA

#' @title addpops adds the populations from a single replicate together
#'
#' @description addpops adds the populations from a single replicate together
#'     to form zone totals. It can only be used on matureB, exploitB, catch,
#'     acatch, and recruits. the sum of acatch is already available as TAC in
#'     the projection dynamics object 'zoneDP'
#'
#' @param invar the summable variable eg catch, recruits, matureB
#' @param nyrs the number of years of data
#' @param reps the number of replicates
#'
#' @return a 2D matrix of yrs x reps
#' @export
#'
#' @examples
#' print("wait on data")
addpops <- function(invar,nyrs,reps) {  # invar=invar; glb=glb
  result <- array(0,dim=c(nyrs,reps),dimnames=list(1:nyrs,1:reps)) # 3D to 2D
  for (iter in 1:reps) result[,iter] <- rowSums(invar[,,iter],na.rm=TRUE)
  return(result)
} # end of addpops

#' @title catchweightCE uses historical catch-by-sau to make weighted zone cpue
#'
#' @description catchweightCE is used when characterizing the historical
#'     fisheries data. It uses matching catch-by-sau to generate catch-weighted
#'     estimates of zone-wide CPUE through time.
#'
#' @param cedat historical cpue from condC
#' @param cdat historical catches from condC
#' @param nsau the number of SAU, from glb
#'
#' @return a vector of catch-weghted cpue for the zone from SAU cpue estimates
#' @export
#'
#' @examples
#' print("wait on suitable data sets")
catchweightCE <- function(cedat,cdat,nsau) {
  ceyrs <- as.numeric(rownames(cedat))
  nyr <- length(ceyrs)
  totCE <- numeric(nyr)
  names(totCE) <- ceyrs
  for (yr in 1:nyr) { #  yr = 1
    picks <- which(cedat[yr,] > 0)
    wts <- numeric(nsau)
    totC <- sum(cdat[yr,picks])
    wts[picks] <- cdat[yr,picks]/totC
    totCE[yr] <- sum(cedat[yr,] * wts,na.rm=TRUE)
  }
  return(totCE)
} # end of catchweightCE

#' @title changevar can alter the value in a single line in, eg, the control file
#'
#' @description changevar is a DANGEROUS function for lazy people. I say that
#'     because it can damage your control files if not used carefully. It speeds
#'     changing a single value within, say, the control.csv file. For example,
#'     if one wanted to conduct a retrospective analysis on what would have
#'     happened should we have introduced the HS sooner, one could use
#'     something like the following to change the value to 44 implying to start
#'     projections in 2017, in stead of 2021, which would have3 been 48.
#'     changevar(varname="CATCHES",newvalue=44,filename="controlsau.csv",
#'               rundir=rundir)
#'      This can only change a single value on a single line but can be used to
#'      alter the values of the variables listed in 'goodnames'. Currently, this
#'      only works for Windows machine. Let me know if others are required.
#'
#' @param filename name of the file to be changed. eg 'control.csv'
#' @param rundir name of the scenario directory in which the file is to be found
#' @param varname the name of the object whose value is to be changed
#' @param newvalue the new value of the object
#' @param prompt should allowable names first be listed with current values?
#' @param verbose should concluding remarks to console be made. default=TRUE
#'
#' @seealso{
#'  \link{findlinenumber}, \link{changeline}
#' }
#'
#' @return nothing but a file is changed - be careful
#' @export
#'
#' @examples
#' print("Wait on some time passing")
changevar <- function(filename,rundir,varname,newvalue,prompt=FALSE,verbose=TRUE) {
  #  varname = "CATCHES"; newvalue=57; filename="controlsau.csv";rundir=rundir
  filen <- paste0(rundir,"/",filename)
  dat <- readLines(filen)
  goodnames <- c("larvdisp","PROJECT","CATCHES","replicates",
                 "withsigR","withsigB","withsigCE","datafile")
  ngood <- length(goodnames)
  if (prompt) {
    cat("\n Allowable names and current values \n")
    for (i in 1:ngood) {
      cat(goodnames[i],"  ",getsingleNum(goodnames[i],dat),"\n")
    }
    if (Sys.info()["sysname"] == "Windows") {
      ans <- winDialog(type = "yesno", message = "Make the value change?")
      ans <- ifelse(ans == "YES", 1L, 2L)
    }
    if (ans == 2) {
      if (verbose) {
        label <- paste0("No change made to ",varname)
      } else { label <- ""}
      return(cat(label," \n"))
    }
  }
  if (varname %ni% goodnames) {
    print(goodnames)
    stop("Input varname in changevar, not in allowable names \n")
  }
  pick <- grep(varname,dat)
  txtbits <- removeEmpty(unlist(strsplit(dat[pick],",")))
  newtext <- paste0(txtbits[1],",",as.character(newvalue))
  if (length(txtbits) > 2) newtext <- paste0(newtext,",",txtbits[3])
  dat[pick] <- newtext
  writeLines(dat,con=filen)
  if (verbose) { label <- paste0("Set ",varname," to ",newvalue)
  } else { label <- "" }
  return(cat(label," \n"))
} # end of changevar

#' @title checksizecompdata is a utility to allow a check of their sizecomp data
#'
#' @description checksizecompdata is a utility which allows a user to check the
#'     coverage and level of sampling of size-composition data prior to
#'     analysis in the MSE (or sizemod) so that a selection can be made to
#'     potentially excluded inadequate samples. This is not a required function
#'     for every analysis it is designed simply as a tool for quality control of
#'     the input data.
#'
#' @param rundir the directory in which all files relating to a
#'     particular run are to be held.
#' @param controlfile default="control.csv", the filename of the control
#'     file present in rundir containing information regarding the run.
#' @param verbose Should progress comments be printed to console, default=TRUE#'
#'
#' @return nothing but it does plot a graph
#' @export
#'
#' @examples
#' print("wait on data sets")
checksizecompdata <- function(rundir,controlfile,verbose=TRUE) {
  zone1 <- readctrlfile(rundir,infile=controlfile,verbose=verbose)
  setuphtml(rundir)
  compdat <- zone1$condC$compdat$lfs
  if (is.null(compdat))
    stop(cat("No size-composition file found in copntrol file \n"))
  palfs <- zone1$condC$compdat$palfs
  saunames <- zone1$SAUnames
  nsau <- length(saunames)
  histyr <- zone1$condC$histyr
  addtable(palfs,"sizecomp_obs_by_year_by_SAU.csv",rundir,
           category="predictedcatchN",
           caption="Number of sizecomp observations by year and SAU.")
  for (plotsau in 1:nsau) {
    lfs <- preparesizecomp(compdat[,,plotsau],mincount=0)
    yrsize <- as.numeric(colnames(lfs))
    pickyr <- match(yrsize,histyr[,"year"])
    LML <- histyr[pickyr,"histLML"]
    plotsizecomp(rundir=rundir,incomp=lfs,SAU=saunames[plotsau],lml=LML,
                 catchN=NULL,start=NA,proportion=TRUE,
                 console=FALSE)
  }
  replist <- NULL
  projy <- 0
  runnotes <- c(controlfile,paste0("SAU = ",nsau))
  make_html(replist = replist,  rundir = rundir,
            controlfile=controlfile, datafile=zone1$ctrl$datafile,
            width = 500, openfile = TRUE,  runnotes = runnotes,
            verbose = verbose, packagename = "aMSE",  htmlname = "sizecomp")
  if (verbose) cat("finished  \n")
} # end of checksizecompdata

#' @title confirmdir checks to see if a directory exists and makes it if not
#'
#' @description confirmdir enables one to be sure a selected directory exists.
#'     If it has not been created then confirmdir will create it if it does not
#'     already exist. This is useful when defining output directories on the
#'     hard drive where large objects may be stored.
#'
#' @param x the directory to be checked and created if necessary
#' @param make should the directory be created if it does not already exist?
#'     default=TRUE
#' @param verbose should responses be sent to the console? default=TRUE
#' @param ask should confirmdir ask whether to create the directory or not?
#'    default = TRUE. Set to FALSE if running a script rather than working
#'    interactively.
#'
#' @return nothing but it can create a directory
#' @export
#'
#' @examples
#' x <- tempdir()
#' confirmdir(x)
confirmdir <- function(x,make=TRUE,verbose=TRUE,ask=TRUE) {
  if (dir.exists(x)) {
    if (verbose) cat(x," already exists  \n")
  } else {
    if (verbose) cat(x," did not exist  \n")
    if (make) {
      if (ask) {
        label <- paste0("Create directory: ",x," [Y, y, N, n]: ")
        goahead <- readline(prompt=label)
      } else {
        goahead <- "Y"
      }
      if (goahead %in% c("y","Y")) {
        dir.create(x, recursive = TRUE)
        if (verbose) cat(x," has been created  \n")
      }
    }
  }
} # end of confirmdir

#' @title copyto copies a vector of files from one scenario directory to another
#'
#' @description copyto copies a vector of files (see examples) from one
#'     scenario's directory to another. If a filename includes the
#'     name of the scenario, eg controlM15h75.csv, is found in the scenario
#'     M15h75, and it is to be copied to, say, M15h5, then the filename will
#'     be changed automatically. Care should be taken when using this
#'     function as it writes files and potentially new directories to your
#'     storage drives.
#'
#' @param prefixdir the directory (path) containing the different scenarios
#' @param fromdir just the name of the directory from which to source the
#'     files listed in filelist (no path)
#' @param todir just the name of the new destination directory (no path)
#' @param filelist a vector of filenames to be copied, as character.
#' @param makenew if the 'todir' does not exist should it be created using
#'     dir.create? default = TRUE
#' @param verbose should details be printed to the console, default=TRUE
#'
#' @return Nothing but it will transfer files and change their names. This
#'     information will be printed to the console if verbose = TRUE
#' @export
#'
#' @examples
#' \dontrun{
#'  prefixdir <- "c:/Users/User/DropBox/A_codeUse/aMSEUse/scenarios/"
#'  fromdir <- "M15h75"
#'  todir <- "M15h5"
#'  vectfiles <- c("controlM15h75.csv","saudataM15h75.csv",
#'                   "lf_WZ90-20.csv","run_aMSE_M15h75.R","TasHS1_Tas.R")
#'  copyto(prefixdir=prefixdir,fromdir=fromdir,todir=todir,filelist=vectfiles)
#' }
copyto <- function(prefixdir,fromdir, todir, filelist,
                   makenew = TRUE,verbose=TRUE) {
  fdir <- filenametopath(prefixdir,fromdir)
  if (!dir.exists(fdir)) stop(cat(fdir," does not exist!   \n\n"))
  nfile <- length(filelist)
  for (i in 1:nfile) {  # check if all files in list exist
    filen <- filenametopath(fdir,filelist[i])
    if (!file.exists(filen)) stop(cat(filen," does not exist \n"))
  }
  tdir <- filenametopath(prefixdir,todir)
  if (!dir.exists(tdir)) {
    if (verbose) cat(tdir," did not exist  \n")
    if (makenew) {
      dir.create(tdir, recursive = TRUE)
      if (verbose) cat(tdir," has been created  \n")
    }
  }
  for (i in 1:nfile) { # i = 1
    filen <- filenametopath(fdir,filelist[i])
    present <- grep(fromdir,filelist[i])
    fileout <- filenametopath(tdir, filelist[i])
    if (length(present) > 0) {
      newfile <- sub(fromdir,todir,filelist[i]) # alter filenames
      fileout <- filenametopath(tdir, newfile)
    }
    file.copy(filen, fileout, overwrite = TRUE, copy.date = TRUE)
  }
  if (verbose)
    for (i in 1:nfile) cat(filelist[i], " has been copied to ",todir,"\n")
  newfilelist <- dir(tdir)  # change saudata file name in controlfile
  pickC <- grep("control",newfilelist)
  if (length(pickC) == 0)
    stop(cat("No recognizable control file in ",tdir,"\n"))
  filename <- filenametopath(tdir,newfilelist[pickC])
  indat <- readLines(filename)
  pickD <- grep("saudata",newfilelist)
  if (length(pickD) == 0)
    stop(cat("No recognizable saudata file in ",tdir,"\n"))
  pick <- grep("datafile",indat)
  indat[pick] <- paste0("datafile, ",newfilelist[pickD],
                        " , name of saudata file,")
  writeLines(indat,con=filename)
  if (verbose)
    cat("\n\n")
    cat("Be sure to change the control, data, and run files where necessary \n")
  return(tdir)
} # end of copyto


#' @title findlinenumber prints out the contents of a text file with line numbers
#'
#' @description findlinenumber solves the problem of finding the line number in
#'     a given text file when one wants to change a specific line using the
#'     function changeline. After inc lines are printed the function stops and
#'     waits for any character to be input (a space will suffice) before printing
#'     the next inc lines
#'
#' @param rundir the directory path in which to find the text file.
#' @param filename the full name of the text file in quotations.
#' @param inc the number of lines on screen
#'
#' @return nothing, but it does print the contents of a text file with the
#'     respective line numbers on the console.
#' @export
#'
#' @examples
#' print("wait on a suitable example")
findlinenumber <- function(rundir,filename,inc=20) {  # rundir=rundir; filename="controlsau.csv"
  filen <- filenametopath(rundir,filename)
  dat <- readLines(filen)
  nline <- length(dat)
  cat(nline," lines long  \n")
  for (i in 1:nline) {
    cat(i,"|  ",dat[i],"\n")
    if ((i %% inc) == 0)  x <- scan(what="character",n=1,quiet=TRUE)
  }
} # end of findlinenumber

#' @title makewidedat converts long data to a wide data format
#'
#' @description makewidedat takes the output of the commlf function makelongdat
#'     and converts the columns of year x length x count into a matrix of counts
#'     with axes of lengths x years, ready for comparisons with predicted length
#'     composition of catches from the operating model.
#'
#' @param inlong the long form data.frame derived for commlf's makelongdat
#'     containing at least 'year', 'length', and 'propcounts'
#' @param mids the range of lengths to be found within the input data file.
#'     These are matched with the size distribution for each year of data and
#'     the appropriate matrix cells filled with the propcounts.
#' @param counts default = FALSE, which means that the wide format will contain
#'     the proportion of counts for each year. If TRUE, then the actual counts
#'     by size-class will be output by year
#'
#' @return a matrix of proportional counts for each size-class x year
#' @export
#'
#' @examples
#' print("wait on data sets")
makewidedat <- function(inlong,mids,counts=FALSE) { # inlong=lf; mids=mids
  columns <- colnames(inlong)
  pickcol <- match(c("year","length","propcounts"),columns)
  label <- paste0("Data input to makewidedat does not contain ",
                  "all of year, length, and propcounts")
  if (length(pickcol) != 3) stop(label)
  years <- sort(unique(inlong[,"year"]))
  nyr <- length(years)
  nc <- length(mids)
  answer <- matrix(0,nrow=nc,ncol=nyr,dimnames=list(mids,years))
  for (yr in 1:nyr) { #  yr=1
    pickyr <- which(inlong[,"year"] == years[yr])
    pickloc <- match(inlong[pickyr,"length"],mids)
    if (counts) {
       answer[pickloc,yr] <- inlong[pickyr,"counts"]
    } else {
       answer[pickloc,yr] <- inlong[pickyr,"propcounts"]
    }
  }
  return(answer)
} # end of makewidedat


#' @title poptosau converts projected population dynamics to SAU scale results
#'
#' @description poptosau the MSE dynamics are run at the population level but
#'     all management decisions are made at the SAU level, hence the results of
#'     the projections need to be translated into results at the SAU level.
#'     poptosau uses the sauindex to sum the variables that can be summed, which
#'     include matureB, exploitB, catch, and recruit, and combines their
#'     population results to form their respective SAU results. Thus, for an
#'     input array of dimensions [projyrs,numpop,reps] one receives back an
#'     array of [projyrs, nSAU,reps]
#'
#' @param invar either zoneDP$ matureB, exploitB, catch, or recruit
#' @param glb the global constants object
#'
#' @return a results array of dimension [projyrs, nSAU,reps]
#' @export
#'
#' @examples
#' print("wait on appropriate internal datasets")
poptosau <- function(invar,glb) {  # invar=zoneDP$matureB; glb=glb
  numpop <- glb$numpop
  nsau <- glb$nSAU
  nyrs <- dim(invar)[1]
  reps <- dim(invar)[3]
  sauindex <- glb$sauindex
  saunames <- glb$saunames
  result <- array(0,dim=c(nyrs,nsau,reps),
                  dimnames=list(1:nyrs,saunames,1:reps)) #aspirational catches
  for (iter in 1:reps)
    for (yr in 1:nyrs)
      result[yr,,iter] <- tapply(invar[yr,,iter],sauindex,sum,na.rm=TRUE)
  return(result)
} # end of poptosau

#' @title poptosauCE combines population cpue into sau as catch weighted sums
#'
#' @description poptosauCE combines cpue from separate populations into their
#'     respective sau using a catch-weighted strategy. The sauindex is used to
#'     identify which populations to apply the sau total catches to.
#'
#' @param catvect the vector of catches x population for a given year
#' @param cpuevect the vector of cpue x population for a given year
#' @param sauindex the sau indices of each population
#'
#' @return a list of saucpue and saucatch
#' @export
#'
#' @examples
#' print("wait on appropriate built-in data files")
#' # catvect=zoneDD$catch[1:finalyr,]; cpuevect=zoneDD$cpue[1:finalyr,]
poptosauCE <- function(catvect,cpuevect,sauindex) {
  saucatch <- tapply(catvect,sauindex,sum,na.rm=TRUE)
  wts <- catvect/saucatch[sauindex]
  saucpue <- tapply((cpuevect * wts),sauindex,sum,na.rm=TRUE)
  return(list(saucpue=saucpue,saucatch=saucatch))
} # end of poptosauCE


#' @title poptozone translates the zone_pop objects to a single zone object
#'
#' @description poptozone combines the dynamic results for each variable so
#'     that results by population become results by zone. matureB, exploitB,
#'     catch, recruit, catchN, and Nt are simple summations of the totals for
#'     each population into their respective Zone The harvest rate would be
#'     end of year or beginning of year estimates derived from dividing the
#'     catch x zone by the exploitable biomass x zone. Similarly the deplsB and
#'     depleB are the end of year matureB and exploitB divided by their
#'     respective unfished estimated by Zone obtained using getvar(zoneC,"B0").
#'
#' @param inzone one of the zone dynamics objects containing replicates, made
#'     up of populations
#' @param NAS the numbers-at-size 4D arrays from doprojection; default=NULL so
#'     it can be ignored during conditioning
#' @param glb the object containing the global constants
#' @param B0 the sum of B0 across all populations, use getvar(zoneC,"B0")
#' @param ExB0 the sum of ExB0 across all populations use getvar(zoneC,"ExB0")
#'
#' @return a list of dynamics variables by zone
#' @export
#'
#' @examples
#' print("wait on suitable internal data sets ")
poptozone <- function(inzone,NAS=NULL,glb, B0, ExB0) {
  # inzone=zoneDP; glb=glb; B0=sum(getvar(zoneC,"B0")); ExB0=sum(getvar(zoneC,"ExB0"))
  N <- glb$Nclass
  invar <- inzone$matureB
  nyrs <- dim(invar)[1]
  reps <- dim(invar)[3]
  matureB <- addpops(invar,nyrs,reps)
  deplsB <- matureB/B0
  exploitB <- addpops(inzone$exploitB,nyrs,reps)
  depleB <- exploitB/ExB0
  midyexpB <- addpops(inzone$midyexpB,nyrs,reps)
  recruit <- addpops(inzone$recruit,nyrs,reps)
  catch <- addpops(inzone$catch,nyrs,reps)
  TAC <- addpops(inzone$acatch,nyrs,reps)
  catchN <- array(data=0,dim=c(N,nyrs,reps), # define some arrays
                  dimnames=list(glb$midpts,1:nyrs,1:reps))
  Nt <- array(data=0,dim=c(N,nyrs,reps),
              dimnames=list(glb$midpts,1:nyrs,1:reps))
  cpue <- catch
  harvestR <- catch
  popcat <- inzone$catch
  popce <- inzone$cpue
  for (iter in 1:reps) {
    harvestR[,iter] <- catch[,iter]/midyexpB[,iter]
    for (yr in 1:nyrs) {
      totC <- sum(popcat[yr,,iter],na.rm=TRUE)
      cpue[yr,iter] <- sum(popce[yr,,iter] * (popcat[yr,,iter]/totC),na.rm=TRUE)
      if (is.null(NAS)) {
        catchN[,yr,iter] <- rowSums(inzone$catchN[,yr,,iter])
        Nt[,yr,iter] <- rowSums(inzone$Nt[,yr,,iter])
      } else {
        catchN[,yr,iter] <- rowSums(NAS$catchN[,yr,,iter])
        Nt[,yr,iter] <- rowSums(NAS$Nt[,yr,,iter])
      }
    }
  }
  harvestR[1,] <- 0   # to account for the missing first year of midyexpB
  outzone <- list(matureB=matureB,exploitB=exploitB,midyexpB=midyexpB,
                  catch=catch,acatch=inzone$acatch,TAC=TAC,harvestR=harvestR,
                  cpue=cpue,recruit=recruit,deplsB=deplsB,depleB=depleB,
                  catchN=catchN,Nt=Nt)
  return(outzone)
} # end of poptozone

#' @title prepareDDNt converts the population based historical Nt to SAU-based
#'
#' @description prepareDDNt processes the zoneDD, following the conditioning on
#'     the historical fishery data, so that the population numbers-at-size are
#'     converted to SAU-based numbers-at-size. It does this for both Nt and
#'     catchN. This is to allow analysis, tabulation, and plotting at the SAU
#'     scale. Input data is 3D  Nclass x years x numpop
#'
#' @param inNt the 3 dimensional array of population numbers-at-size
#' @param incatchN the 3 dimensional array of catch numbers-at-size
#' @param glb the global constants object
#'
#' @return a list of Nt and catchN x SAU
#' @export
#'
#' @examples
#' print("wait on suitable internal data-sets")
prepareDDNt <- function(inNt,incatchN,glb) { # Nt=zoneDD$Nt; catchN=zoneDD$catchN; glb=glb
  sauindex <- glb$sauindex
  nsau <- glb$nSAU
  hyrs <- glb$hyrs
  Nc <- glb$Nclass
  catchN <- array(data=0,dim=c(Nc,hyrs,nsau), # define some arrays
                  dimnames=list(glb$midpts,1:hyrs,glb$saunames))
  Nt <- array(data=0,dim=c(Nc,hyrs,nsau),
              dimnames=list(glb$midpts,1:hyrs,glb$saunames))
  for (sau in 1:nsau) { #  yr=1; sau = 1
    pick <- which(sauindex %in% sau)
    for (yr in 1:hyrs) {
      catchN[,yr,sau] <- rowSums(incatchN[,yr,pick])
      Nt[,yr,sau] <- rowSums(inNt[,yr,pick])
    }
  }
  return(list(Nt=Nt,catchN=catchN))
} # end of prepareDDNt


#' @title rnormz when sd=0 it returns rep(mean,n) but also uses n random numbers
#'
#' @description when rnorm is used with sd=0 it returns the mean but it does not
#'     use up any of the random numbers produced by the pseudo-random number
#'     generator. Hence, if one just uses rnorm with sd=0 sometimes, then the
#'     results would differ if it was repeated with sd > 0 because a slightly
#'     different pseudo-random number sequence would eventuate. rnormz is a
#'     wrapper function that returns the mean if sd = 0, but also uses up the
#'     required number of random values.
#'
#' @param n number of random numbers needed
#' @param mean the mean of the distribution, could be a vector, default=0
#' @param sd the standard deviation of the distribution, could be a vector,
#'     default = 1
#'
#' @seealso{
#'   \link{rlnormz}
#' }
#'
#' @return a vector of length n or Normal random numbers of mu = mean and stdev
#'     = sd, if sd=0 it returns the mean, but uses n random numbers.
#' @export
#'
#' @examples
#' set.seed(12345)
#' rnorm(5, mean=5, sd=1)
#' rnorm(1,mean=5,sd=1)
#' set.seed(12345)
#' rnorm(5,mean=5,sd=0)
#' rnorm(1,mean=5,sd=1)
#' set.seed(12345)
#' rnormz(5,mean=5,sd=0)
#' rnorm(1,mean=5,sd=1)
rnormz <- function(n,mean=0,sd=1) {
  if (sd > 0) {
    return(rnorm(n,mean=mean,sd=sd))
  } else {
    tmp <- rnorm(n,mean=mean,sd=1) # use up n random numbers
    return(rep(mean,n))
  }
} # end of rnormz

#' @title rlnormz when sd=0 it returns rep(mean,n) but also uses n random numbers
#'
#' @description when rlnorm is used with sd=0 it returns the meanlog but it does
#'     not use up any of the random numbers produced by the pseudo-random number
#'     generator. Hence, if one just uses rlnorm with sd=0 sometimes, then the
#'     results would differ if it was repeated with sd > 0 because a slightly
#'     different pseudo-random number sequence would eventuate. rlnormz is a
#'     wrapper function that returns the meanlog if sdlog = 0, but also uses up
#'     the required number of random values.
#'
#' @param n number of random numbers needed
#' @param meanlog the mean of the distribution, on the log scale, default=0
#' @param sdlog the standard deviation of the distribution, on the log scale,
#'     default = 1
#'
#' @seealso{
#'   \link{rnormz}
#' }
#'
#' @return a vector of length n of Log-Normal random numbers, if sd=0 it returns
#'     the meanlog, but uses n random numbers.
#' @export
#'
#' @examples
#' set.seed(12345)
#' rlnorm(5, meanlog=5, sdlog=1)
#' rlnorm(1,meanlog=5,sdlog=1)
#' set.seed(12345)
#' rlnorm(5,meanlog=5,sdlog=0)
#' rlnorm(1,meanlog=5,sdlog=1)
#' set.seed(12345)
#' rlnormz(5,meanlog=5,sdlog=0)
#' rlnorm(1,meanlog=5,sdlog=1)
rlnormz <- function(n,meanlog=0,sdlog=1) {
  if (sdlog > 0) {
    return(rlnorm(n,meanlog=meanlog,sdlog=sdlog))
  } else {
    tmp <- rlnorm(n,meanlog=meanlog,sdlog=1) # use up n random numbers
    return(rep(exp(meanlog),n))
  }
} # end of rlnormz

#' @title sautopop translates a vector of SAU properties into population properties
#'
#' @description sautopop uses the sauindex object to distribute a set of SAU
#'     properties (for example recruitment deviates) into a set of population
#'     properties.
#'
#' @param x the vector of length nsau, of properties
#' @param sauindex the index of which populations are in which SAU
#'
#' @return a vector of length numpop with each population having the property
#'     of the respective SAU
#' @export
#'
#' @examples
#' sauprop <- c(1,2,3,4,5)  # 5 SAU each with 3 populations
#' sauind <- c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5) # the sauindex
#' sautopop(sauprop,sauind)
sautopop <- function(x,sauindex) {
  nsau <- length(x)
  numpop <- length(sauindex)
  ans <- numeric(numpop)
  for (i in 1:nsau) {
    pick <- which(sauindex == i)
    ans[pick] <- x[i]
  }
  return(ans)
} # end of sautopop

#' @title save_hsargs sends a copy of the hs arguments to rundir
#'
#' @description save_hsargs saves a copy of the HS arguments (held in hsargs)
#'     to the rundir directory ready to be printed as text into the HSperf
#'     tab of the output webpage. They are stored in a file called hsargs.txt.
#'     The function expects hsargs to be a list, which it prints component
#'     by component to a txt file.
#'
#' @param rundir the data and results direcotry for the scenario
#' @param hsargs the harvest strategy arguments object.
#'
#' @return it saves a text file into the rundir sub-directory and modifies
#'     the resultTable.csv file for the web-page summary of scenario results
#' @export
#'
#' @examples
#' print("wait on suitable internal data sets")
save_hsargs <- function(rundir,hsargs) {
  filen <- filenametopath(rundir,"hsargs.txt")
  cat(names(hsargs)[1],hsargs[[1]],"\n",file=filen,append=FALSE)
  for (i in 2:length(hsargs))
    cat(names(hsargs)[i],hsargs[[i]],"\n",file=filen,append=TRUE)
  cat("\n\n",file=filen,append=TRUE)
  resfile <- filenametopath(rundir,"resultTable.csv")
  cat(c(filen,"HSperf",type="txtobj",as.character(Sys.time()),
        caption="HS argument Settings "," \n"),
      file=resfile,sep=",",append=TRUE)
} # end of save_hsargs

#' @title saveobject is used to save RData files of particular objects
#'
#' @description saveobject is used to save RData files of particular objects.
#'     Currently, after projecting 56 populations for 100 replicates the final
#'     'out' object is over 500Mb, and zoneDP, the projection object, makes up
#'     430Mb of that. Save that if you wish, but if, say, one wanted only to
#'     save the outzone object (which summarizes outcomes at the zone level),
#'     that is only 22kb. This function facilitates the saving process.
#'
#' @param obname the name of the object to be saved, as a character string
#' @param object the parent object from which the object$obname is to extracted
#' @param postfix any postfix addition you want for the name default=""
#' @param rundir the run directory for the scenario
#'
#' @return nothing but it does save a file to the rundir
#' @export
#'
#' @examples
#' print("wait on tempdir use")
#' # obname="outzone"; postfix="test"; object=out; rundir=rundir
saveobject <- function(obname,object,postfix="",rundir) {
  if (nchar(postfix) == 0) {
    outfile <- paste0(rundir,"/",obname,".RData")
  } else {
    outfile <- paste0(rundir,"/",obname,postfix,".RData")
  }
  x <- object[[obname]]
  save(x,file=outfile)
} # end of save object

#' @title summarizeprod generates a summary of the productivity properties
#'
#' @description summarizeprod generates a summary of the productivity properties
#'     by examining the productivity matrix for each SAU/population and
#'     extracting the Bmsy, annualH, MSY, Depletion, and RelCE at the maximum
#'     catch level (which approximates the MSY). It summarizes the total zone
#'     by summing all the productivity matrices and search for the largest
#'     catch again. It generates estimates of the annualH, depletion and RelCE
#'     by using a weighted average of those values from the separate SAU or
#'     populations, where the weighting is the proportion of the sum of the
#'     MSYs taken in each sau or population. This latter is only an
#'     approximation but provides at least an indication.
#'
#' @param product The productivity array from doproduction containing the
#'     range of imposed harvest rates, and the resulting outputs for each
#'     population
#' @param saunames the names of the different SAU
#'
#' @return a matrix containing the approximate productivity matrix for the zone
#' @export
#'
#' @examples
#' \dontrun{
#' data(zone)
#' product <- zone$product
#' zoneprod <- summarizeprod(product,saunames=zone$zone1$SAUnames)
#' round(zoneprod,3)
#' }
summarizeprod <- function(product,saunames) {
  # product=out$production; saunames=out$glb$saunames
  numrow <- dim(product)[1]
  numpop <- dim(product)[3]
  prodrows <- rownames(product[,,1])
  prodcols <- colnames(product[,,1])
 # nSAU <- length()  # here numpop = nSAU
  columns <- c("Bmsy","AnnH","MSY","Deplet","RelCE")
  rows <- c(saunames,"zone")
  nrows <- length(rows)
  ans <- matrix(0,nrow=nrows,ncol=length(columns),dimnames=list(rows,columns))
  for (pop in 1:numpop) { # sau=1
    sauprod <- product[,,pop]
    pick <- which(sauprod[,"Catch"] == max(sauprod[,"Catch"],na.rm=TRUE))
    ans[pop,] <- sauprod[pick,2:6]
  }
  zoneprod <- matrix(0,nrow=numrow,ncol=length(prodcols),
                     dimnames=list(prodrows,prodcols))
  for (i in 1:numpop) zoneprod <- zoneprod + product[,,i]
  pick <-  which(zoneprod[,"Catch"] == max(zoneprod[,"Catch"],na.rm=TRUE))
  ans[nrows,] <- zoneprod[pick,2:6]
  pmsy <- ans[1:numpop,"MSY"]/sum(ans[1:numpop,"MSY"],na.rm=TRUE)
  ans[nrows,"AnnH"] <- sum(ans[1:numpop,"AnnH"] * pmsy,na.rm=TRUE)
  ans[nrows,"Deplet"] <- sum(ans[1:numpop,"Deplet"] * pmsy,na.rm=TRUE)
  ans[nrows,"RelCE"] <- sum(ans[1:numpop,"RelCE"] * pmsy,na.rm=TRUE)
  return(ans)
} # end of summarizeprod


#' @title sumpop2sau gathers population data into sau data using sauindex
#'
#' @param invect a vector of population values for a given variable
#' @param sauindex the indices of each sau for each population
#'
#' @return a vector of length nsau containing the sum of population values
#'     for each sau
#' @export
#'
#' @examples
#' vect <- c(5.8,6.2,13.2,23.8,3.3,3.7,29.7,26.3,38.9,9.1)
#' sauind <- c(1,1,2,2,3,3,4,4,5,5)
#' sumpop2sau(vect,sauind)   # should be 12 37 7 56 48
sumpop2sau <- function(invect,sauindex) {
  return(tapply(invect,sauindex,sum,na.rm=TRUE))
} # end of sumpop2sau

#' @title zonetosau translates the zonexpop objects to zonexsau objects
#'
#' @description zonetosau combines the dynamic results for each variable so
#'     that results by population become results by SAU. matureB, exploitB,
#'     catch, recruit, catchN, and Nt are simple summations of the totals for
#'     each population into their respective SAU. The harvest rate would be
#'     end of year or beginning of year estimates derived from dividing the
#'     catchxsau by the exploitable biomass x sau. Similarly the deplsB and
#'     depleB are the end of year matureB and exploitB divided by their
#'     respective unfished estimated by SAU obtained using getvar(zoneC,"B0").
#'
#' @param inzone the projected zone dynamics objects made up of populations
#' @param NAS the numbers-at-size 4D arrays from doprojection; default=NULL so
#'     it can be ignored during conditioning
#' @param glb the object containing the global constants
#' @param B0 the estimate B0 for each population use getvar(zoneC,"B0")
#' @param ExB0 the estimate ExB0 for each population use getvar(zoneC,"ExB0")
#'
#' @return a list of dynamics variables by SAU
#' @export
#'
#' @examples
#' print("wait on suitable internal data sets ")
zonetosau <- function(inzone,NAS=NULL,glb, B0, ExB0) { # inzone=zoneDP; NAS=NAS; glb=glb; B0=B0; ExB0=ExB0
  nsau <- glb$nSAU
  sauindex <- glb$sauindex
  saunames <- glb$saunames
  N <- glb$Nclass
  catN <- dim(NAS$catchN)[1]
  invar <- inzone$matureB
  nyrs <- dim(invar)[1]
  reps <- glb$reps
  sauB0 <- tapply(B0,sauindex,sum)
  sauExB0 <- tapply(ExB0,sauindex,sum)
  matureB <- poptosau(invar,glb)
  exploitB <- poptosau(inzone$exploitB,glb)
  midyexpB <- poptosau(inzone$midyexpB,glb)
  recruit <- poptosau(inzone$recruit,glb)
  catlab <- glb$midpts[(N-catN+1):N]
  catchN <- array(data=0,dim=c(catN,nyrs,nsau,reps), # define some arrays
                  dimnames=list(catlab,1:nyrs,saunames,1:reps))
  Nt <- array(data=0,dim=c(N,nyrs,nsau,reps),
              dimnames=list(glb$midpts,1:nyrs,saunames,1:reps))
  # NumNe <- array(data=0,dim=c(N,nyrs,nsau,reps),
  #                dimnames=list(glb$midpts,1:nyrs,saunames,1:reps))
  catch <- inzone$acatch
  cpue <- catch
  harvestR <- catch
  deplsB <- catch
  depleB <- catch
  popcat <- inzone$catch
  popce <- inzone$cpue
  for (iter in 1:reps) {
    for (yr in 1:nyrs) { # iter=1; yr=1
      saudyn <- poptosauCE(popcat[yr,,iter],popce[yr,,iter],sauindex)
      cpue[yr,,iter] <- saudyn$saucpue
      catch[yr,,iter] <- saudyn$saucatch
      deplsB[yr,,iter] <- matureB[yr,,iter]/sauB0
      depleB[yr,,iter] <- exploitB[yr,,iter]/sauExB0

      for (sau in 1:nsau) { #  iter=1; yr=1; sau = 1
        pick <- which(sauindex %in% sau)
        if (length(pick) > 1) {
          if (is.null(NAS)) {
            catchN[,yr,sau,iter] <- rowSums(inzone$catchN[,yr,pick,iter])
            Nt[,yr,sau,iter] <- rowSums(inzone$Nt[,yr,pick,iter])
            #NumNe[,yr,sau,iter] <- rowSums(inzone$NumNe[,yr,pick,iter])
          } else {
            catchN[,yr,sau,iter] <- rowSums(NAS$catchN[,yr,pick,iter])
            Nt[,yr,sau,iter] <- rowSums(NAS$Nt[,yr,pick,iter])
           # NumNe[,yr,sau,iter] <- rowSums(NAS$NumNe[,yr,pick,iter])
          }
        } else {
          if (is.null(NAS)) {
            catchN[,yr,sau,iter] <- inzone$catchN[,yr,pick,iter]
            Nt[,yr,sau,iter] <- inzone$Nt[,yr,pick,iter]
           # NumNe[,yr,sau,iter] <- inzone$NumNe[,yr,pick,iter]
          } else {
            catchN[,yr,sau,iter] <- NAS$catchN[,yr,pick,iter]
            Nt[,yr,sau,iter] <- NAS$Nt[,yr,pick,iter]
           # NumNe[,yr,sau,iter] <- NAS$NumNe[,yr,pick,iter]
          }
        }
      } # end of dealing with numbers-at-size
    }
    harvestR[,,iter] <- catch[,,iter]/midyexpB[,,iter]
  }
  outsau <- list(matureB=matureB,exploitB=exploitB,midyexpB=midyexpB,
                 catch=catch,acatch=inzone$acatch,harvestR=harvestR,cpue=cpue,
                 recruit=recruit,deplsB=deplsB,depleB=depleB,
                 catchN=catchN,Nt=Nt) #,NumNe=NumNe
  return(outsau)
} # end of zonetosau



