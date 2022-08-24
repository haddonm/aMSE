
#' @title getavrec pulls out just the AvRec values for each sau for a scenario
#'
#' @description getavrec extracts the AvRec values for each SAU for a given
#'     scenario. It does this by reading in the saudata file and using grep
#'     to search for the correct line and returning the line as is. Care is
#'     required to only take the first record of 'AvRec' so as not to
#'     consider the variability in sAvRec.
#'
#' @param rundir the directory in which one finds the saudata file for the
#'     given scenario
#' @param datafile the exact name of the saudata file
#' @param nsau the number of SAU to be found
#'
#' @return a numeric vector of the average recruitment for each SAU
#' @export
#'
#' @examples
#' \dontrun{
#'   rundir <- "c:/Users/User/DropBox/A_codeUse/aMSEUse/scenarios/MhLML/M1h5/"
#'   getavrec(rundir,"saudataM1h5.csv",8)
#' }
getavrec <- function(rundir,datafile,nsau) {
  filen <- filenametopath(rundir,datafile)
  dat <- readLines(filen)
  pickA <- grep("AvRec",dat)[1] # ignore sAvRec
  avrec <- getConst(dat[pickA],nsau)
  return(avrec)
} # end of getavrec


#' @title getdynamics extracts the sau dynamics for a sau
#'
#' @description getdynamics assists with conditioning by generating the sauZone
#'     and extracting the fishery dynamics for the historical conditioning
#'     period ready for printing.
#'
#' @param sau which sau to extract, default = 1
#' @param zoneC the constants object for all populations
#' @param zoneDD the dynamic object for all populations
#' @param condC the conditioning fishery data includeding the catches and cpue
#' @param glb the globals object
#'
#' @return amatrix of the dynamics for the given sau
#' @export
#'
#' @examples
#' print("wait on data sets")
getdynamics <- function(zoneC,zoneDD,condC,glb,sau=1) {
  sauindex <- glb$sauindex
  popB0 <- getlistvar(zoneC,"B0")
  B0 <- tapply(popB0,sauindex,sum,na.rm=TRUE)
  popExB0 <- getlistvar(zoneC,"ExB0")
  ExB0 <- tapply(popExB0,sauindex,sum,na.rm=TRUE)
  sauZone <- getsauzone(zoneDD,glb,B0=B0,ExB0=ExB0)
  columns <- c("year","matB","expB","harvest","recruit",
               "predC","predce","depl")
  dyn <- matrix(0,nrow=glb$hyrs,ncol=length(columns),
                dimnames=list(glb$hyrnames,columns))
  dyn[,"year"] <- glb$hyrnames
  dyn[,"matB"] <- sauZone$matB[,sau]
  dyn[,"expB"] <- sauZone$expB[,sau]
  dyn[,"harvest"] <- sauZone$harvestR[,sau]
  dyn[,"recruit"] <- sauZone$recruit[,sau]
  dyn[,"predC"] <- sauZone$catch[,sau]
  dyn[,"predce"] <- sauZone$cpue[,sau]
  dyn[,"depl"] <- sauZone$deplsB[,sau]
  return(invisible(dyn))
} # end of getdynamics

#' @title getmaxCE identifies peak median CPUE during projections for each SAU
#'
#' @description getmaxCE is one of the HS performance statistics. Because it
#'     uses the targetCE from the TasmanianHS it is obviously Tasmanian
#'     specific.
#'
#' @param ceCI is the sauout$outCI$cpue object which contains the CI for the
#'     replicate projections
#' @param glb the globals object
#' @param targetCE the target CE from the hsargs
#'
#' @return a matrix of nSAU x 5: 5p, 50p, 95p year, and difference between
#'     the target and the maximum median value, all for each SAU
#' @export
#'
#' @examples
#' \dontrun{
#'   ceCI <- out$sauout$outCI$cpue
#'   glb <- out$glb
#'   targCE <- hsargs$maxtarg
#'   getmaxCE(ceCI=ceCI,glb=glb,targetCE=targCE)
#' }
getmaxCE <- function(ceCI,glb,targetCE) {
  nsau <- glb$nSAU
  label <- names(ceCI)
  yrs <- glb$pyrnames
  columns <- c("5p","50p","95p","year","diff")
  result <- matrix(0,nrow=nsau,ncol=5,dimnames=list(label,columns))
  for (i in 1:nsau) { # i = 1
    sauce <- ceCI[[label[i]]]
    pick <- which.max(sauce["50%",])
    result[i,] <- c(sauce[,pick],yrs[pick],sauce["50%",pick] - targetCE[i])
  }
  pyrs <- glb$pyrs
  tail(pyrs,n=10)
  return(result)
} # end of getmaxCE


#' @title getprojyrC returns cumulative catch from selected projection years
#'
#' @description getprojyrC is used to calculate the cumulative catch taken
#'     in each SAU and across the zone for selected years of each replicate.
#'     Thus, the output would be a replicates x nSAU matrix of total catches
#'     across the first 'period' years of the projection under the given
#'     harvest strategy.
#'
#' @param catsau the time series of catches summed across populations within
#'     each SAU. This array (reps x SAU x Allyears), is found within the 'out'
#'     object from the do_MSE function (as in out$catsau).
#' @param glb the globals object. Again found in 'out', as in out$glb
#' @param period how many years to cumulate; default = 10
#'
#' @return a reps x SAU matrix of summed catches
#' @export
#'
#' @examples
#' print("wait on suitable internal data sets")
#' # catsau=catch[[1]]; glb=glbc[[1]];period=5
getprojyrC <- function(catsau,glb,period=10) {
  hyrs <- glb$hyrs
  totyrs <- hyrs + glb$pyrs
  nsau <- glb$nSAU
  reps <- glb$reps
  columns <- c(glb$saunames,"zone")
  result <- matrix(0,nrow=reps,ncol=(nsau+1),dimnames=list(1:reps,columns))
  pmcat <- catsau[(hyrs+1):(hyrs+period),,]
  for (i in 1:nsau) { # i = 1
    saucat <- pmcat[,i,]
    result[,i] <- colSums(saucat)
  }
  result[,(nsau+1)] <- rowSums(result)
  return(result)
} # end of getprojyrC

#' @title getaav calculates annual absolute variation in catch
#'
#' @description getaav calculates the annual absolute change in catch
#'     for an input vector of catches, which could be across a series
#'     of years or even across different spatial units for a single
#'     year (an unusual use).
#'     The equation used is aav = 100 x sum(|Ct - Ct-1|)/(sum(Ct).
#'
#' @param invect a vector of catches
#'
#' @return a single scalar value the AAV of the input catches
#' @export
#'
#' @examples
#'   catch <- c(1,2,3,4,5,4,3,2,1)
#'   getaav(catch)  # should equal 0.32
getaav <- function(invect) { # invect=x
  nyr <- length(invect)
  totC <- sum(invect,na.rm=T)
  aac <- sum(abs(invect[2:nyr] - invect[1:(nyr-1)]))
  aav <- 0.0
  if (totC > 0.0) aav <- aac/totC
  return(aav)
} # end of getaav

#' @title getConst extracts 'nb' numbers from a line of text
#'
#' @description getConst parses a line of text and extracts 'nb' pieces of
#'     text as numbers
#'
#' @param inline text line to be parsed, usually obtained using readLines
#' @param nb the number of numbers to extract
#' @param index which non-empty object to begin extracting from?
#'
#' @return a vector of length 'nb'
#' @export
#'
#' @examples
#'   txtline <- "MaxDL , 32,32,32"
#'   getConst(txtline,nb=3,index=2)
getConst <- function(inline,nb,index=2) { # parses lines containing numbers
  ans <- numeric(nb)
  tmp <- unlist(strsplit(inline,","))
  if (length(tmp) < (nb+(index-1)))
    warning(paste("possible problem with data",tmp[1],
                  "missing comma?",sep=" "),"\n")
  count <- 0
  for (j in index:(nb+index-1)) {
    count <- count + 1
    ans[count] <- as.numeric(tmp[j])
  }
  return(ans)
}   # end getConst

#' @title getLFdata reads the observed LF-composition data
#'
#' @description getLFdata reads the observed LF-composition data, which needs
#'     to be stored in a csv file with the following format. Each SAU must use
#'     the same sequence of length-classes as rows, and the same sequence of
#'     years as columns, even if there are empty rows and columns. In addition,
#'     the first column must list the length-class and the second column lists
#'     the SAU. See the 'data-file-description.docx', and the example data files
#'     to see examples of the required format.
#'
#' @param rundir the directory containing the LF data files
#' @param filename the complete filename for the length-composition data file
#'
#' @seealso{
#'  \link{getzoneLF}, \link{getnas}
#' }
#'
#' @return a list containing the lfs array (lengths x years x sau) and palfs
#'     matrix of the count of observations in each year for each sau. This is
#'     then contained in the condC object.
#' @export
#'
#' @examples
#' print("wait on suitable internal data files")
#' # rundir=rundir; filename="lf_WZ90-20.csv"
getLFdata <- function(rundir,filename) {
  filen <- filenametopath(rundir,filename)
  lfdat <- read.csv(filen,header=TRUE)
  colnames(lfdat) <- tolower(colnames(lfdat))
  sau <- sort(unique(lfdat[,"sau"]))
  nsau <- length(sau)
  lcl <- sort(unique(lfdat[,"length"]))  # lcl = length classes
  ncl <- length(lcl)
  yrs <- colnames(lfdat)[3:dim(lfdat)[2]]
  yrs <- as.numeric(gsub("x","",yrs))
  nyrs <- length(yrs)
  lfs <- array(0,dim=c(ncl,nyrs,nsau),dimnames=list(lcl,yrs,sau))
  palfs <- matrix(0,nrow=nyrs,ncol=nsau,dimnames=list(yrs,sau))
  for (i in 1:nsau) { # i = 1
    pick <- which(lfdat[,"sau"] == sau[i])
    x <- lfdat[pick,3:(nyrs+2)]
    lfs[,,i] <- as.matrix(x)       # length frequencies
    palfs[,i] <- colSums(lfs[,,i]) # presence absence LF data
  }
  return(list(lfs=lfs,palfs=palfs))
} # end of getLFdata

#' @title getline extracts a vector of numbers from a txt or csv files
#'
#' @description getline can be used to extract the AvRec and MaxCEpars values
#'     for each SAU for a given scenario. It does this by reading in the
#'     saudata file and using grep to search for the varname, which gives the
#'     first linenumber with that name. One can extract from the saudata file
#'     as well as from the control file. The idea is to use values from that
#'     line when conditioning the operating model.
#'
#' @param rundir the directory in which one finds the saudata or control file
#'     for the given scenario
#' @param filen the exact name of the file being examined
#' @param varname what variable name is required? Take care with spelling
#' @param nobs how many numbers to return, usually nSAU
#'
#' @return a numeric vector for each SAU or of length nobs
#' @export
#'
#' @examples
#' \dontrun{
#'   rundir <- "c:/Users/User/DropBox/A_codeUse/aMSEUse/scenarios/MhLML/M1h5/"
#'   getline(rundir,"saudataM1h5.csv","AvRec",8)
#' }
getline <- function(rundir,filen,varname,nobs) {
  # rundir=rundir;filen="saudataM125h6.csv";varname="AvRec"; nobs=nsau
  filename <- filenametopath(rundir,filen)
  dat <- readLines(filename)
  pickA <- grep(varname,dat)[1] # ignore repeats
  if (length(pickA) > 0) {
    tmp <- dat[pickA]
    tmp <- removeEmpty(unlist(strsplit(tmp,",")))
    outval <- as.numeric(tmp[2:(nobs+1)])
  } else {
    warning(cat(varname," not found in ",filen," \n"))
    outval <- ""
  }
  return(outval)
} # end of getline

#' @title getlistvar extracts a vector or matrix from zoneC
#'
#' @description getlistvar extracts a vector or matrix from zoneC.
#'    If a vector of scalars, the names relate to populations, if a
#'    matrix the columns relate to populations. Only Me, R0, B0, effB0,
#'    ExB0, effExB0, MSY, MSYDepl, bLML, scalece, SaM, SAU, popdef,
#'    LML, Maturity, WtL, Emergent, and MatWt are currently valid
#'    choices. The indexvar = popdef would generate a listing of all
#'    the constants. If you only want a single constant from popdefs
#'    then use indexvar2.
#'
#' @param zoneC the constants components of the simulated zone
#' @param indexvar the name of the variable to be extracted; character
#' @param indexvar2 the name of the variable within popdef to extract
#'
#' @return either a vector or matrix of values depending on variable
#' @export
#'
#' @examples
#'  data(zone)
#'  zoneC <- zone$zoneC
#'  getlistvar(zoneC,"MSY")
#'  getlistvar(zoneC,"B0")
#'  getlistvar(zoneC,"popdef","AvRec")
getlistvar <- function(zoneC,indexvar,indexvar2="") {
  if (is.character(indexvar)) {
    fields <- c("Me","R0","B0","effB0","ExB0","effExB0","MSY",
                "MSYDepl","bLML","scalece","SaM","SAU","qest",
                "popdef","LML","Maturity","WtL","Emergent","MatWt")
    vects<- c("popdef","LML","Maturity","WtL","Emergent","MatWt")
    pick <- match(indexvar,fields)
    numpop <- length(zoneC)
    if (is.na(pick)) {
      warning(paste0(fields,"  "))
      stop("Invalid variable name attempted in getlistvar")
    }
    x <- sapply(zoneC,"[[",indexvar)
    if (fields[pick] %in% vects) {
      colnames(x) <- paste0("p",1:numpop)
    } else {
      names(x) <- paste0("p",1:numpop)
    }
    if ((indexvar == "popdef") & (nchar(indexvar2) > 0)) {
      x <- sapply(zoneC,"[[",indexvar)[indexvar2,]
      names(x) <- paste0("p",1:numpop)
    }
  } else {
    stop("Non-character variable in indexvar within getlistvar")
  }
  return(x)
} # End of getlistvar

#' @title getLogical extracts nb logicals from an input line of text
#'
#' @description getLogical obtains nb logicals from an input line
#'
#' @param inline text line to be parsed, usually obtained using readLines
#' @param nb the number of logicals to extract, if nb is longer than the
#'     number of logicals within inline the vector will contain NAs
#'
#' @return a vector of length nb
#' @export
#'
#' @examples
#'  txtline <- "Depleted, TRUE"
#'  getLogical(txtline,nb=1)
#'  txtline2 <- "calcthis, TRUE, FALSE"
#'  getLogical(txtline2,nb=2)
getLogical <- function(inline,nb) {  #inline <- txtline; nb=2
  tmp <- unlist(strsplit(inline,","))
  tmp <- removeEmpty(tmp)
  outtmp <- as.logical(as.character(tmp[2:(nb+1)]))
  return(outtmp)
}

#' @title getnas gets the numbers-at-size for all populations and the zone
#'
#' @description getnas extracts numbers-at-size for all populations,
#'     SAUs, and the zone. There will therefore be numpop + 1
#'     columns of numbers-at-size fr the particular year given.
#'
#' @param zoneD the dynamic portion of the zone
#' @param yr which yr from the range available should be summarized
#' @param glob the global variables object
#'
#' @seealso{
#'  \link{getLFdata}, \link{getzoneLF}, \link{prepareprojection},
#'  \link{doprojections}
#' }
#'
#' @return a matrix of numpop + 1 columns of numbers-at-size for each population
#'     and for the zone
#' @export
#'
#' @examples
#' print("wait on suitable data sets")
getnas <- function(zoneD,yr,glob) {# zoneD=zoneD;yr=1;glob=glb;
  numpop <- glob$numpop
  nSAU <- glob$nSAU
  columns <- c(paste0("p",1:numpop),"zone") # ,zone$SAUnames
  nas <- matrix(0,nrow=glob$Nclass,ncol=length(columns),
                dimnames=list(glob$midpts,columns))
  for (pop in 1:numpop) nas[,pop] <- zoneD$Nt[,yr,pop]
  nas[,"zone"] <- rowSums(nas[,1:numpop],na.rm=TRUE)
  return(nas)
} #end of getnas

#' @title gerprojyraavc applies the aavc function to the input catches
#'
#' @description gerprojyraavc takes in a 3D array of catches (or whatever) and
#'     applies the aavc function to obtain the average annual variation of
#'     catches across the second array dimension. The expected dimensions are
#'     years x sau x replicates. The function uses 'apply' to obtain the aavc
#'     along the year dimension for every replicate. To contrain the number of
#'     years analysed one needs to constrain the input array, eg catch[1:10,,]
#'
#' @param catches the three dimensional array of catches from the projection
#'     dynamics, years x sau x replicates
#' @param glb the gloabls object for the scenario
#'
#' @return a replicates x nsau matrix of aavc for whatever number of years were
#'     input to the function, the default is all years.
#' @export
#'
#' @examples
#' print("Wait on data sets")
getprojyraavc <- function(catches,glb) {
  nsau <- glb$nSAU
  reps <- glb$reps
  result <- matrix(0,nrow=reps,ncol=nsau,dimnames=list(1:reps,glb$saunames))
  for (i in 1:nsau) result[,i] <- apply(catches[,i,],2,getaav)
  return(result)
} # end of getprojyaavc

#' @title getprojyrs i sused to truncate an array to projected years
#'
#' @description getprojyrs is used when calculating the harvest strategy
#'     performance measures. To do that requires us to focus on the projection
#'     years, and possibly include the last year of observations for continuity.
#'     If examining the cpue, in Tasmania it is necessary to include years from
#'     1992, or at least the startCE year from hsargs = year 30, do not put a
#'     year eg 1992 in here, it needs to be an index.
#'
#' @param sauarr the projection data summarized to the sau scale. This might be
#'     matureB, exploitB, midyexpB, catch, acatch, harvestR, cpue, recruit,
#'     deplsB, or depleB. CatchN and Nt need different treatment
#' @param hyrs the number of harvested years used in teh conditioning, from
#'     the global object = out$glb$hyrs
#' @param pyrs the number of projected years used, from the global object =
#'     out$glb$hyrs
#' @param startyr enables data selection from any indexed year. This is
#'     especially useful when examining cpue as one would start from 1992 which
#'     in Tasmania, = index 30. The default = hyrs, which means the last year of
#'     observations will be included from the simulations. If hyrs+1 is used,
#'     then only the projection years would be used. if startyr=30 then in
#'     Tasmania the years from 1992 to the end of the projections will be used.
#'
#'
#' @return an array of the input data truncated to the projection years with,
#'     optionally, just the last year of the observed period included.
#' @export
#'
#' @examples
#' print("getprojyrs(out$sauout$zonePsau$catch, hyrs=out$glb$hyrs, pyrs=out$glb$pyrs, startyr=out$glb$hyrs)")
getprojyrs <- function(sauarr,hyrs,pyrs,startyr=hyrs) {
  yrs <- dim(sauarr)[1]
  if (yrs != hyrs + pyrs)
    stop(cat("Input error in getprojyrs, array in has incorrect dimensions \n"))
  return(sauarr[startyr:yrs,,])
} # end of getprojyrs



#' @title getsingleNum find a line of text and extracts a single number
#'
#' @description getsingleNum uses grep to find an input line. If the variable
#'     being searched for fails then NULL is returned
#'
#' @param varname the name of the variable to get from intxt
#' @param intxt text to be parsed, usually obtained using readLines
#'
#' @return a single number or, if no value is in the data file a NULL
#' @export
#'
#' @examples
#' \dontrun{
#'  txtlines <- c("replicates, 100","Some_other_text, 52")
#'  getsingleNum("replicates",txtlines)
#'  getsingleNum("eeplicates",txtlines)
#'  getsingleNum("other",txtlines)
#' }
getsingleNum <- function(varname,intxt) {
  begin <- grep(varname,intxt)
  if (length(begin) > 0) {
    return(as.numeric(getConst(intxt[begin],1)))
  } else {
    return(NULL)
  }
}

#' @title getsauzone summarizes zoneD into SAU and zone
#'
#' @description getsauzone rowsums the matrices of matureB, exploitB,
#'     catch, and recruit from the input zoneD into the separate SAUs
#'     and the total into the zone. The harvestR is simply the
#'     respective catch divided by the exploitB, and the cpue are the
#'     individual population cpue weighted relative to the catch taken
#'     from each population in either the SAU or the complete zone.
#'
#' @param zoneD the zoneD after nyrs of dynamics
#' @param glb the globals object
#' @param B0 is the sau based B0
#' @param ExB0 is the sau based ExB0
#'
#' @return a list of six matrices of nSAU columns of SAU summaries,
#'     and one column for the zone
#' @export
#'
#' @examples
#' print("wait on an example")
getsauzone <- function(zoneD,glb,B0,ExB0) { # zoneD=zoneDD; glb=glb; B0=B0;ExB0=ExB0
  iSAU <- glb$sauindex
  SAU <- unique(iSAU)
  nSAU <- length(SAU)
  hyrnames <- glb$hyrnames
  matB <- getsum(zoneD$matureB,iSAU); rownames(matB) <- hyrnames
  deplsB <- matB
  nsau <- glb$nSAU
  for (i in 1:nsau) deplsB[,i] <- deplsB[,i]/B0[i]
  expB <- getsum(zoneD$exploitB,iSAU); rownames(expB) <- hyrnames
  midyexpB <-  getsum(zoneD$midyexpB,iSAU); rownames(midyexpB) <- hyrnames
  depleB <- expB
  for (i in 1:nsau) depleB[,i] <- depleB[,i]/ExB0[i]
  catch <- getsum(zoneD$catch,iSAU); rownames(catch) <- hyrnames
  recruit <- getsum(zoneD$recruit,iSAU); rownames(recruit) <- hyrnames
  harvestR <- catch/midyexpB
  cpue <- catch # just to have a labelled matrix ready
  wtzone <- zoneD$catch/catch[,(nSAU+1)]
  wtsau <- zoneD$catch
  for (mu in 1:nSAU) { # mu=1
    pick <- which(iSAU == mu)
    wtsau[,pick] <- zoneD$catch[,pick]/catch[,mu]
    if (length(pick) > 1) {
      cpue[,mu] <- rowSums(zoneD$cpue[,pick] * wtsau[,pick])
    } else {
      cpue[,mu] <- zoneD$cpue[,pick] * wtsau[,pick]
    }
    cpue[1,mu] <- NA  # cpue[2,mu] # to allow for no catches in first year
  }
  cpue[,(nSAU+1)] <- rowSums(zoneD$cpue * wtzone)
  ans <- list(matB=matB,expB=expB,midyexpB=midyexpB,catch=catch,recruit=recruit,
              harvestR=harvestR,cpue=cpue,deplsB=deplsB,depleB=depleB)
  return(ans)
}  # end of getsauzone


#' @title getStr obtains a string from an input text line
#'
#' @description  getStr obtains a string from an input text line in
#'     which any parts are separated by ','. Then, after ignoring the
#'     first component, assumed to be a label, it returns the first
#'     nb parts.
#'
#' @param inline input text line with components separated by ','
#' @param nb number of parts to return
#'
#' @return a vector of character string(s)
#' @export
#'
#' @examples
#'   txt <- "runlabel, development_run, label for this particular run"
#'   getStr(txt,1)
getStr <- function(inline,nb) {
  tmp <- unlist(strsplit(inline,","))
  tmp <- removeEmpty(tmp)
  outconst <- as.character(tmp[2:(nb+1)])
  return(outconst)
} # end of getStr

#' @title getsum sums each of the main dynamics within zoneD
#'
#' @description getsum is only used by getsauzone to sum each of the main
#'     dynamic variables within zoneD, and hence is not exported.
#'
#' @param inmat what matrix within zoneD to sum into SAU and zone
#' @param index a vector containing an index of populations within SAU
#'
#' @return an nSAU+1 column matrix summarizing each SAU and the zone
#' @export
#'
#' @examples
#' print("wait on an example")
getsum <- function(inmat,index) { # inmat=zoneDD$matureB; index=zoneDD$SAU
  nSAU <- length(unique(index))
  nyr <- nrow(inmat)
  matO <- matrix(0,nrow=nyr,ncol=(nSAU+1),
                 dimnames=list(1:nyr,c(paste0("SAU",1:nSAU),"zone")))
  for (mu in 1:nSAU) { # mu = 1
    pick <- which(index == mu)
    dat <- inmat[,pick]
    if (length(pick) > 1) {
       matO[,mu] <- rowSums(dat)
    } else {
      matO[,mu] <- dat
    }
  }
  matO[,(nSAU+1)] <- rowSums(inmat)
  return(matO)
}

#' @title getunFished - extracts all data relating to year 1; unfished.
#'
#' @description getunFished - extracts all data relating to year 1;
#'    unfished. Of course, this will only work as long as movezoneYear
#'    hasn't been called, which would disrupt the contents of the first
#'    year. But if applied to zone rather than zone1 etc, this should be
#'    fine. The outputs include: R0, B0, B0crypt, popdef, MSY,
#'    MSYDepl, LML, ExB0, and Nt
#'
#' @param zoneC the constant part of the zone
#' @param zoneD the dynamic part of the zone
#' @param glb contains the global variables used everywhere
#'
#' @return a list of multiple components relating to the unfished stock this
#'     includes a list of each of the equilibrium populations and the zone as
#'     the last element of the list
#' @export
#' @examples
#' \dontrun{
#'   data(zone) # would normally use zone <- makeequilzone(rundir,"control.csv")
#'   unfish <- getunFished(zone$zoneC,zone$zoneD,zone$glb)
#'   str(unfish,max.level=2)
#' }
getunFished <- function(zoneC,zoneD,glb) {  # inzone=zone
  if((zoneC[[1]]$B0 - zoneD$matureB[1,1]) > 0.001) {
    warning("Unfished characters taken from a depleted zone in getunFished \n")
  }
  numpop <- glb$numpop
  nofished <- vector("list",(numpop+1))
  label <- c("R0","B0","B0crypt","popdef","MSY","MSYDepl","LML","ExB0",
             "Nt")
  for (pop in 1:numpop) {
    Nt <- zoneD$Nt[,1,pop]
    R0 <- zoneC[[pop]]$R0
    B0 <- zoneC[[pop]]$B0
    B0crypt <- B0 -
            sum(zoneC[[pop]]$Emergent * Nt * zoneC[[pop]]$MatWt)/1e06
    popdef <- zoneC[[pop]]$popdef
    MSY <- zoneC[[pop]]$MSY
    MSYDepl <- zoneC[[pop]]$MSYDepl
    LML <- zoneC[[pop]]$LML[1]
    ExB0 <- zoneC[[pop]]$ExB0
    ans <- list(R0,B0,B0crypt,popdef,MSY,MSYDepl,LML,ExB0,Nt)
    names(ans) <- label
    nofished[[pop]] <- ans
  }
  R0 <- 0
  B0 <- 0
  B0crypt <- 0
  MSY <- 0
  ExB0 <- 0
  Nt <- numeric(glb$Nclass)
  for (pop in 1:numpop) { #sums across all pops to go into place numpop+1
    R0 <- R0 + nofished[[pop]]$R0
    B0 <- B0 + nofished[[pop]]$B0  ## or use: sum(sapply(zone,"[[","B0"))
    B0crypt <- B0crypt + nofished[[pop]]$B0crypt
    MSY <- MSY + nofished[[pop]]$MSY
    ExB0 <- ExB0 + nofished[[pop]]$ExB0
    Nt <- Nt + nofished[[pop]]$Nt
  }
  zoneT <- list(R0,B0,B0crypt,MSY,ExB0,Nt)
  names(zoneT) <- c("R0","B0","B0crypt","MSY","ExB0","Nt")
  nofished[[numpop+1]] <- zoneT
  return(nofished)
}  # End of getUnfished



#' @title getvar a replacement for sapply to obtain scalar constants
#'
#' @description getvar is a replacement for sapply to obtain scalar
#'     constants from zoneC and is significantly faster. It should
#'     be used to obtain things like B0, R0, MSY, scalece, etc. Still
#'     need to use sapply to pull out vectors.
#'
#' @param zoneC the constants object for the zone
#' @param invar a character variable eg. "B0" or "R0"
#'
#' @return a numpop vector of the invar constants from zoneC
#' @export
#'
#' @seealso getlistvar
#'
#' @examples
#' data(zone)
#' zoneC <- zone$zoneC
#' getvar(zoneC,"MSY")
#' getvar(zoneC,"B0")
getvar <- function(zoneC,invar) {
  npop <- length(zoneC)
  recs <- numeric(npop)
  for (i in 1:npop) recs[i] <- zoneC[[i]][[invar]]
  return(recs)
} # end of getvar

#' @title getvect extracts invar from the popdef vector in zoneC
#'
#' @description getvect extracts a numpop vector of invar from the
#'     popdef vector in zoneC. Still need to use sapply to pull out
#'     complete vectors such as popdef or maturity etc.
#'
#' @param zoneC the constants object for the zone
#' @param invar a character variable eg. "steeph", "DLMax"
#'
#' @return a numpop vector of invar from the numpop popdefs in zoneC
#' @export
#'
#' @seealso getlistvar
#'
#' @examples
#' data(zone)
#' zoneC <- zone$zoneC
#' getvect(zoneC,"steeph")
getvect <- function(zoneC,invar) {
  npop <- length(zoneC)
  ans <- numeric(npop)
  for (i in 1:npop) ans[i] <- zoneC[[i]]$popdef[invar]
  return(ans)
} # end of getvect

#' @title getyr2maxce gives a nsau x reps matrix of years to max(cpue)
#'
#' @description getyr2maxce uses apply and 'which.closest' from codeutils to
#'     identify the year of maximum cpue in each replicate for each sau in the
#'     cpue array out of the projected dynamics.
#'
#' @param cpue this is the 3D array of the out$sauout$zonePsau$cpue produced by
#'     each scenario. The out objects can be from a saved RData file from
#'     each scenario. The idea being that each cpue array should be processed in
#'     turn and the matrices of the reuslts processed in whatever function calls
#'     this one.
#' @param glb the globals object for the scenario
#'
#' @seealso {
#'     \link{which.closest}
#' }
#'
#' @return a reps x nsau matrix of the year in which the maximum cpue occurs in
#'     the projections
#' @export
#'
#' @examples
#' \dontrun{
#'   glb <- out$glb  # alternative might be ans[[i]]$glb
#'   cpue <- out$sauout$zonePsau$cpue
#'   ce <- getprojyrs(ce,glb$hyrs,glb$pyrs,startyr=glb$indexCE)
#'   maxce <- getyr2maxce(ce,glb)
#'   str1(maxce)
#' }
getyr2maxce <- function(cpue,glb) {
  # cpue=cpue; glb=glb
  yrs <- as.numeric(dimnames(cpue)[[1]])
  nsau <- glb$nSAU
  reps <- glb$reps
  result <- matrix(0,nrow=reps,ncol=nsau,dimnames=list(1:reps,glb$saunames))
  for (i in 1:nsau)
    result[,i] <- apply(cpue[,i,],2,function(x) which.closest(max(x),x))
  return(result)
} # end of getyr2maxce


#' @title getzoneLF extracts all LF data from a zone across pops and years
#'
#' @description getzoneLF extracts all LF data from a zone across
#'     all populations for each year. Thus an Nclass x Nyrs X numpop
#'     matrix is compressed into a Nclass x Nyrs matrix
#'
#' @param zoneD the dynamic part of the zone
#' @param glb the global constants
#'
#' @return an Nclass x Nyrs matrix containing LF data across all
#'     populations by year
#' @export
#'
#' @examples
#' data(zone)
#' zonelf <- getzoneLF(zone$zoneD,zone$glb)
#' str(zonelf)
getzoneLF <- function(zoneD,glb) { # need to define years
  numpop <- glb$numpop
  hyrs <- glb$hyrs
  storeLF <- matrix(0,nrow=glb$Nclass,ncol=hyrs,
                    dimnames=list(glb$midpts,1:hyrs))
  for (yr in 1:hyrs)
    storeLF[,yr] <- rowSums(zoneD$Nt[,yr,])
  return(storeLF)
} # end of getzoneLF

#' @title getzoneprod zone scale summary of product matrix from doproduction
#'
#' @description getzoneprod takes in the product matrix from doproduction and
#'     sums the mature and exploitable biomassand the total catches. Then it
#'     recalculates, for the zone, the approximate annual harvest rate, the
#'     depletion level, and the relative catch rate. The annual harvest rate and
#'     relative catch rate are simply the arithmetic mean of the values for all
#'     populations, which the depletion is the column of mature biomass divided
#'     by the first value = B0
#'
#' @param product The productivity array from doproduction containing the
#'     range of imposed harvest rates, and the resulting outputs for each
#'     population
#'
#' @return a matrix containing the approximate productivity matrix for the zone
#' @export
#'
#' @examples
#' data(zone)
#' outprod <- getzoneprod(zone$product)
#' print(outprod[1:25,])
getzoneprod <- function(product) {
  numrow <- dim(product)[1]
  rows <- rownames(product[,,1])
  columns <- colnames(product[,,1])
  numpop <- dim(product)[3]
  zoneprod <- matrix(0,nrow=numrow,ncol=length(columns),
                     dimnames=list(rows,columns))
  for (i in 1:numpop) zoneprod <- zoneprod + product[,,i]
  zoneprod[,"AnnH"] <- apply(product[,"AnnH",],1,mean,na.rm=TRUE)
  zoneprod[,"Deplet"] <- zoneprod[,"MatB"]/zoneprod[1,"MatB"]
  zoneprod[,"RelCE"] <- apply(product[,"RelCE",],1,mean,na.rm=TRUE)
  return(zoneprod)
} # end of getzoneprod

#' @title getzoneprops extracts the depletion level for a given year
#'
#' @description getzoneprops extracts a set of properties for each
#'     population and summarizes them for the zone as well. These
#'     properties include effB0, matureB, legalmatB, propprot, MSY,
#'     effexB0, exploitB, SpBDepl, ExBDepl, legalDepl, MSYDepl, LML,
#'     and bLML. See the model documentation for the meaning of each
#'     of these.
#'
#' @param zoneC the constants for the zone being simulated
#' @param zoneD the dynamic part of the zone being simulated
#' @param glb the global constants
#' @param year which years information is wanted, default=1
#'
#' @return a matrix of population properties with a column of totals
#' @export
#'
#' @examples
#' data(zone)
#' ans <- getzoneprops(zone$zoneC,zone$zoneD,zone$glb)
#' print(ans)
getzoneprops <- function(zoneC,zoneD,glb,year=1) {
  numpop <- glb$numpop
  Nclass <- glb$Nclass
  sau <- getvar(zoneC,"SAU")
  B0 <- getvar(zoneC,"B0")
  ExB0 <- getvar(zoneC,"ExB0")
  blml <- getvar(zoneC,"bLML")
  msy <- getvar(zoneC,"MSY")
  msydepl <- getvar(zoneC,"MSYDepl")
  matB <- zoneD$matureB[year,]
  ExB <- zoneD$exploitB[year,]
  harvestR <- zoneD$harvestR[year,]
  catch <- zoneD$catch[year,]
  legalmatB <- numeric(numpop); names(legalmatB) <- 1:numpop
  lml <- sapply(zoneC,"[[","LML")[year,]
  for (pop in 1:numpop) { #   pop=1
    MatWt <- zoneC[[pop]]$MatWt
    first <- which(glb$midpts >= lml[pop])[1]
    legalmatB[pop] <- sum(zoneC[[pop]]$MatWt[first:Nclass] *
                            zoneD$Nt[first:Nclass,year,pop])/1000000.0
  }
  deplet <- matB/B0
  depletEx <- ExB/ExB0
  legaldepl <- legalmatB/B0
  propprot <- (matB - legalmatB)/matB
  pqest <- getlistvar(zoneC,"qest")
  label <- c("B0","matureB","legalmatB","propprot","MSY",
             "exB0","exploitB","SpBDepl","ExBDepl","legalDepl",
            "MSYDepl","LML","bLML","harvestR","catch","qest")
  ans <- rbind(B0,matB,legalmatB,propprot,msy,ExB0,ExB,deplet,depletEx,
               legaldepl,msydepl,lml,blml,harvestR,catch,pqest)
  rownames(ans) <- label
  tot <- numeric(length(rownames(ans))); names(tot) <- label
  tot[c(1:3,5:7)] <- rowSums(ans)[c(1:3,5:7)]
  tot["propprot"] <- (tot["matureB"] - tot["legalmatB"])/tot["matureB"]
  wgts <- ans["B0",]/tot["B0"]
  tot["SpBDepl"] <- sum(wgts * ans["SpBDepl",])
  tot["ExBDepl"] <- sum((ans["exB0",]/tot["exB0"]) * ans["ExBDepl",])
  tot["legalDepl"] <- sum(wgts * ans["legalDepl",])
  tot["MSYDepl"] <- sum(wgts * ans["MSYDepl",])
  tot["LML"] <- mean(ans["LML",])
  tot["bLML"] <- sum(wgts * ans["bLML",])
  tot["harvestR"] <- mean(ans["harvestR",])
  tot["catch"] <- sum(ans["catch",])
  tot["qest"] <- NA
  ans <- as.data.frame(round(cbind(ans,tot),5))
  colnames(ans) <- c(paste0("p",1:numpop),"zone")
  return(ans)
}  # end of getzoneprops
