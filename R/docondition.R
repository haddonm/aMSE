
#' @title adjustavrec re-fits the AvRec saudatafile values in the MSE
#'
#' @description adjustavrec is used during the conditioning of the MSE to
#'     improve the fit to the time-series of CPUE by adjusting the SAU initial
#'     value for AvRec. Once the different SAU have been characterized and the
#'     data file values settled (even when sizemod has been used), when
#'     do_condition is run the fit of the predicted cpue and size-composition
#'     data can often be improved by re-fitting the AvRec constants in each SAU
#'     to adjust their values to better suit the population characteristics once
#'     random variation is added to each population's properties. adjustavrec
#'     does this automatically. Be wary of this function because it directly
#'     alters the values of AvRec inside the saudata.csv file. If the
#'     differences remain too great it may be necessary to adhust the recdevs
#'     as well (see optimizerecdevs).
#'
#' @param rundir the full path to the directory in which all files relating to a
#'     particular run are to be held.
#' @param glb the globals object
#' @param ctrl the control object - contains datafile and controlfile names
#' @param calcpopC a function that takes the output from hcrfun and generates
#'     the actual catch per population expected in the current year.
#' @param verbose Should progress comments be printed to console, default=TRUE
#' @param lowmult multiplier to determine low range of search, default=0.6
#' @param highmult multiplier to determine the high range of search, default=1.4
#' @param iterlim the maximum number of iterations in the search for a solution,
#'     default=30
#'
#' @seealso{
#'     \link{optimizerecdevs}
#' }
#'
#' @return nothing though it will write results to the console if verbose=TRUE
#' @export
#'
#' @examples
#' \dontrun{
#'   # needs at least libraries, aMSE, hutils, makehtml, TasHS
#'   dropdir <- getDBdir()
#'   prefixdir <- paste0(dropdir,"A_codeUse/aMSEUse/scenarios/")
#'   startime <- Sys.time()
#'   postfixdir <- "S21"
#'   verbose <- TRUE
#'   rundir <- filenametopath(prefixdir,postfixdir)
#'   controlfile <- paste0("control",postfixdir,".csv")
#'   confirmdir(rundir)
#'   data(zone1)
#'   data(saudat)
#'   rewritecontrolfile(rundir,zone1,controlfile=controlfile)
#'   rewritedatafile(rundir,zone1,saudat)
#'   rewritecompdata(rundir,zone1)
#'   hsargs <- list(mult=0.1,wid = 4,targqnt = 0.55,
#'                  maxtarg = c(150,150,150,150,150,150,150,150),
#'                  pmwts = c(0.65,0.25,0.1),
#'                  hcr = c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2),
#'                  startCE = 1992)
#'   out <- do_condition(rundir,controlfile,calcpopC=calcexpectpopC,
#'                       verbose = TRUE,doproduct = TRUE,dohistoric=TRUE,
#'                       mincount=120)
#'   makeoutput(out,rundir,postfixdir,controlfile,hsfile="TasHS Package",
#'              doproject=FALSE,openfile=TRUE,verbose=FALSE)
#'   adjustavrec(rundir,out$glb,out$ctrl,calcpopC=calcexpectpopC,verbose=TRUE)
#'  # Now repeat the do_condition and makeoutput steps to see the final result
#' }
adjustavrec <- function(rundir,glb,ctrl,calcpopC,verbose=TRUE,
                        lowmult=0.6,highmult=1.4,iterlim=30) {
  startime <- Sys.time()
  datafile <- ctrl$datafile
  controlfile <- ctrl$controlfile
  nsau <- glb$nSAU
  saunames <- glb$saunames
  final <- numeric(nsau)
  original <- getavrec(rundir,datafile,nsau=nsau)
  for (sau in 1:nsau) { #  sau=3
    initial <- getavrec(rundir,datafile,nsau=nsau)
    param <- initial[sau]
    low <- param * lowmult
    high <- param * highmult
    extra <- initial[-sau]
    if (verbose) cat("Running for sau ",saunames[sau],"\n")
    origssq <- sauavrecssq(param,rundir,controlfile,
                           datafile=datafile,linenum="AvRec",
                           calcpopC=calcpopC,extra=extra,picksau=sau,
                           nsau=nsau)
    if (verbose) cat("oldssq ",origssq,"   orig param ",param,"\n")
    ans <- optim(param,sauavrecssq,method="Brent",lower=low,upper=high,
                 rundir=rundir,
                 controlfile=controlfile,datafile=datafile,linenum="AvRec",
                 calcpopC=calcpopC,extra=extra,picksau=sau,nsau=nsau,
                 control=list(maxit=iterlim))
    final[sau] <- round(ans$par,1)
    if (verbose) {
     # cat("ssq = ",ans$value,"  new R0 value = ",ans$par,"\n")
      if (((ans$par - low) < 1) | ((high - ans$par) < 0))
        warning(cat("Boundary reached for param ",sau,low,high,ans$par,"\n"))
     # cat(sau,"   ",low,"    ",param,"   ",trunc(ans$par),"   ",high,"\n")
      cat("old ssq",round(origssq,1),"     new ",round(ans$value,1),"\n")
      cat("old par",round(param,4),"     new ",round(ans$par,4),"\n\n")
    }
  }
  endtime <- Sys.time()
  if (verbose) {
    cat("Initial AvRec ",original,"\n");
    cat("Final   AvRec ",final,"\n")
    print(endtime - startime)
  }
} # end of adjustavrec

#' @title changecolumn alters a selected column of recdevs in the control file
#'
#' @description changecolumn is designed to be used when conditioning the OM
#'     when optimizing ad hoc recruitment deviates by calculating the SSQ
#'     between the predicted and observed CPUE for a selected range of years.
#'     If the length of new values is not the same as the selected linerange
#'     the function will stop with a warning. When a line is changed then all
#'     -1 values are changed to +1. Used in optimizerecdevs
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
#'     is a fine way to mess up a data file so use with care. Used in
#'     optimizerecdevs and sauavrecssq
#'
#' @param indir the directory path in which to find the text file. Usually,
#'     rundir
#' @param filename the full name of the text file in quotations.
#' @param linenumber either the line number within the text file to be changed,
#'     or, the character name of the variable to be changed, e.g. 'AvRec'
#' @param newline the character string with which to replace the line
#' @param verbose should confirmation be output to the console. default=FALSE
#'
#' @seealso{
#'  \link{findlinenumber}, \link{changevar}, \link{optimizerecdevs},
#'  \link{sauavrecssq}
#' }
#'
#' @return nothing but it does alter a line in a text file. Optionally it may
#'     confirm the action to the console
#' @export
#'
#' @examples
#' print("wait on an example")
changeline <- function(indir, filename, linenumber, newline,verbose=FALSE) {
  # indir=rundir; filename="saudataM15h5.csv"; linenumber="AvRec"; newline=chgevals
  filen <- filenametopath(indir,filename)
  dat <- readLines(filen)
  if (inherits(linenumber,"numeric")) {
    origtext <- dat[linenumber]
    dat[linenumber] <- newline
    writeLines(dat,filen)
  }
  if (inherits(linenumber,"character")) { # ie the name of the variable
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
#'     particular run are to be held.
#' @param controlfile the filename of the control file present in rundir
#'     containing information regarding the particular run.
#' @param calcpopC a function that takes the output from hcrfun (either the
#'     aspirational catch x SAU or TAC x zone) and generates the actual catch
#'     per population in each SAU expected in the current year.
#' @param verbose Should progress comments be printed to console, default=FALSE
#' @param doproduct should production estimates be made. default=FALSE
#' @param dohistoric should the historical catches be applied. Default=TRUE
#' @param matureL is a vector of 2, default = c(70,200), that define the x-axis
#'     bounds on the biology maturity-at-Length plots
#' @param wtatL is a vector of 2, default = c(80,200), that define the x-axis
#'     bounds on the biology weight-at-Length plots
#' @param mincount determines the minimum sample size for a size-composition
#'     sample to be included in plots and analyses. Default = 100
#' @param uplimH defines the upper limit of harvest rate used when estimating
#'     the productivity (also important when initial depletion is not 1.0). The
#'     default = 0.4
#' @param incH defines the interval between H steps when estimating productivity
#'     default = 0.005
#' @param deleteyrs default = 0, meaning delete no years from the sizecomp data.
#'     if there are years that are to be removed then this should be a matrix
#'     of years to delete vs sau, ie years as rows and sau as columns. All
#'     sau need to be included. to compelte hte matrix 0 values can be used.
#' @param prodpops a vector of population numbers whose productivity
#'     characteristics will be plotted into a tab called popprod, default=NULL
#'     which implies no plots will be produced.
#'
#' @seealso{
#'  \link{makeequilzone}, \link{dohistoricC}, \link{prepareprojection},
#'  \link{doprojections}, \link{prodforpop}
#' }
#'
#' @return a large list containing tottime, projtime, starttime, glb, ctrl,
#'     zoneDD, zoneDP, projC, condC, sauout, and outzone
#' @export
#'
#' @examples
#' print("wait on suitable data sets in data")
#' # rundir=rundir; controlfile=controlfile;calcpopC=calcexpectpopC
#' # verbose=TRUE; doproduct=FALSE; dohistoric=TRUE; mincount=120
#' # matureL=c(70,200);wtatL=c(80,200);mincount=120; uplimH=0.35;incH=0.005
do_condition <- function(rundir,controlfile,calcpopC,
                         verbose=FALSE,doproduct=FALSE,dohistoric=TRUE,
                         matureL=c(70,200),wtatL=c(80,200),mincount=100,
                         uplimH=0.4,incH=0.005,deleteyrs=0,prodpops=NULL) {
  starttime <- Sys.time()
  setuphtml(rundir)
  zone <- makeequilzone(rundir,controlfile,doproduct=doproduct,uplimH=uplimH,
                        incH=incH,verbose=verbose)
  equiltime <- (Sys.time()); if (verbose) print(equiltime - starttime)
  # declare main objects
  glb <- zone$glb
  ctrl <- zone$ctrl
  zone1 <- zone$zone1
  projC <- zone$zone1$projC
  condC <- zone$zone1$condC
  zoneC <- zone$zoneC
  zoneD <- zone$zoneD
  biology_plots(rundir, glb, zoneC, matL=matureL,Lwt=wtatL)
  filen <- "population_defined_properties.csv"
  addtable(condC$poprec,filen=filen,rundir=rundir,category="popprops",
           caption=paste0("Specific population properties defined rather than ",
                          "randomly allocated away from a mean."))
  production <- NULL
  if (doproduct) {
    production <- zone$product
    nprod <- length(prodpops)
    if (nprod > 0) {
      for (pop in 1:nprod) {
        pickpop <- prodpops[pop]
        prodforpop(rundir=rundir,inprod=production[,,pickpop],pop=pickpop,
                   console=FALSE)
      }
    }
  }
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
    dyn <- getdynamics(zoneC,zoneDD,condC,glb,sau=1)
    addtable(round(dyn,4),"dynamics_sau1.csv",rundir,category="zoneDD",caption=
               "dynamics of sau following conditioning on historical catches.")
  } else {
    zoneDD <- zoneD
  }
  #pop properties
  popprops <- as.data.frame(sapply(zoneC,"[[","popdef"))
  columns <- NULL
  for (sau in 1:glb$nSAU)
    columns <- c(columns,paste0(glb$saunames[sau],"-",1:glb$SAUpop[sau]))
  colnames(popprops) <- columns
  msy <- sapply(zoneC,"[[","MSY")
  bLML <- sapply(zoneC,"[[","bLML")
  B0 <- sapply(zoneC,"[[","B0")
  pops <- rbind(1:sum(glb$SAUpop),popprops,msy,B0,bLML)
  rownames(pops) <- c("popnum",rownames(popprops),"msy","B0","bLML")
  pops <- t(pops)
  label <- paste0(ctrl$runlabel," population properties. This is a bigtable",
                  " so, if you scroll across or down then use the topleft",
                  "  arrow to return to the home page.")
  addtable(round(pops,4),filen="popprops.csv",rundir=rundir,category="poptable",
           caption=label,big=TRUE)
  saucompdata(allcomp=condC$compdat$lfs,glb=glb,horizline=140,console=FALSE,
              rundir=rundir,ylabel="Size-Composition of Catches",
              tabname="OrigComp")
  # What size comp data is there
  palfs <- condC$compdat$palfs
  label <- paste0(ctrl$runlabel," Number of observations of numbers-at-size in ",
                  "the catch for each SAU.")
  addtable(palfs,filen="sizecompnumbers.csv",rundir=rundir,category="OrigComp",
           caption=label)
  # plot initial equilibrium size-comp by population
  Nt <- zoneDD$Nt[,1,]
  rownames(Nt) <- glb$midpts
  colnames(Nt) <- paste0(glb$saunames[glb$sauindex],"_",1:glb$numpop)
  draftnumbersatsize(rundir, glb, Nt, ssc=5)

  #properties of zoneDD
  hyrs <- glb$hyrs
  propD <- t(getzoneprops(zoneC,zoneDD,glb,year=hyrs))
  addtable(round(propD,4),"propertyDD.csv",rundir,category="zoneDD",caption=
           "Properties of zoneD after conditioning on historical catches.")
  addtable(round(t(zoneDD$harvestR[(hyrs-9):hyrs,]),4),"final_harvestR.csv",
           rundir,category="zoneDD",
           caption="Last ten years of population vs harvest rate.")
  addtable(round(zone$saudat,5),"saudat.csv",rundir,category="zoneDD",
           caption="SAU constant definitions")
  popdefs <- getlistvar(zone$zoneC,"popdef")
  addtable(t(popdefs),"popdefs.csv",rundir,category="zoneDD",
           caption="Population vs Operating model parameter definitions")
  if (dohistoric)
    condout <- plotconditioning(zoneDD,glb,zoneC,condC$histCE,condC$histCatch,
                                rundir,recdevs=condC$recdevs,console=FALSE)
    # plot predicted size-comp of catch vs observed size-comps
  if (!is.null(condC$compdat)) {
    catchN <- zoneDD$catchN
    sauCt <- popNAStosau(catchN,glb)
    lfs <- condC$compdat$lfs
    for (plotsau in 1:glb$nSAU) { # plotsau=1
      if (is.null(lfs)) {
        lfs <- NULL
       } else {
         lfs <- preparesizecomp(condC$compdat$lfs[,,plotsau],mincount=mincount,
                                deleteyears=deleteyrs[,plotsau])
         yrsize <- as.numeric(colnames(lfs))
         pickyr <- match(yrsize,condC$histyr[,"year"])
      }
      LML <- condC$histyr[pickyr,]
      plotsizecomp(rundir=rundir,incomp=lfs,SAU=glb$saunames[plotsau],lml=LML,
                   catchN=sauCt[,,plotsau],start=NA,proportion=TRUE,
                   console=FALSE)
    }
  } else {
    condout <- NULL
  }
  # plot the implied growth
  popgrowth(rundir=rundir,zoneC=zoneC,glb=glb,console=FALSE,maxage=30,
            startsize= 2.0)
  condtime <- Sys.time()
  tottime <- round((condtime - starttime),3)
  out <- list(tottime=tottime,runtime=condtime,starttime=starttime,glb=glb,
              ctrl=ctrl,zoneC=zoneC,zoneD=zoneD,zoneDD=zoneDD,
              condC=condC,condout=condout,projC=projC,production=production,
              saudat=zone$saudat,constants=zone$constants,pops=pops)
  return(out)
} # end of do_condition tottime=tottime,projtime=projtime,starttime=starttime

#' @title getCPUEssq calculates ssq between historical and conditioned cpue
#'
#' @description getCPUEssq takes in the historical CPUE and the conditioned
#'     zone's cpue and calculates the sum-of squared differences between them
#'     to assist with the conditioning. This is a greatly simplified version
#'     of plotcondCPUE, which plots a graph as well as estimating the ssq.
#'     For speed it is better to use getCPUEssq.
#'
#' @param histCE the matrix of historical cpue by SAU
#' @param saucpue the predicted cpue from the conditioning on historical data
#' @param glb the globals object
#'
#' @return a vector of length nsau containing the ssq for each SAU
#' @export
#'
#' @examples
#' print("wait on suitable internal data sets")
getCPUEssq <- function(histCE,saucpue,glb) {
  years <- as.numeric(rownames(histCE))
  hyrs <- glb$hyrnames
  pick <- match(years,hyrs)
  nsau <- glb$nSAU
  cpue <- saucpue[pick,1:nsau]
  rownames(cpue) <- years
  label <- glb$saunames
  colnames(cpue) <- label
  ssq <- numeric(nsau)
  for (sau in 1:nsau) ssq[sau] <- sum((histCE[,sau] - cpue[,sau])^2,na.rm=TRUE)
  return(ssq)
} # end of getCPUEssq

#' @title getLFlogL calculates the multi-nomial log-likelihood for sizecomp data
#'
#' @description getLFlogL when comparing the observed size-composition data from
#'     commercial catches with those rpedicted during the historical fishing
#'     period one can quantify any differences using the multi-nomial log-
#'     likelihood. This function calculates this, and the value can then be
#'     combined with am ssq value from a comparison of the predicted and observed
#'     CPUE ready to be minimized if conditioning the OM on a fisheries real
#'     data.  As always, the operating model is running at a population level
#'     while data arrives at an SAU level, so this is only an approximate
#'     approach, but it should improve over using purely the results from the
#'     sizemod R package.
#'
#' @param catchN the predicted numbers-at-size in the catch for all years
#' @param obsLFs the observed NAS, after cleaning by the preparesizecomp
#'     function
#' @param glb the globals list
#' @param wtsc the weighting given to the size-composition data
#' @param sau which sau is the subject of analysis
#'
#' @return the weighted ssq for the filtered observed size-composition data
#' @export
#'
#' @examples
#' print("wait on data sets")
getLFlogL <- function(catchN,obsLFs,glb,wtsc,sau) {
  # catchN=zoneDD$catchN; obsLFs=obsLFs; glb=out$glb;wtsc=scwt;sau=sau
  sauCN <- popNAStosau(catchN,glb)[,,sau]
  mids <- as.numeric(rownames(obsLFs))
  yrsize <- as.numeric(colnames(obsLFs))
  years <- as.numeric(colnames(sauCN))
  midpts <- glb$midpts
  pickyr <- match(yrsize,years)
  picksize <- match(mids,midpts)
  nyrs <- length(pickyr)
  LLsc <- 0
  for (yr in 1:nyrs) {
    obs <- obsLFs[,yr]
    pred <- sauCN[picksize,pickyr[yr]]
    LLsc <- LLsc + sum((obs - pred)^2,na.rm=TRUE)
  }
  return(LLsc*wtsc)
} # end of getLFlogL

#' @title getssqparts calculates the ssq for the cpue and size-comp data
#'
#' @description getssqparts is used during the conditioning of the operating
#'     model to improve the fit of the model parameters to the observed SAU
#'     scale data on cpue and size-composition. It only uses sum-of squared
#'     deviations so an arbitrary weight is given to the size-composition data,
#'     and that value should be such that changes in the ssq for the cpue are of
#'     the same order of magnitude for the size-composition data.
#'
#' @param rundir the directory in which all the files and outputs are kept
#' @param controlfile the name of the controlfile for the MSE run
#' @param calcpopC the function from the harvest strategy that is used to
#'     calculate the expected catch by population within each say each year. The
#'     annual catch by sau is known but in the MSE that needs to be distributed,
#'     in this case without error, across each population.
#' @param mincount the minimum count used in the size-composition data selection
#'     default = 100, this can be a vector of length nsau with a value for each
#'     sau if you want to give a different value to each sau.
#' @param wtsc the weighting given to the size-composition data ssq,
#'     default= 0.02. This can be a vector of length nsau if you want to give
#'     each sau a different weighting.
#'
#' @return a vector of the total ssq, the cpue ssq, and the size-comp ssq
#' @export
#'
#' @examples
#' print("wait on internal data sets")
getssqparts <- function(rundir,controlfile,calcpopC,mincount=100,wtsc=0.02) {
  # rundir=rundir; controlfile=controlfile; calcpopC=calcexpectpopC; mincount=100;wtsc=0.00001
  zone1 <- readctrlfile(rundir,infile=controlfile,verbose=TRUE)
  ctrl <- zone1$ctrl
  glb <- zone1$globals     # glb without the movement matrix
  constants <- readsaudatafile(rundir,ctrl$datafile)
  out <- setupzone(constants$constants,zone1,doproduct=FALSE,verbose=FALSE)
  zoneC <- out$zoneC
  zoneD <- out$zoneD
  glb <- out$glb             # glb now has the movement matrix
  zone1$globals <- glb
#  zoneC <- resetexB0(zoneC,zoneD) # rescale exploitB to avexplB after dynamics
  projC <- zone1$projC
  condC <- zone1$condC
  zoneDD <- dohistoricC(zoneD,zoneC,glob=glb,condC,calcpopC=calcpopC,
                        sigR=1e-08,sigB=1e-08)
  hyrs <- glb$hyrs
  sauindex <- glb$sauindex
  popB0 <- getlistvar(zoneC,"B0")
  B0 <- tapply(popB0,sauindex,sum)
  popExB0 <- getlistvar(zoneC,"ExB0")
  ExB0 <- tapply(popExB0,sauindex,sum)
  sauZone <- getsauzone(zoneDD,glb,B0=B0,ExB0=ExB0)
  ssq <- getCPUEssq(condC$histCE,sauZone$cpue,glb)
  names(ssq) <- glb$saunames
  nsau <- glb$nSAU
  LFlog <- numeric(nsau); names(LFlog) <- glb$saunames
  minsccount <- mincount
  scwt <- wtsc
  for (sau in 1:nsau) {
    if (length(mincount) > 1) minsccount <- mincount[sau]
    if (length(wtsc) > 1) scwt <- wtsc[sau]
    obsLFs <- preparesizecomp(condC$compdat$lfs[,,sau],mincount=minsccount)
    LFlog[sau] <- getLFlogL(zoneDD$catchN,obsLFs,out$glb,wtsc=scwt,sau=sau)
  }
  totssq <- ssq + LFlog
  ans <- rbind(totssq,ssq,LFlog)
  return(ans)
} # end of getssqparts


#' @title gettasdevssq calculates the SSQ between observed and predicted CPUE
#'
#' @description gettasdevssq is used when conditioning the aMSE operating model
#'     using ad hoc recruitment deviates after conditioning on AvRec. It focuses
#'     on a single SAU at a time. log transformed recdevs are used to avoid the
#'     possibility of negative recdevs. NOte this does not test for equilibrium
#'     of the initial operating model, it just assumes it reaches equilibrium.
#'
#' @param param the log transformed sequence of recdevs for the selected years
#' @param rundir the rundir for the scenario
#' @param ctrlfile the character name of the control file
#' @param calcpopC the function used to distribute catches among populations.
#'     This is from the external JurisdictionHS.R file that is 'sourced' in.
#' @param locyrs the rows within the matrix of recruitment deviates
#'     corresponding to the selected years.
#' @param sau which sau (in TAS 1 - 8) is to be worked on
#' @param obsLFs the observed length-composition of the catches for the sau,
#'     usually from condC$compdat$lfs[,,sau]
#' @param wtsc what weighting should the total Multi-Nomial ikelihood be
#'     multiplied by to scale to the range of change in thw CPUE ssq?
#'     Default=0.1
#' @param verbose should console reports be made? default = FALSE
#' @param outplot should a plot be generated for output to the webpage.
#'     default=FALSE
#' @param console should the plot go to the console. If FALSE and outplot is
#'     TRUE then a plot will go to rundir. If TRUE theplot will go to the
#'     console. DEfault=FALSE
#' @param full should the ssq components be returned together or separately.
#'     default=FALSE
#'
#' @return a scalar value which is the SSQ for the selected years of CPUE
#' @export
#'
#' @examples
#' print("wait on suitable internal data sets")
gettasdevssq <- function(param,rundir,ctrlfile,calcpopC,locyrs,sau,
                         obsLFs=obsLFs,wtsc=0.1,verbose=FALSE,outplot=FALSE,
                         console=FALSE,full=FALSE) {
# ctrlfile=controlfile;
  zone1 <- readctrlfile(rundir,infile=ctrlfile,verbose=verbose)
  zone1$condC$recdevs[locyrs,sau] <- exp(param)
  ctrl <- zone1$ctrl
  glb <- zone1$globals     # glb without the movement matrix
  saudata <- readsaudatafile(rundir,ctrl$datafile)    # make operating model
  constants <- saudata$constants
  out <- setupzone(constants,zone1,doproduct=FALSE,verbose=verbose)
  zoneC <- out$zoneC
  zoneD <- out$zoneD
  glb <- out$glb             # glb now has the movement matrix
  zone1$globals <- glb
#  zoneC <- resetexB0(zoneC,zoneD) # rescale exploitB to avexplB after dynamics
#  projC <- zone1$projC
  condC <- zone1$condC
  zoneDD <- dohistoricC(zoneD,zoneC,glob=glb,condC,calcpopC=calcpopC,
                        sigR=1e-08,sigB=1e-08)
  hyrs <- glb$hyrs
  sauindex <- glb$sauindex
  popB0 <- getlistvar(zoneC,"B0")
  B0 <- tapply(popB0,sauindex,sum)
  popExB0 <- getlistvar(zoneC,"ExB0")
  ExB0 <- tapply(popExB0,sauindex,sum)
  sauZone <- getsauzone(zoneDD,glb,B0=B0,ExB0=ExB0)
  if (outplot) {
    if (console) {
      filename <- ""
    } else {
      filename="SAU_compareCPUE.png"
    }
    ssq <- compareCPUE(condC$histCE,sauZone$cpue,glb,rundir,filen=filename)
  } else {
    ssq <- getCPUEssq(condC$histCE,sauZone$cpue,glb)
  }
  LFlog <- getLFlogL(zoneDD$catchN,obsLFs,out$glb,wtsc=wtsc,sau=sau)
  totssq <- ssq[sau] + LFlog
  if (full) {
    return(c(ssqce=ssq[sau],LFlog=LFlog,totssq=totssq))
  } else {
    return(totssq)
  }
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

#' @title optimizerecdevs improves the fit to the recdevs for a given sau
#'
#' @description optimizerecdevs an important part of optimizing the match
#'     between the dynamics observed in the fishery and those of the
#'     operating model is to adjust the annual recruitment deviates so that
#'     they improve the fit to the observed CPUE and observed numbers-at-size.
#'     optimizerecdevs attempts to do that.
#'
#' @param rundir the rundir for the scenario
#' @param sau which sau to apply this to
#' @param controlfile the name of the control file
#' @param calcpopC the funciton, from the HS file that allocates catches
#'     across each of the populations within each SAU
#' @param wtsc what weight to give to the size-composition data,default=0.1
#' @param maxiter The maximum number of iterations to run before stopping.
#'     Default = 100, one needs at least 400 to see a real difference
#' @param yearrange which recdev years should be 'conditioned'. Default=
#'     1980:2016
#' @param verbose text output to console, default=FALSE
#' @param mincount the minimum number of iterations in the solver, default=100
#' @param plottoconsole should the final plot be sent to the console = TRUE
#'     of the rundir = FALSE. default=FALSE
#' @param optimmethod which optim method to use? default='Nelder-Mead'
#'
#' @seealso{
#'     \link{adjustavrec}
#' }
#'
#' @return a scalar value which is the total SSQ for the selected years for
#'     the CPUE and sizecomps. It also alters the recdevs in the controlfile!
#' @export
#'
#' @examples
#' print("wait on data sets")
optimizerecdevs <- function(rundir,sau,controlfile,calcpopC,wtsc=0.02,
                            maxiter=100,yearrange=1980:2016,verbose=FALSE,
                            mincount=100,plottoconsole=FALSE,
                            optimmethod="Nelder-Mead") {
# rundir=rundir;sau=1;controlfile=controlfile;calcpopC=calcexpectpopC;wtsc=0.05
# maxiter=100;yearrange=1980:2016;verbose=TRUE;mincount=100;plottoconsole=FALSE
  outdevs <- getrecdevcolumn(rundir=rundir,filename=controlfile,
                             yearrange=yearrange,sau=sau)
  pickyr <- outdevs$yearrange
  param <- log(outdevs$recdevs)
  linerange <- outdevs$linerange
  starttime <- Sys.time()
  zone1 <- readctrlfile(rundir,infile=controlfile,verbose=FALSE)
  obsLFs <- preparesizecomp(zone1$condC$compdat$lfs[,,sau],mincount=mincount)
  recdevyrs <- as.numeric(rownames(zone1$condC$recdevs))
  locyrs <- match(pickyr,recdevyrs)
  gettasdevssq(param,rundir,ctrlfile=controlfile,calcpopC=calcpopC,
               locyrs,sau,obsLFs=obsLFs,wtsc=wtsc,outplot=TRUE,console=TRUE)
  starttime <- Sys.time()
  ans <- optim(param,gettasdevssq,method=optimmethod,rundir=rundir,
               ctrlfile=controlfile,calcpopC=calcpopC,
               locyrs=locyrs,sau=sau,obsLFs=obsLFs,wtsc=wtsc,outplot=FALSE,
               control=list(maxit=maxiter,trace=4,reltol=1e-06))
  outfit(ans)
  print(Sys.time() - starttime)
  changecolumn(rundir=rundir,filename=controlfile,linerange=linerange,
               column=(sau+1),newvect=round(exp(ans$par),5),verbose=verbose)
  ssq <- gettasdevssq(ans$par,rundir,ctrlfile=controlfile,
                      calcpopC=calcpopC,locyrs,sau,obsLFs=obsLFs,
                      wtsc=wtsc,outplot=TRUE,console=plottoconsole,full=TRUE)
  return(ssq)
} # end of optimizerecdevs


#' @title plotcondCPUE plots the historical cpue against conditioned cpue
#'
#' @description plotcondCPUE generates a plot of the historical cpue and compares
#'     it with the predicted cpue from the conditioning. The predicted is black
#'     and the observed is green. This is primarily there are an aid to
#'     conditioning the operating model. It also calculates the simple sum of
#'     squared differences between the observed and predicted, again to aid in
#'     the conditioning.
#'
#' @param histCE the matrix of historical cpue by SAU
#' @param saucpue the predicted cpue from the conditioning on historical data
#' @param glb the globals object
#' @param rundir the rundir for the given scenario
#' @param filen the file name of the plot if it is to be saved
#'
#' @return a vector of length nsau containing the ssq for each SAU
#' @export
#'
#' @examples
#' print("wait on suitable internal data sets")
plotcondCPUE <- function(histCE,saucpue,glb,rundir,filen="") {
  if (nchar(filen) > 0) filen <- filenametopath(rundir,filen)
  years <- as.numeric(rownames(histCE))
  hyrs <- glb$hyrnames
  pick <- match(years,hyrs)
  nsau <- glb$nSAU
  cpue <- saucpue[pick,1:nsau]
  rownames(cpue) <- years
  label <- glb$saunames
  colnames(cpue) <- label
  ssq <- numeric(nsau)
  doplots=pickbound(nsau)
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,verbose=FALSE)
  parset(plots=doplots,margin=c(0.3,0.3,0.05,0.05),outmargin=c(0,1,0,0),
         byrow=FALSE)
  for (sau in 1:nsau) {
    ssq[sau] <- sum((histCE[,sau] - cpue[,sau])^2,na.rm=TRUE)
    ymax <- getmax(c(cpue[,sau],histCE[,sau]))
    plot(years,cpue[,sau],type="l",lwd=2,col=1,xlab="",ylab="",panel.first=grid(),
         ylim=c(0,ymax),yaxs="i")
    lines(years,histCE[,sau],lwd=2,col=3)
    lab2 <- paste0(label[sau],"  ",round(ssq[sau],1))
    mtext(lab2,side=1,line=-1.3,cex=1.25)
  }
  mtext("CPUE",side=2,line=-0.3,outer=TRUE,cex=1.25)
  if (nchar(filen) > 0) {
    caption <- "Comparison of Observed CPUE with that predicted by operating model."
    addplot(filen,rundir=rundir,category="condition",caption)
  }
  invisible(ssq)
} # end of plotcondCPUE


#' @title printline literally prints a selected line from a text file
#'
#' @description printline is used to check that one has selected the
#'     correct line for modification inside a datafile. This is done when
#'     conditioning the MSE operating model on the parameter AvRec
#'
#' @param rundir the scenario directory holding all files for a scenario
#' @param datafile the name of the datafile for the scenario, from the ctrl
#'     object
#' @param linenumber which linenumber to print, default =29
#'
#' @return nothing but it does print a line to the console
#' @export
#'
#' @examples
#' print("wait on suitable example")
printline <- function(rundir, datafile, linenumber=29) {
  dat <- readLines(filenametopath(rundir,datafile))
  if (length(dat) <= linenumber) {
    print(dat[linenumber])
  } else {
    print("You are trying to print past the length of the input file!")
  }
} #end of printline


#' @title prodforpop plots the productivity characteristics for a population
#'
#' @description prodforpop plots the productivity characteristics for a
#'     population. This is especially useful for a detailed consideration of
#'     an SAU's dynamics. This also allows one to check that the harvest rate
#'     range used to estimate production runs past the MSY level.
#'
#' @param rundir the full path to the scenario directory
#' @param inprod the production matrix for a single population, for example,
#'     out$production[,,2] is the second population's production properties.
#' @param pop which population is being examined? A single integer
#' @param console should the plot be saved to a png file or go to the console.
#'     the default = TRUE plotting to theconsole
#'
#' @return nothing but it does generate a plot
#' @export
#'
#' @examples
#' # syntax prodforpop(rundir=rundir,inprod=produciton[,,2],pop=2,console=FALSE)
#' # rundir=""; inprod=out$production[,,2]; pop=2; console=TRUE
prodforpop <- function(rundir,inprod,pop,console=TRUE) {
  yield <- inprod[,"Catch"]
  spb <- inprod[,"MatB"]
  Ht <- rownames(inprod)
  deplet <- inprod[,"Deplet"]
  pickmsy <- which.max(yield)
  msy <- yield[pickmsy]
  maxy <- getmax(yield)
  filen <- ""
  if (!console) {
    label <- paste0("production_plots_for_population_",pop,".png")
    filen <- pathtopath(rundir,label)
  }
  plotprep(width=7,height=6,newdev=FALSE,filename=filen,verbose=FALSE)
  parset(plots=c(3,2),cex=0.9,outmargin=c(1,0.25,0,0))
  plot(spb,yield,type="l",lwd=2,col=1,xlab="Spawning Biomass t",
       ylab="",panel.first = grid(),
       ylim=c(0,maxy),yaxs="i")
  abline(v=spb[pickmsy],col=2,lwd=2)
  plot(Ht,spb,type="l",lwd=2,xlab="Annual Harvest Rate",
       ylab="Spawning Biomass t",panel.first = grid(),
       ylim=c(0,getmax(spb)),yaxs="i")
  abline(h=spb[pickmsy],col=2,lwd=2)
  abline(v=Ht[pickmsy],col=2,lwd=2)
  plot(Ht,yield,type="l",lwd=2,col=1,xlab="Annual Harvest Rate",
       ylab="",panel.first = grid(),
       ylim=c(0,maxy),yaxs="i")
  abline(v=Ht[pickmsy],col=2,lwd=2)
  plot(spb,deplet,type="l",lwd=2,ylab="Total Depletion Level",
       xlab="Spawning Biomass t",panel.first = grid(),
       ylim=c(0,1.05),yaxs="i")
  abline(h=deplet[pickmsy],col=2,lwd=2)
  abline(v=spb[pickmsy],col=2,lwd=2)
  plot(deplet,yield,type="l",lwd=2,col=1,xlab="Total Depletion Level",
       ylab="",panel.first = grid())
  abline(v=deplet[pickmsy],col=2,lwd=2)
  plot(Ht,deplet,type="l",lwd=2,col=1,xlab="Annual Harvest Rate",
       ylab="Total Depletion Level",panel.first = grid(),
       ylim=c(0,1.05),yaxs="i")
  abline(h=deplet[pickmsy],col=2,lwd=2)
  abline(v=Ht[pickmsy],col=2,lwd=2)
  mtext(paste0("Population ",pop),side=1,line=-0.1,outer=TRUE,cex=1.1)
  mtext("Production t",side=2,line=-1,outer=TRUE,cex=1.1)
  if (!console) {
    caption <- "The production curves for the population."
    addplot(filen,rundir=rundir,category="popprod",caption)
  }
} # end of prodforpop

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
#' @param controlfile the name of the control file, no path
#' @param datafile the name of the data file, no path
#' @param linenum the linenumber containing the AvRec vector
#' @param calcpopC the HS function that calculates the aspirational catches
#'     for each SAU
#' @param extra the vector of SAU AvRecs minus the selected SAU value
#' @param picksau which SAU is to be worked on as in 1:nsau
#' @param nsau the total number of SAU
#' @param outplot shoudla plot be generated for output to the webpage.
#'     default=FALSE
#'
#' @return the SSQ for the selected SAU in a comparison of predicted and
#'     observed CPUE
#' @export
#'
#' @examples
#' print("Wait on suitable internal data files")
#' #  param=param;rundir=rundir;controlfile=controlfile;datafile=datafile;linenum=29,
#' #  calcpopC=calcexpectpopC;extra=extra;picksau=sau;nsau=nsau
sauavrecssq <- function(param,rundir,controlfile,datafile,linenum,
                        calcpopC,extra,picksau,nsau,outplot=FALSE) {
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
  changeline(rundir,datafile,linenum,replacetxt)
  zone <- makeequilzone(rundir,controlfile,doproduct=FALSE,verbose=FALSE)
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
  if (outplot) {
    filename=paste0("compareCPUE_",picksau,"_AvRec_cond.png")
    ssq <- plotcondCPUE(condC$histCE,sauZone$cpue,glb,rundir,filen=filename)
  } else {
    ssq <- getCPUEssq(condC$histCE,sauZone$cpue,glb)
  }
  return(ssq[picksau])
} # end of sauavrecssq
