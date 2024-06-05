



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
#' @param openfile should the HTML outputs be opened, default=TRUE
#'
#' @return nothing but it does plot a graph
#' @export
#'
#' @examples
#' print("wait on data sets")
checksizecompdata <- function(rundir,controlfile,verbose=TRUE,openfile=TRUE) {
  zone1 <- readctrlfile(rundir,infile=controlfile,verbose=verbose)
  setuphtml(rundir)
  compdat <- zone1$condC$compdat$lfs
  if (is.null(compdat))
    stop(cat("No size-composition file found in copntrol file \n"))
  palfs <- zone1$condC$compdat$palfs
  saunames <- zone1$SAUnames
  nsau <- length(saunames)
  histyr <- zone1$condC$histyr
  sauLML <- zone1$globals$sauLML
  addtable(palfs,"sizecomp_obs_by_year_by_SAU.csv",rundir,
           category="predictedcatchN",
           caption="Number of sizecomp observations by year and SAU.")
  for (plotsau in 1:nsau) {
    lfs <- preparesizecomp(compdat[,,plotsau],mincount=0,deleteyears=0)
    yrsize <- as.numeric(colnames(lfs))
    pickyr <- match(yrsize,histyr[,"year"])
#    if (sauLML) {
      LML <- histyr[pickyr,]
    # } else {
    #   LML <- histyr[pickyr,2]
    # }
    plotsizecomp(rundir=rundir,incomp=lfs,SAU=saunames[plotsau],lml=LML,
                 catchN=NULL,start=NA,proportion=TRUE,
                 console=FALSE)
  }
  replist <- NULL
  projy <- 0
  runnotes <- c(controlfile,paste0("SAU = ",nsau))
  make_html(replist = replist,  rundir = rundir,
            controlfile=controlfile, datafile=zone1$ctrl$datafile,
            width = 500, openfile = openfile,  runnotes = runnotes,
            verbose = verbose, packagename = "aMSE",  htmlname = "sizecomp")
  if (verbose) cat("finished  \n")
} # end of checksizecompdata

#' @title draftnumbersatsize plots details of initial numbers-at-size
#'
#' @description draftnumbersatsize plots up the initial unfished numbers-
#'     at-size distribution, including the totals for each sau, and the
#'     numbers-at-size for each population within each sau.
#'
#' @param rundir the results directory output files will be sent here
#' @param glb the globals list
#' @param Nt the numbers-at-size at equilibrium = zoneD$Nt[,1,]
#' @param ssc index for starting size class. thus 1 = 2, 2 = 4, 5 = 10, etc.
#'     default = 5 for it plots size classes from 10mm up
#'
#' @return nothing but it does add numpops + 1 plots to the results directory
#' @export
#'
#' @examples
#' print("this will be quite long when I get to it")
draftnumbersatsize <- function(rundir, glb, Nt, ssc=5) {
  mids <- glb$midpts
  nc <- glb$Nclass
  nsau <- glb$nSAU
  saunames <- glb$saunames
  Nt <- Nt/1000.0
  filen <- pathtopath(rundir,"Initial_N-at-Size_by_SAU.png")
  plotprep(width=7,height=6,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  parset(plots=pickbound(nsau))
  for (plotsau in 1:nsau) {
    pickcol <-  which(glb$sauindex == plotsau)
    sNt <-
      Ntt <- rowSums(Nt[,pickcol],na.rm=TRUE)
    xlabel <- paste0("Shell Length mm (",mids[ssc]," - 210mm)")
    plot(mids[ssc:nc],Ntt[ssc:nc],type="l",lwd=2,xlab=xlabel,
         ylab="Numbers-at_size '000s",panel.first=grid())
    mtext(saunames[plotsau],side=3,line=-1.25,outer=FALSE,cex=1.25,font=7)
  }
  caption <- paste0("The numbers-at-size for each SAU separately. Recruitment",
                    " numbers are omitted for clarity.")
  addplot(filen,rundir=rundir,category="NumSize",caption)
  # Now plot individual populations
  for (plotsau in 1:nsau) { # plotsau=1
    filename <- paste0("Initial_N-at-Size_by_pop_for_",saunames[plotsau],".png")
    filen <- pathtopath(rundir,filename)
    pickcol <-  which(glb$sauindex == plotsau)
    npop <- length(pickcol)
    sNt <- Nt[,pickcol]
    xlabel <- paste0("Shell Length mm (",mids[ssc]," - 210mm)")
    plotprep(width=12,height=10,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
    parset(plots=pickbound(npop),outmargin=c(1,1,0,0),byrow=FALSE)
    for (pop in 1:npop) {
      plot(mids[ssc:nc],Ntt[ssc:nc],type="l",lwd=2,xlab="",
           ylab=colnames(sNt)[pop],panel.first=grid())
    }
    mtext(xlabel,side=1,line=-0.2,outer=TRUE,cex=1.25,font=7)
    mtext("Numbers-at_size '000s",side=2,line=-0.2,outer=TRUE,cex=1.25,font=7)
    caption <- paste0("The numbers-at-size for each SAU and population ",
                    " separately. Recruitment numbers are omitted for clarity.")
    addplot(filen,rundir=rundir,category="NumSize",caption)
  }
} # end of draftnumbersatsize

#' @title getdraftpops generates a draft copy of the population properties
#'
#' @description getdraftpops is used to generate a first draft of the population
#'     properties from the sau data using a standard sau based control and data
#'     file. In am sau based run of aMSE the input properties for each sau are
#'     expanded in a random manner across the populations, all that is except
#'     the AvRec, which is determined by the use of proprec in the datafile.
#'     However, all other properties are assigned randomly. The outcome is
#'     characterized in terms of each population's productivity, which can then
#'     be edited to more closely sorrespond to the productivity as observed in
#'     the GPS data.
#'
#' @param rundir the directory containing the control csv files. It
#'     can/will also act to store results in a manner that will allow them
#'     to be displayed using makehtml.
#' @param ctrlfile the main file that controls the particular run. It contains
#'     the name of the data file that is used to biologically condition the
#'     numpop populations
#' @param verbose Should progress comments be printed to console, default=TRUE
#' @param uplim defines the upper limit of harvest used when estimating the
#'     productivity (also important when initial depletion is not 1.0). The
#'     default = 0.35
#' @param incH defines the interval between H steps when estimating productivity
#'     default = 0.005
#'
#' @return a list of the draft population properties needed by makezoneC,
#'     an object containing zoneC, zoneD, product, and glb, and the condC
#' @export
#'
#' @examples
#' print("wait on example") # verbose=TRUE;uplim=0.35;incH=0.005
getdraftpops <- function(rundir,ctrlfile,verbose=TRUE,uplim=0.35,incH=0.005) {
  zone1 <- readctrlfile(rundir,infile=ctrlfile,verbose=verbose)
  ctrl <- zone1$ctrl
  glb <- zone1$globals     # glb without the movement matrix
  bysau <- zone1$ctrl$bysau
  opar <- NULL
  parsin <- zone1$condC$parsin
  if (parsin) opar <- as.matrix(zone1$condC$optpars)
  if (is.null(bysau)) bysau <- 0
  if (bysau) {
    saudata <- readsaudatafile(rundir,ctrl$datafile,optpar=opar)
    constants <- saudata$constants
    saudat <- saudata$saudat
    zone1$condC$poprec <- saudata$poprec
  } else {
    constants <- readpopdatafile(rundir,ctrl$datafile)
    saudat <- constants
  }
  if (verbose) cat("Files read, now making zone \n")
  out <- setupzone(constants,zone1,doproduct=TRUE,uplim=uplim,inc=incH,
                   verbose=verbose) # make operating model
  popprops <- sapply(out$zoneC,"[[","popdef")
  columns <- NULL
  for (sau in 1:glb$nSAU)
    columns <- c(columns,paste0(glb$saunames[sau],"-",1:glb$SAUpop[sau]))
  colnames(popprops) <- columns
  popprops2 <- rbind(popnum=1:glb$numpop,popprops)
  return(invisible(list(pops=popprops2,out=out,condC=zone1$condC,ctrl=ctrl)))
} # end of getdraftpops



