



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
  plotprep(width=7,height=6,newdev=TRUE,filename=filen,cex=0.9,verbose=FALSE)
  parset(plots=pickbound(nsau))
  for (plotsau in 1:nsau) { # plotsau=1
    pickcol <-  which(glb$sauindex == plotsau)
    if (length(pickcol) > 1) {
      Ntt <- rowSums(Nt[,pickcol],na.rm=TRUE)
    } else {
      Ntt <- Nt[,pickcol]
    }
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


