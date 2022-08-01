

#' @title comparevar generates the quantiles for each of a set of input scenarios
#'
#' @description comparevar is used when post-processing results and comparing
#'     different scenarios. For a given variable within the projected dynamics,
#'     the function estimates the quantiles for each of the set of input
#'     scenarios. It takes the full timeline of dynamics and outputs just the
#'     projections and the quantiles of those projections. The variables it can
#'     work with include: matureB, exploitB, midyrexpB, catch, acatch, harvestR,
#'     cpue, recruit, deplsB, and depleB.
#'
#' @param dyn this is a list of the out$sauout$zonePsau produced by each
#'     scenario. The out object can be from a saved RData file from
#'     each scenario
#' @param glbc a list of the global objects from each scenario being compared
#' @param scenes a list of the out$ctrl$runlabel from each scenario
#' @param var what variable from the dynamics to summarize, valid names include
#'     catch, acatch, cpue, harvestR, desplsB, depleB, recruit, matureB, and
#'     exploitB, default = 'cpue'
#'
#' @seealso {
#'   \link{plotscene}
#' }
#'
#' @return a list of the quantiles of the input var for all sau, with each
#'     scenario being a different component of the quantscen list, and a list of
#'     three dimensional arrays of the actual values of var in the varc list.
#' @export
#'
#' @examples
#' print("wait on internal datasets")
comparevar <- function(dyn,glbc,scenes,var="cpue") {
  # dyn=dyn; glbc=glbc; scenes=scenes; var="catch"
  nscen <- length(dyn)
  glb1 <- glbc[[1]]
  nsau <- glb1$nSAU
  reps <- glb1$reps
  hyrs <- glb1$hyrs
  pyrs <- glb1$pyrs
  projy <- (hyrs+1):(hyrs+pyrs)
  pyrnames <- glb1$pyrnames
  saunames <- glb1$saunames
  availvar <- names(dyn[[1]])
  varc <- vector(mode="list",length=nscen)
  names(varc) <- scenes
  pickV <- grep(var,availvar)[1]
  if (is.na(pickV)) stop(cat(var," is not a variable in the dynamics \n"))
  for (i in 1:nscen) varc[[i]] <- dyn[[i]][[pickV]][projy,,]
  quantres <- vector(mode="list",length=nsau)
  names(quantres) <- saunames
  quantscen <- vector(mode="list",length=nscen)
  names(quantscen) <- scenes
  for (i in 1:nscen) { # i = 1
    proj <- varc[[i]]
    for (j in 1:nsau) quantres[[j]] <-  apply(proj[,j,],1,quants)
    quantscen[[i]] <- quantres
  }
  return(list(quantscen=quantscen,varc=varc))
} # end of comparevar

#' @title plotscene literally plots up the output from comparevar
#'
#' @description plotscene takes the output from comparevar and plots the
#'     quantiles relative to each other.
#'
#' @param scenquant a list of the quantiles for the given variable from each
#'     scenario
#' @param glbc  a list of the global objects from each scenario being compared
#' @param var  what variable from the dynamics to summarize, valid names include
#'     catch, acatch, cpue, harvestR, desplsB, depleB, recruit, matureB, and
#'     exploitB, default = 'cpue'
#' @param ymin allows one to set the lower limit to the yaxis, default = 0
#' @param filen a filename for use if saving the output to said file,
#'     default = '', which implies the plot goes to the console.
#' @param legloc legend location, default='topleft'
#' @param legplot which plot should contain the legend, default=1
#' @param qnt which quantile to plot? defauilt = 90 percent
#'
#' @return nothing but it does generate a plot with nsau panels
#' @export
#'
#' @examples
#'  print("wait on internal datasets")
plotscene <- function(scenquant,glbc,var="cpue",ymin=0,filen="",
                      legloc="topleft",legplot=1,qnt=90) {
  scenes <- names(scenquant)
  nscen <- length(scenes)
  saunames <- names(scenquant[[1]])
  nsau <- length(saunames)
  qval <- scenquant[[1]]
  yrs <- as.numeric(colnames(qval[[1]]))
  maxy <- matrix(0,nrow=nscen,ncol=nsau)
  for (i in 1:nscen) {
    qval <- scenquant[[i]]
    for (j in 1:nsau) maxy[i,j] <- max(qval[[j]])
  }
  ymax <- apply(maxy,2,max)
  doplots=pickbound(nsau)
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,verbose=FALSE)
  parset(plots=doplots,margin=c(0.3,0.4,0.05,0.05),outmargin=c(0,1,0,0),
         byrow=FALSE)
  for (j in 1:nsau) { # j = 1; i = 2
    y50 <- scenquant[[1]][[j]][3,]
    plot(yrs,y50,type="l",lwd=3,col=1,ylim=c(ymin,ymax[j]),panel.first=grid(),
         ylab=saunames[j],xlab="")
    if (qnt == 90) {
      lines(yrs,scenquant[[1]][[j]][2,],lwd=1,col=1)
      lines(yrs,scenquant[[1]][[j]][4,],lwd=1,col=1)
    } else {
      lines(yrs,scenquant[[1]][[j]][1,],lwd=1,col=1)
      lines(yrs,scenquant[[1]][[j]][5,],lwd=1,col=1)
    }
    if (legplot == j) {
      legend(legloc,legend=scenes,col=1:nscen,lwd=3,bty="n",cex=1.0)
    }
    for (i in 2:nscen) {
      lines(yrs,scenquant[[i]][[j]][3,],lwd=3,col=i)
      if (qnt == 90) {
        lines(yrs,scenquant[[i]][[j]][2,],lwd=1,col=i)
        lines(yrs,scenquant[[i]][[j]][4,],lwd=1,col=i)
      } else {
        lines(yrs,scenquant[[i]][[j]][1,],lwd=1,col=i)
        lines(yrs,scenquant[[i]][[j]][5,],lwd=1,col=i)
      }
    }
  }
  mtext(var,side=2,line=-0.2,outer=TRUE,cex=1.2)
#  if (nchar(filen) > 0) dev.off()
} # end of plotscene


#' @title comparedynamics plots medians and quantiles of dynamic variables
#'
#' @description comparedynamics plots, for each sau, the medians and 90th
#'     quantiles of the distributions of the reps for each scenario for four
#'     dynamic variables: cpue, catch (actual catches not aspirational),
#'     mature biomass, and harvets rate.
#'
#' @param rundir to location for all results of the comparisons
#' @param dyn a list of the dynamics objects from each scenario's out object
#' @param glbc a list of the globals objects from each scenario
#' @param scenes a vector of the names for each scenario
#'
#' @return nothing but it does add four plots to rundir and four lines to the
#'     resultTable.csv file in rundir ready for a makecompareoutput
#' @export
#'
#' @examples
#' \dontrun{
#'   print("wait on an example")
#'
#' }
comparedynamics <- function(rundir,dyn,glbc,scenes) {
  invar <- "cpue"
  compscen <- comparevar(dyn,glbc,scenes,var=invar)
  filename <- filenametopath(rundir,"compare_cpue_scenes.png")
  plotscene(compscen$quantscen,glbc=glbc,var=invar,ymin=0,filen=filename,
            legloc="bottomleft",legplot=1)
  caption <- paste0("The projected CPUE dynamics for each SAU.")
  addplot(filen=filename,rundir=rundir,category="Dynamics",caption)

  invar <- "catch"
  compscen <- comparevar(dyn,glbc,scenes,var=invar)
  filename <- filenametopath(rundir,"compare_catch_scenes.png")
  plotscene(compscen$quantscen,glbc=glbc,var=invar,ymin=0,filen=filename,
            legloc="topleft",legplot=1)
  caption <- paste0("The projected actual catch dynamics for each SAU.")
  addplot(filen=filename,rundir=rundir,category="Dynamics",caption)

  invar <- "matureB"
  compscen <- comparevar(dyn,glbc,scenes,var=invar)
  filename <- filenametopath(rundir,"compare_matureB_scenes.png")
  plotscene(compscen$quantscen,glbc=glbc,var=invar,ymin=0,filen=filename,
            legloc="topleft",legplot=1)
  caption <- paste0("The projected mature biomass dynamics for each SAU.")
  addplot(filen=filename,rundir=rundir,category="Dynamics",caption)

  invar <- "harvestR"
  compscen <- comparevar(dyn,glbc,scenes,var=invar)
  filename <- filenametopath(rundir,"compare_harvestR_scenes.png")
  plotscene(compscen$quantscen,glbc=glbc,var=invar,ymin=0,filen=filename,
            legloc="topleft",legplot=1)
  caption <- paste0("The projected harvest rate dynamics for each SAU.")
  addplot(filen=filename,rundir=rundir,category="Dynamics",caption)
} # end of comparedynamics

#' @title makecompareoutput generates HTML files when comnparing scenarios
#'
#' @description makecompareoutput simplifies taking the output from when
#'     comparing scenarios and producing the HTML files used to display all the
#'     results in rundir. Note that 'runnotes' is a vector of character strings
#'     made up using paste0.
#'
#' @param rundir the full path to the directory holding the all result files
#' @param glbc a list of the globals objects from each scenario
#' @param scenes a vector of the names for each scenario
#' @param postdir the name of the directory holding the results, also used to
#'     name the internal webpage
#' @param openfile should the website be opened automatically? default=TRUE
#' @param verbose should details of producing the HTML files be written to the
#'     console?  default=FALSE
#'
#' @return nothing but it does generate a set of HTML files in rundir. If
#'     verbose=TRUE it also writes text to the console
#' @export
#'
#' @examples
#' print("wait on internal data-sets")
makecompareoutput <- function(rundir,glbc,scenes,postdir,
                              openfile=TRUE,verbose=FALSE) {
  replist <- list(starttime=as.character(Sys.time()),
                  endtime="")
  nscene <- length(scenes)
  reps <- NULL
  for (i in 1:nscene) reps <- paste0(reps,glbc[[i]]$reps,"  ")
  runnotes <- c("Scenario Comparisons",paste0("RunTime = ",0),
                paste0("replicates = ",reps),
                paste0("years projected = ",glbc[[1]]$pyrs),
                paste0("Populations = ",glbc[[1]]$numpop),
                paste0("SAU = ",glbc[[1]]$nSAU))

  make_html(replist = replist,  rundir = rundir,
            controlfile=NULL, datafile=NULL, hsfile="TasHS Package",
            width = 500, openfile = TRUE,  runnotes = runnotes,
            verbose = verbose, packagename = "aMSE",  htmlname = postdir)
  if (verbose) cat("finished  \n")
}

#' @title tabluateproductivity produces tables of the MSY and related statistics
#'
#' @description tabluateproductivity writes out csv files for each scenario
#'     that contain the B0, Bmsy, MSY, Depmsy, and CEmsy for each sau
#'
#' @param rundir the full path to the directory holding the all result files
#' @param prods the sauprod objects from each scenario in a list
#' @param scenes a vector of the names for each scenario
#'
#' @return nothing but it does add nscenario csv files as tables to the rundir
#' @export
#'
#' @examples
#' print("wait on data sets")
tabluateproductivity <- function(rundir,prods,scenes) {
  nscen <- length(scenes)
  for (i in 1:nscen) {
    filen <- paste0("productivity_",scenes[i],".csv")
    caption <- paste0("SAU productivity statistics for ",scenes[i],".")
    addtable(prods[[i]],filen,rundir=rundir,category="productivity",caption)
  }
} # end of tabulateproductivity


