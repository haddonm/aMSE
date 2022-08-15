#' @title boxbysau generates boxplots of sau x aavc or sum5 or sum10
#'
#' @description boxbysau is used to generate boxplots of the aavc, or sum5,
#'     of sum10 for each sau, in each of the scenarios. The number of years over
#'     which the aavc is calculated can be modified using the aavcyrs argument.
#'
#' @param rundir the directory in which comparisons are being made. It is best
#'     to be a separate directory form any particular scenario.
#' @param hspm a list of the harvest strategy performance measures for actual
#'     catches as calculated by calccatchHSPM
#' @param glbc a list of the globals objects from each scenario
#' @param scenes the names of the different scenarios being compared
#' @param compvar which variable should be compared. The options are 'aavc',
#'      'sum5', and 'sum10'. default = 'aavc'
#' @param filen the filename for saving the plot, the default="", which sends
#'     the plot to the console.
#' @param aavcyrs the number of projection years to use for the aavc calculation
#'     the default = 10
#' @param maxval the maximum y-axis value for the aavc, default=0, which means
#'     the maximum for each boxplot is taken from the data for that plot.
#'
#' @seealso{
#'    \link{getprojyrs}, \link{getprojyraavc},  \link{catchHSPM}
#' }
#'
#' @return a list of the boxplot statistics for the aavc x number of sau x the
#'     number of scenarios compared
#' @export
#'
#' @examples
#' print("wait on data sets")
boxbysau <- function(rundir,hspm,glbc,scenes,compvar="aavc",filen="",aavcyrs=10,
                      maxval=0) {
  #  compvar="sum10"; hspm=outcatchHSPM
  nscen <- length(scenes)
  pickV <- grep(compvar,names(hspm))
  boxvar <- hspm[[pickV]]
  nsau <- glbc[[1]]$nSAU
  saunames <- glbc[[1]]$saunames
  boxresult <- vector(mode="list",length=nsau)
  names(boxresult) <- saunames
  reps <- glbc[[1]]$reps
  boxdat <- matrix(0,nrow=reps,ncol=nscen,dimnames=list(1:reps,scenes))
  if (nchar(filen) > 0) {
    filen <- filenametopath(rundir,filen)
    caption <- paste0("Boxplots of ",compvar," for each SAU.")
  }
  plotprep(width=8, height=8,newdev=FALSE,filename = filen,verbose=FALSE)
  parset(plots=pickbound(nsau),margin=c(0.3,0.4,0.05,0.1),outmargin=c(1,0,0.25,0),
         byrow=FALSE)
  for (sau in 1:nsau) { # sau=1; i = 1
    for (i in 1:nscen) boxdat[,i] <- boxvar[[i]][,sau]
    if (maxval > 0) {
      maxy <- maxval
    } else {
      maxy <- getmax(boxdat)
    }
    boxresult[[sau]] <- boxplot(boxdat,ylim=c(0,maxy),yaxs="i",
                                ylab=saunames[sau])
    abline(h=boxresult[[sau]]$stats[3,1],lwd=2,lty=2,col=2)
  }
  label <- switch(compvar,
                  "aavc" = "Average Annual Variation in Catches",
                  "sum5" = "Sum of Catches over first 5 Years",
                  "sum10" = "Sum of Catches over first 10 Years")
  mtext(label,side=1,outer=TRUE,line=-0.3,cex=1.0)
  if (nchar(filen) > 0)
    addplot(filen=filen,rundir=rundir,category="catches",caption)
} # end of boxbysau


#' @title calccatchHSPM calculates all the HSPMs for the actual catch data
#'
#' @description calccatchHSPM calculates all the Harvest Strategy Performance
#'     Measures for the actual catch data. The number of years over which the
#'     aavc is calculated can be modified using the aavcyrs argument. The output
#'     data can be analysed separately or plotted (eg as boxplots)
#'
#' @param dyn a list of the dynamics objects from each scenario's out object
#' @param glbc a list of the globals objects from each scenario
#' @param scenes the names of the different scenarios being compared
#' @param aavcyrs the number of projection years to use for the aavc calculation
#'     the default = 10
#'
#' @seealso{
#'    \link{getprojyrs}, \link{getprojyraavc},  \link{catchHSPM}
#' }
#'
#' @return a list of  the aavc x sau x scene, and the sum5 and sum10 for the
#'     same sau and scenes.
#' @export
#'
#' @examples
#' print("wait on data sets")
calccatchHSPM <- function(dyn,glbc,scenes,aavcyrs=10) {
  # dyn=dyn;glbc=glbc;scenes=scenes;aavcyrs=10
  nscen <- length(scenes)
  catch <- makelist(scenes)
  catches <- makelist(scenes)
  aavc <- makelist(scenes)
  sum5 <- makelist(scenes)
  sum10 <- makelist(scenes)
  saunames <- glbc[[1]]$saunames
  for (i in 1:nscen) {   # i = 1
    catch[[i]] <- dyn[[i]]$catch
    catches[[i]] <- getprojyrs(catch[[i]],glbc[[i]]$hyrs,glbc[[i]]$pyrs,
                             startyr=glbc[[i]]$hyrs)
    aavc[[i]] <- getprojyraavc((catches[[i]][1:aavcyrs,,]),glbc[[i]])
    sum5[[i]] <- getprojyrC(catch[[i]],glbc[[i]],period=5)
    sum10[[i]] <- getprojyrC(catch[[i]],glbc[[i]],period=10)
  }
  return(invisible(list(aavc=aavc,sum5=sum5,sum10=sum10)))
} # end of calccatchHSPM


#' @title catchbpstats saves the boxplot statistics and adds to MSE webpage
#'
#' @description catchbpstats takes the results from plotting the boxplots of the
#'     catch HSPMs aavc, sum5, and sum10 and tabulates the statistical
#'     properties of the boxplots, ie, the lower and upper whiskers, the
#'     inter-quartile bounds and the median. These are all output into a single
#'     table. The aavc, sum5, and sum10 results for each scenario are grouped
#'     together in the final data.frame.
#'
#' @param rundir the directory into which all the comparison results are placed
#' @param outtab the output from catchHSPM
#'
#' @seealso  \link{catchHSPM}
#'
#' @return  invisibly returns the matrix version of the input list. It also
#'     saves a csv file and adds to the webpage
#' @export
#'
#' @examples
#' print("wait on example data sets")
catchbpstats <- function(rundir,outtab) {
  # translate the outtab list into a dataframe
  # rundir=rundir; outtab=outtab
  label <- names(outtab)
  nlab <- length(label)
  col1 <- NULL
  pick <- grep("aavc",label)
  npick <- length(pick)
  for (i in 1:npick) col1 <- c(col1,rep(label[pick[i]],5))
  pick <- grep("sum5",label)
  npick <- length(pick)
  for (i in 1:npick) col1 <- c(col1,rep(label[pick[i]],5))
  pick <- grep("sum10",label)
  npick <- length(pick)
  for (i in 1:npick) col1 <- c(col1,rep(label[pick[i]],5))
  numrow <- length(col1)
  columns <- c("scenarios","bpstat",colnames(outtab[[1]]))
  numcol <- length(columns)
  col2 <- rep(rownames(outtab[[1]]),nlab)
  result <- as.data.frame(matrix(0,nrow=numrow,ncol=numcol,
                                 dimnames=list(1:numrow,columns)))
  result[,1] <- col1
  result[,2] <- col2
  begin <- 1
  finish <- 5
  pick <- grep("aavc",label)
  npick <- length(pick)
  for (i in 1:npick) {
    result[begin:finish,3:numcol] <- outtab[[pick[i]]]
    begin <- begin + 5
    finish <- finish + 5
  }
  pick <- grep("sum5",label)
  npick <- length(pick)
  for (i in 1:npick) {
    result[begin:finish,3:numcol] <- outtab[[pick[i]]]
    begin <- begin + 5
    finish <- finish + 5
  }
  pick <- grep("sum10",label)
  npick <- length(pick)
  for (i in 1:npick) {
    result[begin:finish,3:numcol] <- outtab[[pick[i]]]
    begin <- begin + 5
    finish <- finish + 5
  }
  filen <- paste0("catch_boxplot_stats.csv")
  caption <- "Boxplot statistics for all scenarios."
  addtable(result,filen,rundir=rundir,category="catches",caption)
  return(invisible(result))
} # end of catchbpstats

#' @title catchHSPM generates boxplots of catch HSPMs and output statistics
#'
#' @description catchHSPM is used to generate boxplots of the aavc, the sum5,
#'     and the sum10 for projected catches. The number of years over which the
#'     aavc is calculated can be modified using the aavcyrs argument. In
#'     addition to the boxplots, the function outputs the statistics describing
#'     the statistics from the boxplots.
#'
#' @param rundir the directory in which comparisons are being made. It is best
#'     to be a separate directory form any particular scenario.
#' @param hspm a list of the harvest strategy performance measures for actual
#'     catches as calculated by calccatchHSPM
#' @param glbc a list of the globals objects from each scenario
#' @param scenes the names of the different scenarios being compared
#' @param filen the filename for saving the plot, the default="", which sends
#'     the plot to the console.
#' @param aavcyrs the number of projection years to use for the aavc calculation
#'     the default = 10
#' @param maxvals the maximum y-axis value for the aavc, the sum5, and the sum10,
#'     default=c(0.25,750,1500)
#'
#' @seealso{
#'    \link{getprojyrs}, \link{getprojyraavc},  \link{getprojyrC}
#' }
#'
#' @return a list of the boxplot statistics for the three HSPMs x the number of
#'     scenarios compared
#' @export
#'
#' @examples
#' print("wait on data sets")
catchHSPM <- function(rundir,hspm,glbc,scenes,filen="",aavcyrs=10,
                      maxvals=c(0.25,750,1500)) {
  nscen <- length(scenes)
  aavc <- hspm$aavc
  sum5 <- hspm$sum5
  sum10 <- hspm$sum10
  label <- NULL
  for (i in 1:nscen) label <- c(label,paste0(scenes[i],"_",c("aavc","sum5","sum10")))
  boxresult <- vector(mode="list",length=(3 * nscen))
  names(boxresult) <- label
  nres <- length(boxresult)
  if (nchar(filen) > 0) {
    filen <- filenametopath(rundir,filen)
    caption <- paste0("Boxplots of aavc, sum5, and sum10 for each SAU.")
  }
  plotprep(width=9, height=9,newdev=FALSE,filename = filen,verbose=FALSE)
  parset(plots=c(3,2),margin=c(0.4,0.4,0.05,0.1),outmargin=c(1,0,0,0),byrow=FALSE)
  count <- 0
  for (i in 1:nscen) { # i = 1
    count <- count + 1
    boxresult[[count]] <- boxplot(aavc[[i]],ylim=c(0,maxvals[1]),yaxs="i",
                                  ylab="Average Annual Catch Variation")
    count <- count + 1
    boxresult[[count]] <- boxplot(sum5[[i]][,1:8],ylim=c(0,maxvals[2]),yaxs="i",
                                  ylab="sum first 5 years catch")
    count <- count + 1
    boxresult[[count]] <- boxplot(sum10[[i]][,1:8],ylim=c(0,maxvals[3]),yaxs="i",
                                  ylab="sum first 10 years catch")
    mtext(scenes[i],side=1,outer=FALSE,cex=1,line=1.75)
  }
  if (nchar(filen) > 0)
    addplot(filen=filen,rundir=rundir,category="catches",caption)
  for (i in 1:nres) { # i = 1
    namelab <- boxresult[[i]]$names
    stats <- boxresult[[i]]$stats
    colnames(stats) <- namelab
    rownames(stats) <- c("lowwhisk","25thQ","median","75thQ","highwhisk")
    boxresult[[i]] <- stats
  }
  return(invisible(boxresult))
} #end of catchHSPM

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

#' @title constructTasHSout remakes outputs from the hcrfun for all projections
#'
#' @description constructTasHSout uses the scenario's hcrfun to reconstrcut the
#'     scores and reference points for the Tasmanian HS.
#'
#' @param hcrfun the same as used in do_MSE = mcdahcr from TasHS package
#' @param ce the 3D array of cpue, which, for Tas, will include all years back
#'     to 1992, the start of defensible cpue data.
#' @param acatch the projected aspirational catches required by mcdahcr
#' @param glb the global object for the given scenario  out$glb
#' @param hsargs the hsargs used in teh scenario out$hsargs
#'
#' @return invisibly a list of multTAC, grad1, grad4,score1, score4, scoret,
#'     scoretot, and targrefpts. Each is a 3D array of years x sau x reps,
#'     except for the targrefpts which is the sau x 4 matrix of refpts for the
#'     final year in each replicate.
#' @export
#'
#' @examples
#' print("wait on data sets")
constructTasHSout <- function(hcrfun,ce,acatch,glb,hsargs) {
  startCE <- glb$indexCE
  finalyr <- glb$hyrs + glb$pyrs
  saunames <- glb$saunames
  nsau <- glb$nSAU
  reps <- glb$reps
  ce <- getprojyrs(ce,glb$hyrs,glb$pyrs,startyr=glb$indexCE)
  yrnames <- rownames(ce[,,1])
  nyrs <- length(yrnames)
  multTAC <- array(data=0,dim=c(nyrs,nsau,reps),
                   dimnames=list(yrnames,saunames,1:reps))
  grad1 <- grad4 <- score1 <- score4 <- scoret <- scoretot <- multTAC
  ptlabel <- c("low","trp","high","realtrp")
  targrefpts <- array(data=0,dim=c(nsau,length(ptlabel),reps),
                      dimnames=list(saunames,ptlabel,1:reps))
  for (iter in 1:reps) { # iter=1
    cpue <- ce[,,iter]
    hcrdata <- list(arrce=cpue,yearnames=yrnames,acatches=acatch[finalyr,,iter],
                    fis=NULL,nas=NULL)
    hcrout <- hcrfun(hcrdata,hsargs,saunames=saunames)
    multTAC[,,iter] <- hcrout$multTAC
    details <- hcrout$details
    grad1[,,iter] <- details$grad1
    grad4[,,iter] <- details$grad4
    score1[,,iter] <- details$score1
    score4[,,iter] <- details$score4
    scoret[,,iter] <- details$scoret
    scoretot[,,iter] <- details$scoretot
    targrefpts[,,iter] <- hcrout$refpts
  }
  return(invisible(list(multTAC=multTAC,grad1=grad1,grad4=grad4,score1=score1,
                        score4=score4,scoret=scoret,scoretot=scoretot,
                        targrefpts=targrefpts)))
} # end of constructrTasHSout

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


#' @title comparefinalscores compares the median final HS scores from scenarios
#'
#' @description comparefinalscores takes the 'scores' object from each 'out'
#'     object generated by do_MSE and plots, for each sau, the median final
#'     HS score for each scenario within each sau, using the scenario name as
#'     the legend entry. 'scores' is a list object made up of two list objects
#'     generated in do_MSE using both 'getcpueHS' and 'plotfinalscores'. Here we
#'     use only the outmed list and not the outHS list. The outmed list contains
#'     the medians of each of the replicates scores held in outHS, and this
#'     function plots the median final scores for each scenario against each
#'     other for each sau.
#'
#' @param rundir to location for all results of the comparisons
#' @param scores derived from the out$scores object saved from each scenario
#' @param scenes the scenario names
#' @param legloc the location of the legend, default = bottomright
#' @param filen the name of the png file if it is to be stored, default=""
#' @param category in what tab should the tables be stored. default="Tables"
#'
#' @return nothing but it does generate a plot.
#' @export
#'
#' @examples
#' print("wait on data")
comparefinalscores <- function(rundir,scores,scenes,legloc="bottomright",
                               filen="",category="Tables") {
  #   scores=scores;scenes=scenes;legloc="bottomleft"; filen="";
  nmed <- length(scores)
  meds <- vector(mode="list",length=nmed)
  for (i in 1:nmed) {
    tmp <- scores[[i]]$outmed
    meds[[i]] <-  sapply(tmp,"[[","medsc")
  }
  saunames <- colnames(meds[[1]])
  nsau <- length(saunames)
  makeplot <- function(meds,sau,scenes,legloc) {
    nmed <- length(meds)
    yrs <- as.numeric(rownames(meds[[1]]))
    plot(yrs,meds[[1]][,sau],type="l",lwd=2,col=1,ylim=c(0,10),ylab=sau,
         panel.first=grid())
    if (nmed > 1) for (i in 2:nmed) lines(yrs,meds[[i]][,sau],lwd=2,col=i)
    legend(legloc,legend=scenes,col=1:nmed,lwd=3,bty="n",cex=1.0)
  }
  if (nchar(filen) > 0) filen <- filenametopath(rundir,filen)
  plotprep(width=8, height=9,newdev=FALSE,filename=filen,verbose=FALSE)
  parset(plots=c(4,2),margin=c(0.25,0.4,0.05,0.1),byrow=FALSE,
         outmargin=c(0,1,0,0))
  for (i in 1:nsau) {
    makeplot(meds=meds,sau=saunames[i],scenes=scenes,legloc=legloc)
  }
  mtext("Median Final Scores",side=2,line=-0.1,outer=TRUE, cex=1.0)
  return(invisible(meds))
} # end of comparefinalscores


#' @title tabulatefinalHSscores saves and write the medina final scores to rundir
#'
#' @description tabulatefinalHSscores saves csv files for each scenario and
#'     includes a table of each of the HS final scores in the website output.
#'     The 'meds' object is obtained from the 'comparefinalscores' function.
#'
#' @param rundir the location for all results of the comparisons
#' @param meds a list of each scenarios median final HS scores from
#'     comparefinalscores
#' @param scenes the scenario names
#' @param category in what tab should the tables be stored. default="Tables"
#'
#' @return nothing but it does save two tables and generates the HTML for the
#'     output website
#' @export
#'
#' @examples
#' print("wait on data")
tabulatefinalHSscores <- function(rundir,meds,scenes,category="Tables") {
  for (i in 1:length(meds)) {
    filen <- paste0("finalscores_",scenes[i],".csv")
    caption <- paste0("Final HS scores for ",scenes[i],",")
    addtable(meds[[i]],filen,rundir=rundir,category=category,caption)
  }
} # end of tabulatefinalHSscores

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


#' @title plotfinalscores plots projected catches, cpue, PMs, and final score
#'
#' @description plotfinalscores takes the inidivual sau output from getcpueHS
#'     and generates a plot made up of six sub-plots of the replicates for each
#'     sau. In sequence these plots are the projected catches, the projected
#'     cpue, the gradient 1 scores, the gradient 4 scores, the target cpue, and
#'     the final total score. This enables the relationships between the scores
#'     to be examined. Only used in do_MSE
#'
#' @param outscores each sau's results from getcpueHS, which extracts the HS
#'     scores and other results from the projection dynamics
#' @param minprojC the minimum y-axis value for the projected catches
#' @param minprojCE the minimum y-axis value for the projected cpue
#' @param mintargCE the minimum y-axis value for the projected target cpue
#' @param filen a filename for use if saving the output to said file,
#'     default = '', which implies the plot goes to the console.
#'
#' @seealso{
#'    \link{do_MSE}
#' }
#'
#' @return invisibly a list of the median catches, cpue, targetce, finalscore,
#'     grad1 score, grad4 score, and target cpue score.
#' @export
#'
#' @examples
#' print("wait on data sets")
plotfinalscores <- function(outscores,minprojC=0,minprojCE=60,mintargCE=110,
                            filen="") {
  cpue <- outscores$cpue
  yrs <- as.numeric(rownames(cpue))
  reps <- ncol(cpue)
  cetarg <- outscores$cetarg
  catch <- outscores$catch
  score <- outscores$finalsc
  tyrs <- as.numeric(rownames(cetarg))
  pickyr <- match(tyrs,yrs)
  ymax <- getmax(catch)    # first plot catches
  medcpue <- apply(cpue[pickyr,],1,median)
  medtarg <- apply(cetarg,1,median)
  medcatch <- apply(catch,1,median)
  medsc <- apply(score,1,median)
  plotprep(width=8, height=9,newdev=FALSE,filename=filen,verbose=FALSE)
  parset(plots=c(7,1),margin=c(0.2,0.4,0.025,0.1))
  plot(tyrs,catch[,1],type="l",lwd=1,col="grey",ylim=c(minprojC,ymax),
       ylab="Projected Catch",xlab="",panel.first=grid())
  for (i in 1:reps) lines(tyrs,catch[,i],lwd=1,col="grey")
  lines(tyrs,medcatch,lwd=2,col=2)
  ymax <- getmax(cpue[pickyr,])  # now plot cpue
  plot(tyrs,cpue[pickyr,1],type="l",lwd=1,col="grey",ylim=c(minprojCE,ymax),
       ylab="Projected CPUE",xlab="",panel.first=grid())
  for (i in 1:reps) lines(tyrs,cpue[pickyr,i],lwd=1,col="grey")
  lines(tyrs,medcpue,lwd=2,col=2)
  lines(tyrs,medtarg,lwd=2,col=4)
  g1s <- outscores$g1s  # now the grad1 score
  medg1s <- apply(g1s,1,median)
  plot(tyrs,g1s[,1],type="l",lwd=1,col="grey",ylim=c(0,10),
       panel.first=grid(),ylab="Grad1 Score",xlab="")
  for (i in 1:reps) lines(tyrs,g1s[,i],lwd=1,col="grey")
  lines(tyrs,medg1s,lwd=2,col=2)
  g4s <- outscores$g4s  # now the grad4 score
  medg4s <- apply(g4s,1,median)
  plot(tyrs,g4s[,1],type="l",lwd=1,col="grey",ylim=c(0,10),
       panel.first=grid(),ylab="Grad4 Score",xlab="")
  for (i in 1:reps) lines(tyrs,g4s[,i],lwd=1,col="grey")
  lines(tyrs,medg4s,lwd=2,col=2)
  ts <- outscores$ts  # now the target cpue score
  medts <- apply(ts,1,median)
  plot(tyrs,ts[,1],type="l",lwd=1,col="grey",ylim=c(0,10),
       panel.first=grid(),ylab="Target CE Score",xlab="")
  for (i in 1:reps) lines(tyrs,ts[,i],lwd=1,col="grey")
  lines(tyrs,medts,lwd=2,col=2)
  ymax <- getmax(cetarg,mult=1.01)   # now the target cpue
  plot(tyrs,cetarg[,1],type="l",lwd=1,col="grey",ylim=c(mintargCE,ymax),
       panel.first=grid(),ylab="Target CE",xlab="")
  for (i in 1:reps) lines(tyrs,cetarg[,i],lwd=1,col="grey")
  abline(h=150,lwd=2,col=2,lty=2)
  lines(tyrs,medtarg,lwd=2,col=2)
  ylabel <- "Final TasHS Score"  # now the final scores
  plot(tyrs,score[,1],type="l",lwd=1,col="grey",ylim=c(0,10),yaxs="i",
       panel.first=grid(),ylab=ylabel,xlab="")
  for (i in 1:reps) lines(tyrs,score[,i],lwd=1,col="grey")
  lines(tyrs,medsc,lwd=2,col=2)
  return(invisible(list(medcatch=medcatch,medcpue=medcpue,medtarg=medtarg,
                        medsc=medsc,medg1s=medg1s,medg4s=medg4s,medts=medts)))
} # end of plotfinalscores

#' @title tabulateproductivity produces tables of the MSY and related statistics
#'
#' @description tabulateproductivity writes out csv files for each scenario
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
tabulateproductivity <- function(rundir,prods,scenes) {
  nscen <- length(scenes)
  for (i in 1:nscen) {
    filen <- paste0("productivity_",scenes[i],".csv")
    caption <- paste0("SAU productivity statistics for ",scenes[i],".")
    addtable(prods[[i]],filen,rundir=rundir,category="productivity",caption)
  }
} # end of tabulateproductivity


