#' @title boxbysau generates boxplots of sau x aavc or sum5 or sum10
#'
#' @description boxbysau is used to generate boxplots of the aavc, or sum5,
#'     of sum10 for each sau, in each of the scenarios. The number of years over
#'     which the aavc is calculated can be modified using the aavcyrs argument.
#'     If the replicates differ between scenarios then these boxplots will be
#'     omitted.
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
  allreps <- unlist(lapply(glbc,"[[","reps"))
  if (min(allreps) != max(allreps)) {
    label <- paste0("Cannot compare scenarios with different replicate numbers",
                    " catch box plots omitted \n")
    warning(cat(label))
  } else {
    reps <- max(allreps)
    boxdat <- matrix(0,nrow=reps,ncol=nscen,dimnames=list(1:reps,scenes))
    if (nchar(filen) > 0) {
      filen <- filenametopath(rundir,filen)
      caption <- paste0("Boxplots of ",compvar," for each SAU. Redline = median of first scenario.")
    }
    plotprep(width=8, height=8,newdev=FALSE,filename = filen,verbose=FALSE)
    parset(plots=pickbound(nsau),margin=c(0.3,0.4,0.05,0.1),
           outmargin=c(1,0,0.25,0),byrow=FALSE)
    for (sau in 1:nsau) { # sau=1; i = 1
      for (i in 1:nscen) {
        tmp <- boxvar[[i]][,sau]
        boxdat[1:length(tmp),i] <- tmp
      }
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
      addplot(filen=filen,rundir=rundir,category="catchBoxPlots",caption)
  }
} # end of boxbysau


#' @title calccatchHSPM calculates all the HSPMs for the actual catch data
#'
#' @description calccatchHSPM calculates all the Harvest Strategy Performance
#'     Measures for the actual catch data. The number of years over which the
#'     aavc is calculated can be modified using the aavcyrs argument. The output
#'     data can be analysed separately or plotted (eg as boxplots)
#'
#' @param catch a list of the actual catches from the out$sauout$zonePsau$catch
#'     array from each scenario's out object
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
#'     same sau and scenes. Plus the medians across all replicates for each
#'     sau and zone for each of the three PMs.
#' @export
#'
#' @examples
#' print("wait on data sets")
calccatchHSPM <- function(catch,glbc,scenes,aavcyrs=10) {
  # catch=catch;glbc=glbc;scenes=scenes;aavcyrs=10
  nscen <- length(scenes)
  catches <- makelist(scenes)
  aavc <- makelist(scenes)
  sum5 <- makelist(scenes)
  sum10 <- makelist(scenes)
  saunames <- glbc[[1]]$saunames
  for (i in 1:nscen) {   # i = 1
    catches[[i]] <- projectiononly(catch[[i]],glbc[[i]])
    aavc[[i]] <- getprojyraavc((catches[[i]][1:aavcyrs,,]),glbc[[i]])
    sum5[[i]] <- getprojyrC(catch[[i]],glbc[[i]],period=5)
    sum10[[i]] <- getprojyrC(catch[[i]],glbc[[i]],period=10)
  }
  msePM <- c("aavc","sum5","sum10")
  medians <- makelist(msePM)
  medians[[1]] <- sapply(aavc,function(x) apply(x,2,median))
  medians[[2]] <- sapply(sum5,function(x) apply(x,2,median))
  medians[[3]] <- sapply(sum10,function(x) apply(x,2,median))
  return(invisible(list(aavc=aavc,sum5=sum5,sum10=sum10,medians=medians)))
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
#' @param wid the width of hteplot, default = 12, room for 3 scenarios
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
#' # rundir=rundir;hspm=outcatchHSPM;glbc=glbc;scenes=scenes
#' # filen="compare_catches_boxplots.png";aavcyrs=10
catchHSPM <- function(rundir,hspm,glbc,scenes,filen="",aavcyrs=10,
                      maxvals=c(0.25,750,1500),wid=12) {
  nscen <- length(scenes)
  aavc <- hspm$aavc
  sum5 <- hspm$sum5
  sum10 <- hspm$sum10
  nsau <- glbc[[1]]$nSAU
  label <- NULL
  for (i in 1:nscen) label <- c(label,paste0(scenes[i],"_",c("aavc","sum5","sum10")))
  boxresult <- vector(mode="list",length=(3 * nscen))
  names(boxresult) <- label
  nres <- length(boxresult)
  if (nchar(filen) > 0) {
    filen <- filenametopath(rundir,filen)
    caption <- paste0("Boxplots of aavc, sum5, and sum10 for each SAU and scenario.")
  }
  plotprep(width=wid, height=9,newdev=FALSE,filename = filen,verbose=FALSE)
  parset(plots=c(3,nscen),margin=c(0.3,0.4,0.05,0.1),outmargin=c(0,0,1.5,0),
         byrow=FALSE)
  count <- 0
  maxval <- numeric(nscen)
  for (i in 1:nscen) { # i = 1
    count <- count + 1
    for (j in 1:nscen) maxval[j] <- getmax(aavc[[j]])
    boxresult[[count]] <- boxplot(aavc[[i]],ylim=c(0,max(maxval)),yaxs="i",
                                  ylab="Average Annual Catch Variation")
    mtext(scenes[i],side=3,outer=FALSE,cex=1.2,line=0.5)
    count <- count + 1
    for (j in 1:nscen) maxval[j] <- getmax(sum5[[j]][,1:nsau])
    boxresult[[count]] <- boxplot(sum5[[i]][,1:nsau],ylim=c(0,max(maxval)),yaxs="i",
                                  ylab="sum first 5 years catch")
    count <- count + 1
    for (j in 1:nscen) maxval[j] <- getmax(sum10[[j]][,1:nsau])
    boxresult[[count]] <- boxplot(sum10[[i]][,1:nsau],ylim=c(0,max(maxval)),yaxs="i",
                                  ylab="sum first 10 years catch")
  }
  if (nchar(filen) > 0)
    addplot(filen=filen,rundir=rundir,category="catches",caption)
  for (i in 1:nres) { # i = 1
    namelab <- boxresult[[i]]$names
    bstat <- boxresult[[i]]$stats
    colnames(bstat) <- namelab
    rownames(bstat) <- c("lowwhisk","25thQ","median","75thQ","highwhisk")
    boxresult[[i]] <- bstat
  }
  return(invisible(boxresult))
} #end of catchHSPM

#' @title catchvsMSY by sau divides the projected catches by the sau MSY value
#'
#' @description catchvsMSY is used when calculating the ratio of the projected
#'     catches to the MSY. This uses the estimate of MSY for each sau that is
#'     given in out$sauprod object from each scenario, each of which is placed
#'     into the 'prods' list.
#'
#' @param catch the projected catches array of all years x nsau x replicates
#' @param glbc the globals object for each scenario
#' @param prods a list of the production statistics matrix for each scenario,
#'     this is the out$sauprod matrix stored for each scenario.
#' @param scenes  the names of the different scenarios being compared in the
#'     same order in which the out objects are decomposed into their components
#'
#' @return a list of 3D arrays of the ratio of catch/MSY for each sau in each
#'     scenario
#' @export
#'
#' @examples
#' print("wait on data sets")
#' # catch=catch;glbc=resout$glbc;prods=prods;scenes=scenes
catchvsMSY <- function(catch,glbc,prods,scenes) {
  nscen <- length(scenes)
  cdivmsy <- makelist(scenes)
  yrtomsy <-makelist(scenes)
  glb <- glbc[[1]]
  yrs <- glb$pyrnames
  nyrs <- glb$pyrs
  hyrs <- glb$hyrs
  reps <- glb$reps
  for (i in 1:nscen) {   # i = 1
    glb <- glbc[[i]]
    nsau <- glb$nSAU
    tomsy <- matrix(0,nrow=reps,ncol=nsau,dimnames=list(1:reps,1:nsau))
    cvsmsy <- getprojyrs(catch[[i]],glb$hyrs,glb$pyrs,startyr=glb$hyrs)
    tmp <- cvsmsy
    prodsc <- prods[[i]]
    msy <- prodsc[,"MSY"]
    for (sau in 1:nsau) {
      cvsmsy[,sau,] <- tmp[,sau,]/msy[sau]
      tomsy[,sau] <- apply(tmp[,sau,],2,whenfirst,value=msy[sau])
    }
    cdivmsy[[i]] <- cvsmsy
    yrtomsy[[i]] <- tomsy
  }
  return(invisible(list(cdivmsy=cdivmsy,yrtomsy=yrtomsy)))
} # end of catchvsMSY

#' @title CBmsyphaseplot makes a zone scale phase plot of Cdivmsy vs BdivBmsy
#'
#' @description CBmsyphaseplot is used within do_comparison to compare harvest
#'     strategy scenarios at the zone level by using phase plots of predicted
#'     catches divided by MSY against predicted mature biomass divided by Bmsy.
#'     Some HS can recommend taking catches >= MSY even though the biomass
#'     available is not sufficient for such catches to be sustainable or allow
#'     for rebuilding. The percent of replicates within the dangerous region
#'     where catches >= MSY yet biomass is below Bmsy is given topleft. Top
#'     right is the goodplace, with the  numbers being the minimum, median,
#'     and maximum number of years it takes for both catches >= MSY AND
#'     mature biomass >= Bmsy.
#'
#' @param rundir the directory in which all results are held for a scenario or
#'     comparison of scenarios
#' @param zone the list of zone outputs for each scenario
#' @param prods a list of the production statistics matrix for each scenario,
#'     this is the out$sauprod matrix stored for each scenario.
#' @param glbc the globals object for each scenario (needed in case the hyrs
#'     differ between scenarios)
#' @param category the tab name within the internal website into which to
#'     place he plot, default = 'a' to ensure it has something
#' @param console should each plot go to the console or be saved to rundir.
#'     default = TRUE ie go to console
#'
#' @seealso  \link{do_comparison}
#'
#' @returns an invisible list of lists - Cvsmsy, BvsBmsy, and goodplace
#' @export
#'
#' @examples
#' print("wait on example data sets")
#' # CBmsyphaseplot(rundir=rundir,zone=zone,prods=prods,glbc=glbc,
#' #                category="zone",console=FALSE)
CBmsyphaseplot <- function(rundir,zone,prods,glbc,category="a",console=TRUE) {
  scenes <- names(zone)
  nscen <- length(scenes)
  glb <- glbc[[1]]
  hyrs <- glb$hyrs
  pyrs <- glb$pyrs
  reps <- glb$reps
  BvsspawnB <- makelist(scenes)
  Cvsmsy <- makelist(scenes)
  catch <- makelist(scenes)
  spawnB <- makelist(scenes)
  zonemsy <- numeric(3); names(zonemsy) <- scenes
  zoneBmsy <- zonemsy
  for (scen in 1:nscen) {
    catch[[scen]] <- zone[[scen]]$catch[(hyrs+1):(hyrs+pyrs),]
    spawnB[[scen]] <- zone[[scen]]$matureB[(hyrs+1):(hyrs+pyrs),]
    zonemsy[scen] <- sum(prods[[scen]][,"MSY"])
    zoneBmsy[scen] <- sum(prods[[scen]][,"Bmsy"])
  }
  minC <- 1
  maxC <- 0
  minB <- 1
  maxB <- 0
  for (scen in 1:nscen) { # define the ratios and their ranges
    msy <- zonemsy[scen]
    Bmsy <- zoneBmsy[scen]
    catches <- catch[[scen]]
    Cvsmsy[[scen]] <- catches/msy
    getC <- getmin(Cvsmsy[[scen]])
    if (getC < minC) minC <- getC
    getC <- getmax(Cvsmsy[[scen]])
    if (getC > maxC) maxC <- getC
    matureB <- spawnB[[scen]]
    BvsspawnB[[scen]] <- matureB/Bmsy
    getB <- getmin(BvsspawnB[[scen]])
    if (getB < minB) minB <- getB
    getB <- getmax(BvsspawnB[[scen]])
    if (getB > maxB) maxB <- getB
  }
  filen=""
  if (!console) {
    filename <- "C-vs-MSY_B-vs-Bmsy_by_zone_phaseplot.png"
    filen <- pathtopath(rundir,filename)
    caption <- paste0("Phase plot of Cmsy and Bmsy ratio for the Zone scale",
                      "topleft number is percent replicates catching >= MSY ",
                      "when biomass is too low. Top right is minimum, median, ",
                      "and maximum number of years for both to be >= MSY and ",
                      ">= Bmsy, with the last number being the percent of ",
                      "achieving >= Bmsy at the same time as catch >= MSY.")
  }
  goodplace <- makelist(scenes)
  plotprep(width=9,height=10,newdev=FALSE,filename=filen,verbose=FALSE)
  parset(plots=c(nscen,1),margin=c(0.5,0.5,0.1,0.1),cex=1.0)
  for (scen in 1:nscen) {  # scen = 1
    ylabel <- paste0("C vs MSY ",scenes[scen])
    BM <- BvsspawnB[[scen]]
    CM <- Cvsmsy[[scen]]
    pickG <- which((BM >= 1.0) & (CM >= 1.0))
    pickB <- which((BM < 1.0) & (CM >= 1.0))
    plot(BM,CM,type="p",xlab="B vs Bmsy",
         ylab=ylabel,xlim=c(minB,maxB),ylim=c(minC,maxC),panel.first=grid())
    points(BM[pickG],CM[pickG],pch=16,col=3)
    points(BM[pickB],CM[pickB],pch=16,col=2)
    abline(v = 1, lwd=2,col=4)
    abline(h = 1, lwd=2, col=4)
    both <-  pickG %% pyrs
    projyr <- which(both == 0)
    if (length(projyr) > 0) both[projyr] <- pyrs
    goodplace[[scen]] <- both
    years <- ""
    if (length(both) > 0)
       years <- paste(min(both),median(both),max(both),collapse=" ")
    toplab <- paste0(years,"   ",round(100*length(pickG)/(pyrs * reps),1),"    ")
    mtext(toplab,side=3,outer=FALSE,line=-2,adj=1)
    leftlab <- paste0("   ",round(100*length(pickB)/(pyrs * reps),1))
    mtext(leftlab,side=3,outer=FALSE,line=-2,adj=0)
  }
  if (!console) {
    addplot(filen=filen,rundir=rundir,category=category,caption=caption)
  }
  return(invisible(list(Cvsmsy=Cvsmsy,BvsBmsy=BvsspawnB,goodplace=goodplace)))
} # end of CBmsyphaseplot


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
#' @description constructTasHSout uses the scenario's hcrfun to reconstruct the
#'     scores and reference points for the Tasmanian HS.
#'
#' @param hcrfun the same as used in do_MSE = mcdahcr from TasHS package
#' @param ce the 3D array of cpue, which, for Tas, will include all years back
#'     to 1992, the start of defensible cpue data.
#' @param acatch the projected aspirational catches required by mcdahcr
#' @param glb the global object for the given scenario  out$glb
#' @param hsargs the hsargs used in the scenario out$hsargs
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
#' # rundir=rundir;dyn=dyn;glbc=glbc;scenes=scenes
#' }
comparedynamics <- function(rundir,dyn,glbc,scenes) {
  quantscen <- makelist(c("cpue","catch","matureB","harvestR"))
  invar <- "cpue"
  compscen <- comparevar(dyn,glbc,scenes,var=invar)
  quantscen[[1]] <- compscen$quantscen
  filename <- filenametopath(rundir,"compare_cpue_scenes.png")
  plotscene(compscen$quantscen,glbc=glbc,var=invar,ymin=0,filen=filename,
            legloc="bottomleft",legplot=1)
  caption <- paste0("The projected CPUE dynamics for each SAU.")
  addplot(filen=filename,rundir=rundir,category="Dynamics",caption)
  invar <- "catch"
  compscen <- comparevar(dyn,glbc,scenes,var=invar)
  quantscen[[2]] <- compscen$quantscen
  filename <- filenametopath(rundir,"compare_catch_scenes.png")
  plotscene(compscen$quantscen,glbc=glbc,var=invar,ymin=0,filen=filename,
            legloc="topleft",legplot=1)
  caption <- paste0("The projected actual catch dynamics for each SAU.")
  addplot(filen=filename,rundir=rundir,category="Dynamics",caption)

  invar <- "matureB"
  compscen <- comparevar(dyn,glbc,scenes,var=invar)
  quantscen[[3]] <- compscen$quantscen
  filename <- filenametopath(rundir,"compare_matureB_scenes.png")
  plotscene(compscen$quantscen,glbc=glbc,var=invar,ymin=0,filen=filename,
            legloc="topleft",legplot=1)
  caption <- paste0("The projected mature biomass dynamics for each SAU.")
  addplot(filen=filename,rundir=rundir,category="Dynamics",caption)

  invar <- "harvestR"
  compscen <- comparevar(dyn,glbc,scenes,var=invar)
  quantscen[[4]] <- compscen$quantscen
  filename <- filenametopath(rundir,"compare_harvestR_scenes.png")
  plotscene(compscen$quantscen,glbc=glbc,var=invar,ymin=0,filen=filename,
            legloc="topleft",legplot=1)
  caption <- paste0("The projected harvest rate dynamics for each SAU.")
  addplot(filen=filename,rundir=rundir,category="Dynamics",caption)
  return(quantscen)
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
  # rundir=rundir;scores=scores;scenes=scenes;legloc="bottomright"
  # filen="";category="Scores"
  nmed <- length(scores)
  meds <- vector(mode="list",length=nmed)
  for (i in 1:nmed) { # i=2
    tmp <- scores[[i]]$finalsc
    indim <- dim(tmp)
    dimlabel <- dimnames(tmp)
    nyrs <- indim[1]
    yrnames <- dimlabel[[1]]
    nsau <- indim[2]
    saunames <- dimlabel[[2]]
    medfsc <- matrix(0,nrow=nyrs,ncol=nsau,dimnames=list(yrnames,saunames))
    if (is.null(tmp)) { #do what when tmp == null
      scenes[i] <- paste0(scenes[i]," No Final Score")
    } else {
      for (sau in 1:nsau) {  # sau = 1
        medfsc[,sau] <- apply(tmp[,sau,],1,median)
      }
    }
    meds[[i]] <-  medfsc
  }
  makeplot <- function(meds,sau,scenes,legloc) {
    nmed <- length(meds)
    rows <- sapply(meds,nrow)
    maxrow <- max(rows)
    if (all(rows == maxrow)) {
      yrs <- as.numeric(rownames(meds[[1]]))
    } else {
      pickR <- which(rows == maxrow)
      yrs <- as.numeric(rownames(meds[[pickR]]))
    }
    plot(yrs,meds[[1]][,sau],type="l",lwd=2,col=1,ylim=c(0,10),ylab=sau,
         panel.first=grid())
    if (nmed > 1) for (i in 2:nmed) lines(yrs,meds[[i]][,sau],lwd=2,col=i)
    legend(legloc,legend=scenes,col=1:nmed,lwd=3,bty="n",cex=1.0)
  }
  if (nchar(filen) > 0) filen <- filenametopath(rundir,filen)
  plotprep(width=8, height=9,newdev=FALSE,filename=filen,verbose=FALSE)
  parset(plots=pickbound(nsau),margin=c(0.25,0.4,0.05,0.1),byrow=FALSE,
         outmargin=c(0,1,0,0))
  for (i in 1:nsau) { # i = 1
    makeplot(meds=meds,sau=saunames[i],scenes=scenes,legloc=legloc)
  }
  mtext("Median Final Scores",side=2,line=-0.1,outer=TRUE, cex=1.0)
  return(invisible(meds))
} # end of comparefinalscores

#' @title cpueboxbysau generates boxplots of sau x years to max median cpue
#'
#' @description cpueboxbysau is used to generate boxplots of the years taken to
#'     reach the maximum median cpue for each sau x senario. It also provides
#'     the summary statistics from the boxplots for each sau. If the replicates
#'     differ between scenarios then these boxplots will be omitted.
#'
#' @param rundir the directory in which comparisons are being made. It is best
#'     to be a separate directory form any particular scenario.
#' @param cpue a list of the cpue arrays from the out$sauout$zonePsau$cpue
#'     arrays from each scenario's out object
#' @param glbc a list of the globals objects from each scenario
#' @param scenes the names of the different scenarios being compared
#' @param filen the filename for saving the plot, the default="", which sends
#'     the plot to the console.
#' @param startyr what year to start the search for the maximum median cpue,
#'     default = 0, this means that only the projection years are used.
#' @param maxval the maximum y-axis value for the years to max median cpue,
#'     default=0, which means the maximum for each boxplot is taken from the
#'     data for that plot.
#'
#' @seealso{
#'    \link{getprojyrs}, \link{getprojyraavc},  \link{catchHSPM}
#' }
#'
#' @return a matrix of the boxplot statistics for the years to max median cpue
#'     for each scenario and sau considered.
#' @export
#'
#' @examples
#' print("wait on data sets")
cpueboxbysau <- function(rundir,cpue,glbc,scenes,filen="",startyr=0,maxval=0) {
  #  rundir=rundir;cpue=cpue;glbc=glbc;scenes=scenes;filen="sau_maxceyr_box.png";startyr=0; maxval=0
  nscen <- length(scenes)
  nsau <- glbc[[1]]$nSAU
  saunames <- glbc[[1]]$saunames
  columns <- c("sau","stat",scenes)
  numcol=length(columns)
  cpuebox <- as.data.frame(matrix(0,nrow=(nsau*5),ncol=numcol,
                                  dimnames=list(1:(nsau*5),columns)))
  label <- NULL
  for (sau in 1:nsau) label <- c(label,rep(saunames[sau],5))
  cpuebox[,1] <- label
  cpuebox[,2] <- rep(c("lowwhisk","25thQ","median","75thQ","highwhisk"),nsau)
  allreps <- unlist(lapply(glbc,"[[","reps"))
  if (min(allreps) != max(allreps)) {
    label <- paste0("Cannot compare scenarios with different replicate numbers",
                    " so boxplots of maximum cpue year by sau omitted \n")
    warning(cat(label))
  } else {
    reps <- max(allreps)
    boxdat <- matrix(0,nrow=reps,ncol=nscen,dimnames=list(1:reps,scenes))
    if (startyr == 0) {
      beginyr <- glbc[[1]]$hyrs
    } else {
      beginyr <- startyr
    }
    ce <- makelist(scenes)
    maxceyr <- makelist(scenes)
    for (i in 1:nscen) { # i =1
      glb <- glbc[[i]]
      ce[[i]] <- getprojyrs(cpue[[i]],glb$hyrs,glb$pyrs,startyr=beginyr)
      maxceyr[[i]] <- getyr2maxce(ce[[i]],glb)
    }
    if (nchar(filen) > 0) {
      filen2 <- filenametopath(rundir,filen)
      caption <- paste0("Boxplots of year of maximum cpue by sau.",
                        "Redline is median max-cpue for first scenario.")
    }
    start=1; finish=5
    plotprep(width=8, height=8,newdev=FALSE,filename = filen2,verbose=FALSE)
    parset(plots=pickbound(nsau),margin=c(0.3,0.4,0.05,0.1),outmargin=c(1,0,0.25,0),
           byrow=FALSE)
    for (sau in 1:nsau) { # sau=1; i = 1

      for (i in 1:nscen) {
        tmp <- maxceyr[[i]][,sau]
        boxdat[1:length(tmp),i] <- tmp
      }
      if (maxval > 0) {
        maxy <- maxval
      } else {
        maxy <- getmax(boxdat)
      }
      tmp <- boxplot(boxdat,ylim=c(0,maxy),yaxs="i",ylab=saunames[sau])
      abline(h=tmp$stats[3,1],lwd=2,lty=2,col=2)
      cpuebox[start:finish,3:numcol] <- tmp$stats
      start <- finish + 1
      finish <- finish + 5
    }
    label <- "Years to Maximum Median CPUE by sau"
    mtext(label,side=1,outer=TRUE,line=-0.3,cex=1.0)
    if (nchar(filen) > 0) {
      addplot(filen=filen2,rundir=rundir,category="cpueBoxPlots",caption)
      filen3 <- unlist(strsplit(filen,".",fixed=TRUE))[1]
      filen3 <- paste0(filen3,".csv")
      filen <- filenametopath(rundir,filen3)
      write.csv(cpuebox,file=filen)
      caption <- "Years to maximum median cpue by sau from boxplots."
      addtable(cpuebox,filen3,rundir=rundir,category="cpueBoxPlots",caption)
    }
    return(invisible(cpuebox))
  }
} # end of cpueboxbysau

#' @title cpueHSPM generates boxplots of cpue HSPMs and output statistics
#'
#' @description cpueHSPM is used to generate boxplots of the year of the maximum
#'     median cpue for all projections. The number of years over which the
#'     aavc is calculated can be modified using the startyr argument. In
#'     addition to the boxplots, the function outputs the statistics describing
#'     the statistics from the boxplots.
#'
#' @param rundir the directory in which comparisons are being made. It is best
#'     to be a separate directory form any particular scenario.
#' @param cpue a list of the cpue arrays from the out$sauout$zonePsau$cpue
#'     arrays from each scenario's out object
#' @param glbc a list of the globals objects from each scenario
#' @param scenes the names of the different scenarios being compared
#' @param filen the filename for saving the plot, the default="", which sends
#'     the plot to the console.
#' @param startyr what year to start the search for the maximum median cpue,
#'     default = 0, this means that only the projection years are used.
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
cpueHSPM <- function(rundir,cpue,glbc,scenes,filen="",startyr=0) {
  #  rundir=rundir;cpue=cpue;glbc=glbc;scenes=scenes;filen="sau_maxce_boxplots.png";startyr=0
  nscen <- length(scenes)
  boxresult <- makelist(scenes)
  if (nchar(filen) > 0) {
    filen <- filenametopath(rundir,filen)
    caption <- paste0("Boxplots of year of maximum cpue.")
  }
  plotprep(width=12, height=9,newdev=FALSE,filename = filen,verbose=FALSE)
  parset(plots=pickbound(nscen),margin=c(0.4,0.4,0.05,0.1),outmargin=c(1,0,0,0),byrow=FALSE)
  for (i in 1:nscen) { # i = 1
    glb <- glbc[[i]]
    if (startyr == 0) {
      beginyr <- glb$hyrs
    } else {
      beginyr <- startyr
    }
    ce <- getprojyrs(cpue[[i]],glb$hyrs,glb$pyrs,startyr=beginyr)
    maxceyr <- getyr2maxce(ce,glb)
    qs <- boxplot(maxceyr,ylim=c(0,max(maxceyr)+2),yaxs="i",
                  ylab="Years to maximum CPUE",xlab="")
    mtext(scenes[i],side=1,outer=FALSE,cex=1,line=1.75)
    qsstats <- qs$stats
    colnames(qsstats) <- glb$saunames
    rownames(qsstats) <- c("lower_whisker","lower_hinge","median","upper_hinge",
                           "upper_whisker")
    boxresult[[i]] <- qsstats
  }
  if (nchar(filen) > 0)
    addplot(filen=filen,rundir=rundir,category="cpueBoxPlots",caption)
  return(invisible(boxresult))
} #end of cpueHSPM

#' @title do_comp_outputs extracts components of the output from do_comparison
#'
#' @description do_comp_outputs is a helper function that extracts a number of
#'     the more important components from each scenario into lists of their
#'     own. It does this for 'dyn', which contains the replicate outputs from
#'     the dynamics summarized at the sau scale (this includes Nt and catchN).
#'     There is also 'zone', which does the same thing but at a zone scale.
#'     There is also scenes, which contains the names of each scenario, 'prods',
#'     which contains the productivity characteristics of each scenario, and
#'     'scores', which contains the HS scores for each internal performance
#'     measure.
#'
#' @param result the output from the do_comparison function
#' @param projonly should only the projection years be returned? default=TRUE.
#'
#' @seealso{
#'    \link{do_comparison}
#' }
#'
#' @return a list of one vector of names and 7 potentially large lists. Most
#'     importantly, popout has the dynamics at a population scale, dyn has them
#'     at an sau scale, and zone has them at the zone scale.
#' @export
#'
#' @examples
#' print("wait on datasets")
#' # result=result; projout=TRUE
do_comp_outputs <- function(result,projonly=TRUE) {
  nscen=length(result$ans)
  scenes <- result$scenes
  dyn <- makelist(scenes)
  glbc <- makelist(scenes)
  prods <- makelist(scenes)
  scores <- makelist(scenes)
  zone <- makelist(scenes)
  popout <- makelist(scenes)
  for (i in 1:nscen) { # i = 1
    out <- result$ans[[i]]
    glb <- out$glb
    glbc[[i]] <- glb
    prods[[i]] <- out$sauprod
    scores[[i]] <- out$scores
    if (projonly) {
      nobj <- length(out$zoneDP)
      objnames <- names(out$zoneDP)
      for (obj in 1:nobj)
        popout[[i]][[objnames[obj]]] <- projectiononly(out$zoneDP[[objnames[obj]]],glbc[[i]])
      nobj <- length(out$sauout)
      objnames <- names(out$sauout)
      for (obj in 1:nobj)
        dyn[[i]][[objnames[obj]]] <- projectiononly(out$sauout[[objnames[obj]]],glbc[[i]])
      nobj <- length(out$outzone)
      objnames <- names(out$outzone)
      for (obj in 1:nobj)
        zone[[i]][[objnames[obj]]] <- projectiononly(out$outzone[[objnames[obj]]],glbc[[i]])
    } else {
      popout[[i]] <- out$zoneDP
      dyn[[i]] <- out$sauout
      zone[[i]] <- out$outzone
    }
  }
  return(invisible(list(scenes=scenes,popout=popout,dyn=dyn,zone=zone,
                        glbc=glbc,prods=prods,scores=scores)))
} # end of do_comp_outputs

#' @title doquantplot generates a plot of multiple quantile values
#'
#' @description doquantplot is used when trying to present a summary of multiple
#'     scenarios applied across a sequence of years (for example in the
#'     projection years of an MSE). It assumes that for each scenario quantiles
#'     have been calculated for each year so that the input array if 5 x nyrs x
#'     nscenes, where nyrs is the number of years and nscenes is the number of
#'     scenarios. Five quantiles are assumed with values
#'     probs=c(0.025, 0.05, 0.5, 0.95 0.975), and the option exists of plotting
#'     either the 90th (0.05 and 0.95) or the 95th (0.025, 0.975). If the input
#'     quantiles have different values then the 'q90' argument,if TRUE, will
#'     reference the 2nd and 4th row and the 1st and 5th if FALSE. If polys is
#'     TRUE then filled transparent polygons will be plotted, with colours in
#'     the same sequence as the scenarios, alternatively, if polys=FALSE, then
#'     lines will be plotted instead. This function is currently only used by
#'     plotzonedyn
#'
#' @param varq an array of quantiles with dimensions the 5 quantiles of the
#'     variable in question (2.5, 5, 50, 95 97.5), by years, the number of
#'     years for which quantiles are available
#' @param varname The names of the variable being summarized, used as the
#'     y-label on the plot
#' @param yrnames the numeric values of the years to be plotted. Used as the
#'     x-axis as well as the x-axis labels
#' @param scenes the names given to the different scenarios
#' @param q90 should the 90th or 95th quantiles be plotted, TRUE = 90th
#' @param polys should transparent polygons be plotted = TRUE, or lines = FALSE
#' @param intens if polys=TRUE then intens signifies the intensity of colour on
#'     a scale of 0 - 255. 127 is about 50 percent dense.
#' @param addleg add a legend? default=FALSE
#' @param locate where to place legend if there is one, default='bottomright'
#' @param hlin a vector of values for each scenario or NULL
#' @param scencol what sequence of colours should each scenario have. The
#'     default=NULL which implies the colour will reflect the sequence in
#'     which they are read in.
#'
#' @seealso {
#'    \link{plotzonedyn}, \link{RGB}
#' }
#'
#' @return nothing but it does add polygons or lines to a plot
#' @export
#'
#' @examples
#' print("wait on datasets")
doquantplot <- function(varq,varname,yrnames,scenes,q90,polys,intens=127,
                        addleg=FALSE,locate="bottomright",hlin=NULL,
                        scencol=NULL) {
  maxy <- getmax(varq)
  nscen <- length(scenes)
  if (is.null(scencol)) {
    scencol <- 1:nscen
  }
  plot(yrnames,varq[3,,1],type="l",lwd=2,col=0,xlab="",ylab=varname,
       panel.first=grid(),ylim=c(0,maxy))
  if (polys) {
    for (i in 1:nscen) {
      if (q90) {
        poldat <- makepolygon(varq[2,,i],varq[4,,i],yrnames)
        polygon(poldat,col=RGB(scencol[i],alpha=intens))
      } else {
        poldat <- makepolygon(varq[1,,i],varq[5,,i],yrnames)
        polygon(poldat,col=RGB(scencol[i],alpha=intens))
      }
      if (!is.null(hlin)) abline(h=hlin[i],lwd=2,col=scencol[i])
    }
  } else {
    for (i in 1:nscen) {
      lines(yrnames,varq[3,,i],lwd=2,col=scencol[i])
      if (q90) {
        lines(yrnames,varq[2,,i],lwd=1,col=scencol[i])
        lines(yrnames,varq[4,,i],lwd=1,col=scencol[i])
      } else {
        lines(yrnames,varq[1,,i],lwd=1,col=scencol[i])
        lines(yrnames,varq[5,,i],lwd=1,col=scencol[i])
      }
    }
  } # end of polys if
  if (addleg) {
    legend(locate,legend=scenes,col=scencol,lwd=3,bty="n",cex=1.1)
  }
} # end of doquantplot

#' @title doquantribbon generates a ribbon of multiple quantile values
#'
#' @description doquantribbon is used when trying to present a summary of multiple
#'     scenarios applied across a sequence of years (for example in the
#'     projection years of an MSE). It assumes that for each scenario quantiles
#'     have been calculated for each year so that the input array if 5 x nyrs x
#'     nscenes, where nyrs is the number of years and nscenes is the number of
#'     scenarios. Five quantiles are assumed with values
#'     probs=c(0.025, 0.05, 0.5, 0.95 0.975), and the option exists of plotting
#'     either the 90th (0.05 and 0.95) or the 95th (0.025, 0.975). If the input
#'     quantiles have different values then the 'q90' argument,if TRUE, will
#'     reference the 2nd and 4th row and the 1st and 5th if FALSE. If polys is
#'     TRUE then filled transparent polygons will be plotted, with colours in
#'     the same sequence as the scenarios, alternatively, if polys=FALSE, then
#'     lines will be plotted instead. This function is currently only used by
#'     sauribbon.
#'
#' @param varq an array of quantiles with dimensions the 5 quantiles of the
#'     variable in question (2.5, 5, 50, 95 97.5), by years, the number of
#'     years for which quantiles are available
#' @param varname The names of the variable being summarized, used as the
#'     y-label on the plot
#' @param sauname the name of the sau, to add to the y-axis label
#' @param yrnames the numeric values of the years to be plotted. Used as the
#'     x-axis as well as the x-axis labels
#' @param scenes the names given to the different scenarios
#' @param q90 should the 90th or 95th quantiles be plotted, default=TRUE = 90th
#' @param intens signifies the intensity of colour on a scale of 0 - 255.
#'     127 is about 50 percent dense.
#' @param addleg add a legend? currently can only use topleft, topright,
#'     bottomleft, or bottomright = default
#' @param addmedian should a solid line depicting the median for each scenario
#'     be included on each polygon. default = FALSE
#'
#' @seealso {
#'    \link{sauribbon}, \link{RGB}
#' }
#'
#' @return nothing but it does add polygons to a plot
#' @export
#'
#' @examples
#' print("wait on datasets")
#' # varq=qntvar;varname=varname;yrnames=glbc[[1]]$pyrnames;scenes=scenes
#' # q90=q90;intens=intens;addleg=addleg;addmedian=addmedian
doquantribbon <- function(varq,varname,sauname,yrnames,scenes,q90=TRUE,intens=127,
                          addleg="bottomright",addmedian=0) {
  maxy <- getmax(varq)
  nscen <- length(scenes)
  plot(yrnames,varq[3,,1],type="l",lwd=2,col=0,xlab="",
       ylab=paste0(varname,"_",sauname),panel.first=grid(),ylim=c(0,maxy))
  for (i in 1:nscen) { # i =2
    if (q90) {
      poldat <- makepolygon(varq[2,,i],varq[4,,i],yrnames)
      polygon(poldat,col=RGB(i,alpha=intens))
      if (addmedian > 0) lines(yrnames,varq[3,,i],lwd=addmedian,col=i)
    } else {
      poldat <- makepolygon(varq[1,,i],varq[5,,i],yrnames)
      polygon(poldat,col=RGB(i,alpha=intens))
      if (addmedian > 0) lines(yrnames,varq[3,,i],lwd=addmedian,col=i)
    }
  }
  if (nchar(addleg) > 0) {
    if (addleg %in% c("topleft","topright","bottomleft","bottomright")) {
      legend(addleg,legend=scenes,col=1:nscen,lwd=3,bty="n",cex=1.1)
    } else {
      cat("Inappropriate legend location given, see ?doquantribbon \n")
    }
  }
} # end of doquantribbon

#' @title getcompout extracts a component by scenario from output of do_comparison
#'
#' @description getcompout is used to extract a given variable from the dynamics
#'     included in the output from the do_comparison function. It first isolates
#'     the 'dyn' component from the 'result' object, then extracts and collates
#'     the selected variable. There are ten valid variable names that this
#'     function will accept: matureB, exploitB, midyexpB, catch, acatch,
#'     harvestR, cpue, recruit, deplsB, and depleB. Anything else will throw an
#'     error message. The outputs can be used by plotdynphase.
#'
#' @param dyn a list of the dynamics objects from each scenbario
#' @param glb the globals object. The assumption is that all glb objects are
#'     very similar between scenarios, although different reps can be compared
#' @param scenes the names of the scenarios selected for comparison
#' @param pickvar the name of the variable of interest, must be in quotes.
#' @param projonly should only the projection years be used, default=TRUE
#'
#' @seealso [\link{plotdynphase} for how this may be used]
#'
#' @return a list of the selected variable for each scenario. Each list
#'     component will be a 3-D array of years x sau x replicate
#' @export
#'
#' @examples
#' print("wait on example data")
getcompout <- function(dyn,glb,scenes,pickvar="matureB",projonly=TRUE) { # out=result; pickvar="matureB"
  validvars <- c("matureB","exploitB","midyexpB","catch","acatch","harvestR",
                 "cpue","recruit","deplsB","depleB")
  if (pickvar %in% validvars) {
    nscen <- length(scenes)
    ans <- vector(mode="list",length=nscen)
    names(ans) <- scenes
    for (i in 1:nscen) {
      if (projonly) {
        yrs <- glb$hyrs:(glb$hyrs + glb$pyrs)
        ans[[i]] <- dyn[[i]][[pickvar]][yrs,,]
      } else {
        ans[[i]] <- dyn[[i]][[pickvar]]
      }
    }
  } else {
    cat("Invalid input, pick from \n",validvars," \n")
    ans <- NULL
  }
  return(ans)
} # end of getcompout

#' @title getabdevs absolute deviations between a vector and its moving average
#'
#' @description getabdevs is used to calculate the absolute deviations between
#'     a vector of numbers and either a loess fit to the vector or a moving
#'     average across that vector. This provides a measure of variation akin
#'     to a Average Annual Variation (aav) suitable for oscillating projections.
#'     This function can be used with apply to summarize across replicates.
#'
#' @param x a vector of numbers to be tested, perhaps a single replicate
#'     projection
#' @param years the loess or moving average is relative to a set of years or
#'     some other sequence. It must hve the same length as x
#' @param smooth default = TRUE, which implies a standard loess curve will be
#'     fitted to the curve before finding the deviations. If set FALSE, then
#'     a moving average of length n will be used.
#' @param n default = 7, the length of a moving average if one is used.
#'
#' @return a list of the avdev and sumdev across the vector
#' @export
#'
#' @examples
#' x <- runif(98,min=1,max=25)
#' outdev <- getabdevs(x,1:98)
getabdevs <- function(x, years,smooth=TRUE, n = 7) {
  if (smooth) {
    ma <- loess(x ~ years)
    pick <- which(!is.na(ma))
    sumdev <- sum(abs(x[pick] - ma$fitted[pick]))
  } else {
    ma <- movav(x,n=n)
    pick <- which(!is.na(ma))
    sumdev <- sum(abs(x[pick] - ma[pick]))
  }
  return(sumdev=sumdev)
} # end of getdevs

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
#' @param scenarionames a vector of the original scenarionames
#' @param postdir the name of the directory holding the results, also used to
#'     name the internal webpage
#' @param filesused a vector of the names of the scenarios that were compared,
#'     this would usually be files[pickfiles]. Added to the home page.
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
#' # rundir=rundir;glbc=glbc;scenes=scenes;postdir=postfixdir
#' # filesused=files[c(1,2,3,4)];openfile=TRUE;verbose=FALSE
makecompareoutput <- function(rundir,glbc,scenes,scenarionames,postdir,
                              filesused,openfile=TRUE,verbose=FALSE) {
  replist <- list(starttime=as.character(Sys.time()),
                  endtime="")
  nscene <- length(scenes)
  reps <- NULL
  for (i in 1:nscene) reps <- paste0(reps,glbc[[i]]$reps,"  ")
  if (all.equal.character(scenarionames,scenes) != TRUE) {
    scenarios <- paste0(scenarionames,"_",scenes)
  } else {
    scenarios <- scenes
  }
  runnotes <- c("Scenarios: ",
                scenarios,
                paste0("RunTime = ",0),
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

#' @title makelistobj extracts a list of objects from a bigger list of objects
#'
#' @description makelistobj when making comparisons between scenarios first it
#'     is necessary to gather information from each scenario into large lists.
#'     For example, after running d the do_comparison function the output is
#'     a list which includes an obj called 'ans'. This contains the output from
#'     do_MSE for each scenario. Inside that are objects such as glb, ctrl,
#'     zoneDP, sauout, and so on. For further analysis we need to be able to
#'     extract copies of some of these objects into lists of their own. The
#'     function makelistobj can do that.
#'
#' @param res the list object from which individual objects are wanted from each
#'     scenario
#' @param objname the name of the object wanted from each scenario
#'
#' @seealso{
#'     \link{do_comparison}, \link{do_MSE}
#' }
#'
#' @return a list of the named objects from res
#' @export
#'
#' @examples
#' \dontrun{
#'   # see the help for do_comparison for code lines prior to the next line
#'   result <- do_comparison(rundir,postfixdir,outdir,files,pickfiles=c(1,2,3))
#'   glball <- makelistobj(res=result$ans,objname="glb")
#'   ctrlall <- makelistobj(res=result$ans,objname="ctrl")
#'   str(glball)
#' }
makelistobj <- function(res,objname) { # res=result$ans;objname="glb"
  scenes <- names(res)
  nscen <- length(scenes)
  outlist <- makelist(scenes)
  for (i in 1:nscen) outlist[[i]] <- res[[scenes[i]]][[objname]]
  return(invisible(outlist))
} # end of makelistobj

#' @title movav generates a moving average of a vector
#'
#' @description movav generates a moving average of a vector of numbers. It
#'     allows the length of the moving average to be set using 'n', and
#'     whether the moving average starts at index 1 or is centred across the
#'     vector, using 'centred'. Only odd numbered lengths for the moving
#'     average are permitted to allow centring to be balanced. The output is
#'     originall a time-series but is converted to a numeric vector.
#'
#' @param x a vector of numbers for which a moving average is wanted
#' @param n the length of the moving average, default = 5
#' @param centred should the moving average be centred across the vector or
#'     start at index = 1, default = centred=TRUE
#'
#' @return a numeric vector containing a moving average of the input vector
#'     or a printed warning if the length is an even number.
#' @export
#'
#' @examples
#' x <- runif(25,min=1,max=25)
#' outav <- movav(x,n=5)
#' cbind(x,outav)
movav <- function(x, n = 5,centred=TRUE) {
  if ((n %% 2) == 0) {
    print("Only use odd numbers for moving averages, NULL returned. \n")
    return(NULL)
  } else {
    loc <- ifelse(centred,2,1)
    ma <- as.numeric(stats::filter(x, rep(1/n, n), sides = 2))
    return(ma)
  }
} # end of movav

#' @title plotallphaseplots generate a set of phase plots of dynamics variables
#'
#' @description plotallphaseplots produced a series of phase plots for a set of
#'     different scenarios of how mature biomass and its depletion varies with
#'     aspirational catch along with other combinations of cpue and harvest
#'     rate against catch.
#'
#' @param rundir the compare scenario's directory
#' @param dyn a list of the dynamics objects from each scenario
#' @param prods a list of the sauprod objects containing MSY, B0, Bmsy, etc
#' @param glb the globals object. The assumption is that all glb objects are
#'     very similar between scenarios, although different reps can be compared
#' @param scenes the names of the scenarios selected for comparison
#' @param width width of the plot. default=9
#' @param height height of the plot. default=10
#' @param fnt which font to use. default=7 = bold times
#' @param pntcex the size of the year points along each median line
#' @param zero should the phase plots have a zero origin. default=FALSE
#' @param legloc location of the legend placed into the last plot.
#'     default='bottomright'
#'
#' @seealso [\link{plotallphaseplots} which wraps a set of phase plots]
#'
#' @return nothing
#' @export
#'
#' @examples
#' print("wait on suitable data")
plotallphaseplots <- function(rundir,dyn,prods,glb,scenes,width=9,height=10,
                              fnt=7,pntcex=1.5,zero=FALSE,legloc="bottomright") {
  # rundir=rundir; dyn=dyn;prods=prods;glb=glbc[[1]];scenes=scenes;width=9;height=10;
  # fnt=7; pntcex=1.5; zero=FALSE; legloc="topright"
  matureB <- getcompout(dyn=dyn,glb=glb,scenes=scenes,pickvar="matureB")
  exploitB <- getcompout(dyn=dyn,glb=glb,scenes=scenes,pickvar="exploitB")
  acatch <- getcompout(dyn=dyn,glb=glb,scenes=scenes,pickvar="acatch")
  catch <- getcompout(dyn=dyn,glb=glb,scenes=scenes,pickvar="catch")
  deplsB <- getcompout(dyn=dyn,glb=glb,scenes=scenes,pickvar="deplsB")
  depleB <- getcompout(dyn=dyn,glb=glb,scenes=scenes,pickvar="depleB")
  cpue <- getcompout(dyn=dyn,glb=glb,scenes=scenes,pickvar="cpue")
  harvestR <- getcompout(dyn=dyn,glb=glb,scenes=scenes,pickvar="harvestR")
  msy <- getmatcolfromlist(prods,"MSY")
  Bmsy <- getmatcolfromlist(prods,"Bmsy")
  Bexmsy <- getmatcolfromlist(prods,"Bexmsy")
  catbymsy <- catch
  matBbyBmsy <- matureB
  expBbyBexmsy <- exploitB
  for (scene in 1:length(scenes)) {
    for (sau in 1:glb$nSAU) {
      catbymsy[[scene]][,sau,] <- (catch[[scene]][,sau,])/msy[[scene]][sau]
      matBbyBmsy[[scene]][,sau,] <- matureB[[scene]][,sau,]/Bmsy[[scene]][sau]
      expBbyBexmsy[[scene]][,sau,] <- exploitB[[scene]][,sau,]/Bexmsy[[scene]][sau]
    }
  }
  fileout <- plotdynphase(xlist=depleB,ylist=catbymsy,scenes=scenes,glb=glb,
                          rundir=rundir,xlab="Exploitable Biomass Depletion",
                          ylab="Actual Catches div MSY",legloc=legloc,
                          console=FALSE,width=width,height=height,fnt=fnt,
                          pntcex=pntcex,zero=zero,yline=1.0)
  addplot(filen=fileout$filen,rundir=rundir,category="phaseplots",
          caption=fileout$caption)
  fileout <- plotdynphase(xlist=deplsB,ylist=catbymsy,scenes=scenes,glb=glb,
                          rundir=rundir,xlab="Mature Biomass Depletion",
                          ylab="Actual Catches div MSY",legloc=legloc,
                          console=FALSE,width=width,height=height,fnt=fnt,
                          pntcex=pntcex,zero=zero,yline=1.0)
  addplot(filen=fileout$filen,rundir=rundir,category="phaseplots",
          caption=fileout$caption)
  fileout <- plotdynphase(xlist=matBbyBmsy,ylist=catbymsy,scenes=scenes,glb=glb,
                          rundir=rundir,xlab="Mature Biomass div Bmsy",
                          ylab="Actual Catches div MSY",legloc=legloc,
                          console=FALSE,width=width,height=height,fnt=fnt,
                          pntcex=pntcex,zero=zero,yline=1.0,hline=1.0)
  addplot(filen=fileout$filen,rundir=rundir,category="phaseplots",
          caption=fileout$caption)
  fileout <- plotdynphase(xlist=expBbyBexmsy,ylist=catbymsy,scenes=scenes,glb=glb,
                          rundir=rundir,xlab="Exploitable Biomass div Bexmsy",
                          ylab="Actual Catches div MSY",legloc=legloc,
                          console=FALSE,width=width,height=height,fnt=fnt,
                          pntcex=pntcex,zero=zero,yline=1.0,hline=1.0)
  addplot(filen=fileout$filen,rundir=rundir,category="phaseplots",
          caption=fileout$caption)
  fileout <- plotdynphase(xlist=depleB,ylist=catch,scenes=scenes,glb=glb,
                          rundir=rundir,xlab="Exploitable Biomass Depletion",
                          ylab="Actual Catches (t)",legloc=legloc,
                          console=FALSE,width=width,height=height,fnt=fnt,
                          pntcex=pntcex,zero=zero)
  addplot(filen=fileout$filen,rundir=rundir,category="phaseplots",
          caption=fileout$caption)
  fileout <- plotdynphase(xlist=deplsB,ylist=catch,scenes=scenes,glb=glb,
                          rundir=rundir,xlab="Mature Biomass Depletion",
                          ylab="Actual Catches (t)",legloc=legloc,
                          console=FALSE,width=width,height=height,fnt=fnt,
                          pntcex=pntcex,zero=zero)
  addplot(filen=fileout$filen,rundir=rundir,category="phaseplots",
          caption=fileout$caption)
  fileout <- plotdynphase(xlist=matureB,ylist=acatch,scenes=scenes,glb=glb,
                          rundir=rundir,xlab="Mature Biomass (t)",
                          ylab="Aspirational Catch (t)",legloc=legloc,
                          console=FALSE,width=width,height=height,fnt=fnt,
                          pntcex=pntcex,zero=zero)
  addplot(filen=fileout$filen,rundir=rundir,category="phaseplots",
          caption=fileout$caption)
  fileout <- plotdynphase(xlist=deplsB,ylist=acatch,scenes=scenes,glb=glb,
                          rundir=rundir,xlab="Depletion of Mature Biomass (t)",
                          ylab="Aspirational Catch (t)",
                          legloc=legloc,console=FALSE,width=width,
                          height=height,fnt=fnt,pntcex=pntcex,zero=zero)
  addplot(filen=fileout$filen,rundir=rundir,category="phaseplots",
          caption=fileout$caption)
  fileout <- plotdynphase(xlist=catch,ylist=cpue,scenes=scenes,glb=glb,
                          rundir=rundir,xlab="Actual Catch (t)",
                          ylab="CPUE",legloc=legloc,console=FALSE,
                          width=width,height=height,fnt=fnt,
                          pntcex=pntcex,zero=zero)
  addplot(filen=fileout$filen,rundir=rundir,category="phaseplots",
          caption=fileout$caption)

  fileout <- plotdynphase(xlist=catch,ylist=harvestR,scenes=scenes,glb=glb,
                          rundir=rundir,xlab="Actual Catch (t)",
                          ylab="Annual Harvest Rate",legloc=legloc,
                          console=FALSE,width=width,height=height,fnt=fnt,
                          pntcex=pntcex,zero=zero)
  addplot(filen=fileout$filen,rundir=rundir,category="phaseplots",
          caption=fileout$caption)
  fileout <- plotdynphase(xlist=exploitB,ylist=catch,scenes=scenes,glb=glb,
                          rundir=rundir,xlab="Exploitable Biomass",
                          ylab="Actual Catch",legloc=legloc,
                          console=FALSE,width=width,height=height,fnt=fnt,
                          pntcex=pntcex,zero=zero)
  addplot(filen=fileout$filen,rundir=rundir,category="phaseplots",
          caption=fileout$caption)
} # end of plotphaseplots

#' @title plotdevs plots projection deviates for a variable across scenerios
#'
#' @description plotdevs provides a plot of a scenario Performance Measure. In
#'     this case the PM are the deviations from a loess fit to  whatever input
#'     variable has been selected from the dyn objects in do_comparison(). One
#'     needs first generate a list of the selected variable from the the input
#'     dyn objects (one from each scenario). A standard comparison would be on
#'     the variation of catches through the years of the projection. But this
#'     could be applied to any of the variables within the dyn object, which
#'     includes matureB, exploitB, midyexpB, catch, acatch, harvestR, cpue,
#'     recruit, deplsB, and depleB.
#'
#' @param invar the sub-object from each scenario as a list, see examples
#' @param scenes a vector of names fr each scenario being compared
#' @param saunames a vector of the names of each sau in the projections
#' @param filen the filename and path if the plot is to be saved, default = "".
#'     if a filename is given it should be a .png file
#'
#' @return a list f the mean deviates per sau and scenario, and a list of their
#'     respective standard deviations. It also plots a graph
#' @export
#'
#' @examples
#' \dontrun{
#' # example syntax, assumes dynamics have been extracted from all scenarios
#'    glb <- glbc[[1]]
#'    saunames <- glb$saunames
#'    nsau <- glb$nSAU
#'    nscen <- length(scenes)
#'    x <- makelist(scenes)
#'    for (scen in 1:length(scenes)) x[[scen]] <- dyn[[scen]]$catch[59:88,,]
#'    devout <- plotdevs(x,scenes,saunames,filen="")
#' }
plotdevs <- function(invar,scenes,saunames,filen=""){
  nsau <- length(saunames)
  nscen <- length(scenes)
  meandevs <- matrix(0,nrow=nscen,ncol=nsau,dimnames=list(scenes,saunames))
  sddevs <- meandevs
  tmp <- invar[[1]][,1,]
  yrs <- as.numeric(rownames(tmp))
  plotprep(width=12,height=7,newdev=FALSE,filename=filen,verbose=FALSE)
  parset(plots=c(nscen,nsau),margin=c(0.25,0.25,0.05,0.05),
         outmargin = c(1.2,1.2,0.1,0.1),byrow=TRUE)
  for (scen in 1:nscen) {
    for (sau in 1:nsau) {
      # scen = 1; sau=1
      sauC <- invar[[scen]][,sau,]
      devs <- apply(sauC,2,getabdevs,years=yrs,smooth=TRUE)
      avdevs <- mean(devs,na.rm=TRUE)
      meandevs[scen,sau] <- avdevs
      sddevs[scen,sau] <- sd(devs,na.rm=TRUE)
      hist(devs,breaks=20,main="",xlab="",ylab="")
      mtext(round(avdevs,2),side=3,line=-1.1,cex=0.9,adj=1)
      if (sau == 1) mtext(scenes[scen],side=2,line=1.25,cex=1.0)
      if (scen == nscen) mtext(saunames[sau],side=1,line=1.25,cex=1)
    }
  }
  return(list(meandevs=meandevs,sddevs=sddevs))
} # end of plotdevs

#' @title plotdynphase plots the medians of 2 dynamic variables as a phase diagram
#'
#' @description plotdynphase is used when comparing how the population dynamics
#'      are affected by multiple scenarios. The median of the replicate
#'      trajectories for each selected variable are plotted against each other
#'      in a classical phase plot. Thus, a stable solution should display a
#'      single point or an oscillation around a point of convergence. Each plot
#'      can have an origin at zero or be scaled to display maximum detail using
#'      the 'zero' argument. As with getcompout there are only ten valid variable
#'      names that can be combined: matureB, exploitB, midyexpB, catch, acatch,
#'      harvestR, cpue, recruit, deplsB, and depleB. The results for each SAU
#'      are plotted together.
#'
#' @param xlist x-axis list of each scenes replicates for the selected variable,
#'     generated from result using getcompout - 'get component out'
#' @param ylist y-axis list of each scenes replicates for the selected variable
#' @param scenes names of the difference scenes being compared. result$scenes
#' @param glb the aMSE globals object, inside one of the result$ans lists
#' @param rundir the directory into which all files and results are placed
#' @param xlab label for the x-axis
#' @param ylab label for the y-axis
#' @param legloc location of the legend placed into the last plot.
#'     default='bottomright'
#' @param console should the plot go to the console or a png file. default=TRUE
#' @param width width of the plot. default=9
#' @param height height of the plot. default=10
#' @param fnt which font to use. default=7 = bold times
#' @param pntcex the size of the year points along each median line
#' @param zero should the phase plots have a zero origin. default=FALSE
#' @param yline default = NULL, used to add horizontal lines to each plot
#' @param hline default = NULL, used to add vertical lines to each plot
#'
#' @seealso [\link{getcompout} for how to generate the xlist and ylist]
#'
#' @return an invisible caption to the plot.
#' @export
#'
#' @examples
#' print("wait on suitable data sets")
#' # result <- do_comparison(rundir=rundir,postfixdir=postfixdir,outdir="",
#' #                         files=fullfiles,pickfiles=c(1,4),verbose=TRUE,
#' #                         intensity=100)
#' # matureB <- getcompout(result,pickvar="matureB")
#' # acatch <- getcompout(result,pickvar="acatch")
#' # glb <- result$ans[[1]]$glb
#' # scenes <- result$scenes
#' # plotdynphase(xlist=acatch,ylist=matureB,scenes=scenes,glb=glb,rundir="",
#' #              xlab="Aspirational Catch (t)",ylab="Mature Biomass (t)",
#' #              legloc="bottomright",console=TRUE,width=9,height=10,fnt=7,
#' #              pntcex=1.5,zero=FALSE)
plotdynphase <- function(xlist,ylist,scenes,glb,rundir="",xlab="xlabel",
                         ylab="ylabel",legloc="bottomright",console=TRUE,
                         width=9,height=10,fnt=7,pntcex=1.25,zero=FALSE,
                         yline=NULL,hline=NULL) {
# xlist=depleB;ylist=catbymsy;scenes=scenes;glb=glb;rundir=rundir;
# xlab="Exploitable Biomass Depletion";ylab="Actual Catches div MSY";
# legloc=legloc; console=FALSE;width=width;height=height;fnt=fnt
# pntcex=pntcex;zero=zero;yline=1.0; hline=140
  nscenes <- length(scenes)
  nsau <- glb$nSAU
  saunames <- glb$saunames
  yrnames <- c(tail(glb$hyrnames,1),glb$pyrnames)
  yvar <- xvar <- array(0,c(glb$pyrs+1,nsau,nscenes),
                        dimnames=list(yrnames,glb$saunames,scenes))
  for (scen in 1:nscenes) {
    for (sau in 1:nsau) {
      yvar[,sau,scen] <- apply(ylist[[scen]][,sau,],1,median)
      xvar[,sau,scen] <- apply(xlist[[scen]][,sau,],1,median)
    }
  }
  yrs <- as.numeric(yrnames)
  if (console) { filen <- "" } else {
    filen <- pathtopath(rundir,
                        paste0(removeEmpty(ylab),"_vs_",removeEmpty(xlab),".png"))
  }
  plotprep(width=width,height=height,newdev=FALSE,filename=filen,
           verbose=FALSE,usefont=fnt)
  parset(plots=pickbound(nsau),outmargin=c(1,1,0,0),margin=c(0.3,0.3,0.025,0.1))
  for (sau in 1:nsau) {
    nyr <- length(xvar[,1,1])
    pts <- seq(5,nyr,5)
    maxx <- getmax(xvar[,sau,])
    minx <- getmin(xvar[,sau,])
    maxy <- getmax(yvar[,sau,])
    miny <- getmin(yvar[,sau,])
    if (zero) {
      minx <- 0
      miny <- 0
    }
    plot(xvar[,sau,1],yvar[,sau,1],type="l",lwd=2,xlim=c(minx,maxx),
         ylim=c(miny,maxy),panel.first=grid(),xlab="",ylab="")
    points(xvar[pts,sau,1],yvar[pts,sau,1],pch=16,cex=pntcex,col=1)
    points(xvar[c(1,nyr),sau,1],yvar[c(1,nyr),sau,1],pch=16,cex=1.5*pntcex,
           col=c(4,3))
    points(xvar[nyr,sau,1],yvar[nyr,sau,1],pch=16,cex=1.0,col=1.25)
    for (scen in 2:nscenes) {
      lines(xvar[,sau,scen],yvar[,sau,scen],lwd=2,col=scen)
      points(xvar[pts,sau,scen],yvar[pts,sau,scen],pch=16,cex=pntcex,col=scen)
      points(xvar[c(1,nyr),sau,scen],yvar[c(1,nyr),sau,scen],pch=16,
             cex=1.5*pntcex,col=c(4,3))
      points(xvar[nyr,sau,scen],yvar[nyr,sau,scen],pch=16,cex=1.25,col=scen)
    }
    if (length(yline) > 0) abline(h=yline,lwd=1,col=1)
    if (length(hline) > 0) abline(v=hline,lwd=1,col=1)
    mtext(saunames[sau],side=3,outer=FALSE,cex=1.1,line=-1.1,adj=0.05)
  }
  mtext(xlab,side=1,outer=TRUE,cex=1.1,line=-0.2)
  mtext(ylab,side=2,outer=TRUE,cex=1.1,line=-0.2)
  legend(legloc,legend=scenes,col=1:nscenes,lwd=3,bty="n",cex=1.1)
  caption=""
  if (!console) {
    caption <- paste0("Dynamic variables ",xlab," and ",ylab,
                      " compared for each SAU.")
    dev.off()
  }
  return(invisible(list(caption=caption,filen=filen)))
} # end of plotdynphase

#' @title plotdynvarinyear  histograms of a dyn variable by scenario and sau
#'
#' @description plotdynvarinyear generates a set of histograms of a user
#'     select variable within the dyn object from do_comparison, all with the
#'     same x-axis, for a given year by scenario and SAU. The number of breaks
#'     can be modified but defaults to 15. Only a single year at a time is
#'     characterized.
#'
#' @param rundir the directory in which comparisons are being made. It is best
#'     to be a separate directory form any particular scenario.
#' @param dyn the dynamics objects from each scenario in a single dyn object
#'     from the output of do_comparison
#' @param whichvar which variable within the dyn object to characterize
#' @param whichyr which year number to use
#' @param glb the globals object needed for the number and names of the sau,
#'     and the years and their names.
#' @param bins how make breaks should be used in the histograms, default=15
#' @param console should each plot go to the console or be saved to rundir.
#'     default = TRUE ie go to console
#' @param whichfun each plot has by default the mean of the distribution added
#'     if one wanted the median or geomean then change to that function name
#'
#' @return the filename (no path) used if any, and it plots a graph
#' @export
#'
#' @examples
#' # assumes that do_comparison function has been run. syntax used:
#' # filen <- plotdynvarinyear(rundir=rundir,dyn=dyn,whichvar="catch",
#' #                           whichyr=c(63),glb,console=TRUE)
plotdynvarinyear <- function(rundir,dyn,whichvar,whichyr,glb,bins=15,
                             console=TRUE,whichfun=mean) {
# rundir=rundir; dyn=dyn;whichvar="catch";whichyr=c(63);glb=glb;bins=15;
# console=TRUE; whichfun=mean
  years <- c(glb$hyrnames,glb$pyrnames)
  nsau <- glb$nSAU
  saunames=glb$saunames
  scenes <- names(dyn)
  nscen <- length(scenes)
  indexvar <- match(whichvar,names(dyn[[1]]))
  res <- state_by_year(indyn=dyn,whichvar=whichvar,whichyrs=whichyr,
                       whichfun=whichfun,glb=glb)
  xranges <- matrix(0,nrow=nsau,ncol=2,dimnames=list(saunames,c("low","high")))
  xranges[,1] <- 1e7; xranges[,2] <- -1e7
  for (sau in 1:nsau) {
    for (scen in 1:nscen) {
      tmp <- range(dyn[[scen]][[indexvar]][whichyr,sau,])
      xranges[sau,1] <- ifelse(xranges[sau,1] < tmp[1],xranges[sau,1],tmp[1])
      xranges[sau,2] <- ifelse(xranges[sau,2] > tmp[2],xranges[sau,2],tmp[2])
    }
  }
  xranges[,1] <- floor(xranges[,1])
  xranges[,2] <- ceiling(xranges[,2])
  fname <- paste0(whichvar,"_distribution_in_",years[whichyr])
  fname1 <- paste0(fname,".png")
  filen <- pathtopath(rundir,fname1)
  if (console) filen=""
  plotprep(width=12,height=7,newdev=FALSE,filename=filen,verbose=FALSE)
  parset(plots=c(nscen,nsau),margin=c(0.25,0.25,0.05,0.05),
         outmargin = c(1.2,2.5,0.1,0.1),byrow=TRUE)
  for (scen in 1:nscen) {
    for (sau in 1:nsau) {
      # yr = 1; scen = 1; sau=1
      getyr <- dyn[[scen]][[indexvar]][whichyr,sau,]
      hist(getyr,breaks=bins,main="",xlab="",ylab="",xlim=xranges[sau,])
      abline(v=res[1,scen,sau],lwd=3,col=4)
        if (sau == 1) mtext(scenes[scen],side=2,line=1.25,cex=1.0)
        if (scen == nscen) mtext(saunames[sau],side=1,line=1.25,cex=1)
    }
  }
  mtext(fname,side=2,outer=TRUE,line=1,cex=1.5)
  return(invisible(fname1))
} # end of plotdynvarinyear

#' @title plotrateofchange plot the output from getrateofchange function
#'
#' @description plotrateofchange produces a plot of the percent changes within
#'     the whichvar selected across all selected scenarios within do_comparison.
#'     The outcome is placed within the ScenarioPMs tab. The variables it can
#'     work with from dyn from teh output of do_comparison include: matureB,
#'     exploitB, midyrexpB, catch, acatch, harvestR, cpue, recruit, deplsB,
#'     and depleB.
#'
#' @param rundir the directory in which comparisons are being made. It is best
#'     to be a separate directory form any particular scenario.
#' @param res the output object from the getrateofchange function
#' @param whichvar which variable within the dyn object to characterize
#' @param glb the globals object needed for the number and names of the sau,
#'     and the years and their names.
#' @param console should each plot go to the console or be saved to rundir.
#'     default = TRUE ie go to console
#'
#' @seealso \link{getrateofchange}, \link{do_comparison}
#'
#' @returns invisibly the filename used without the path
#' @export
#'
#' @examples
#' # syntax
#' # pickvar <- "acatch"
#' # res <- getrateofchange(dyn=dyn,whichvar=pickvar,glb=glb)
#' # filen <- plotrateofchange(rundir=rundir,res=res,whichvar=pickvar,glb=glb)
plotrateofchange <- function(rundir,res,whichvar,glb,console=TRUE) {
  scenes <- names(res)
  nscen <- length(scenes)
  saunames <- glb$saunames
  nsau <- glb$nSAU
  projyrs <- glb$pyrnames
  pyrs <- glb$pyrs
  outy <- matrix(0,nrow=nsau,ncol=2,dimnames=list(saunames,c("low","high")))
  outy[,1] <- 1e6
  outy[,2] <- -1e6
  for (sau in 1:nsau) {
    for (scen in 1:nscen) {
      tmp <- res[[scen]]$pdiffer[,sau]
      yrge <- range(tmp)
      outy[sau,1] <- ifelse((yrge[1] < outy[sau,1]),yrge[1],outy[sau,1])
      outy[sau,2] <- ifelse((yrge[2] > outy[sau,2]),yrge[2],outy[sau,2])
    }
  }
  outy[,1] <- floor(outy[,1])
  outy[,2] <- ceiling(outy[,2])
  outy
  yrs <- as.numeric(projyrs[2:pyrs])
  fname <- paste0(whichvar,"_rate_of_change_across_scenarios")
  fname1 <- paste0(fname,".png")
  filen <- pathtopath(rundir,fname1)
  if (console) filen=""
  plotprep(width=9,height=8,newdev=FALSE,filename=filen,verbose=FALSE)
  parset(plots=pickbound(nsau),margin=c(0.2,0.4,0.1,0.1),outmargin=c(0,1,0,0),
         byrow=FALSE)
  for (sau in 1:nsau) {
    tmp <- res[[1]]$pdiffer[,sau]
    plot(yrs,tmp,type="l",lwd=2,col=1,ylab=saunames[sau],ylim=outy[sau,],xlab="")
    abline(h=0,lwd=1,col=1)
    for (scen in 2:nscen) {
      tmp <- res[[scen]]$pdiffer[,sau]
      lines(yrs,tmp,lwd=2,col=scen)
    }
    if (sau == 1) legend("topright",scenes,lwd=3,col=1:nscen,bty="n",cex=1)
  }
  label <- paste0("Percentage Change in median ",whichvar," across Scenarios.")
  mtext(label,side=2,outer=TRUE,cex=1.2,line=-0.5)
  return(invisible(fname1))
} # end of plotrateofchange

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

#' @title plotscenproj plots a 3D array of data by sau in a single plot
#'
#' @description plotsceneproj uses a 3D array of years x sau x replicates to
#'     produce a plot by sau of those projections all in a single plot. Options
#'     exist to include a horizontal dashed line and to include a median line
#'     and either 90 or 95 quantiles across the plot.
#'
#' @param rundir the directory in which all results are held for a scenario or
#'     comparison of scenarios
#' @param inarr the 3D array of projections for a given scenario derived from
#'     catch or cpue or whatever
#' @param glb the global object relating to the particular acenario
#' @param scene the scenario name
#' @param filen the filename in which to store the plot, default="" which draws
#'     the plot to the console
#' @param label what is the name of the variable being plotted, default=""
#' @param maxy what constant maximum y-axis to use? default=0 which uss the
#'     maximum for each sau
#' @param Q which quantile to use, default = 90. If set to 0 no median or
#'     quantiles are plotted, 95 will lead to the 95 percent quantiles being
#'     plotted
#' @param hline should a horizontal dashed line be included. default=NA, which
#'     means that no line is added. Otherwise whichever value is given this
#'     argument will lead to a dashed horizontal black line of width 1.
#'
#' @return invisibly a list of the quantiles for each sau
#' @export
#'
#' @examples
#' print("wait on data sets")
#' #  rundir=rundir; inarr=cdivmsy[[1]];glb=glbc[[1]];scenes=scenes[1];
#' #  filen="";label="Catch / MSY";maxy=0
plotsceneproj <- function(rundir,inarr,glb,scene,filen="",label="",
                          maxy=0,Q=90,hline=NA) {
  nsau <- glb$nSAU
  saunames <- glb$saunames
  outmed <- makelist(saunames)
  reps <- glb$reps
  if (nchar(filen) > 0) {
    filen <- filenametopath(rundir,filen)
    caption <- paste0("Projections of ",label," for ",scene," for each SAU.")
  }
  yrs <- as.numeric(names(inarr[,1,1]))
  plotprep(width=8, height=8,newdev=FALSE,filename = filen,verbose=FALSE)
  parset(plots=pickbound(nsau),margin=c(0.3,0.4,0.05,0.1),outmargin=c(0,1,0.25,0),
         byrow=FALSE)
  for (sau in 1:nsau) { # sau=1; i = 1
    dat <- inarr[,sau,]
    meds <- apply(dat,1,quants)
    outmed[[sau]] <- meds
    ymax <- maxy
    if (ymax == 0) ymax <- getmax(dat)
    plot(yrs,dat[,1],type="l",lwd=1,col="grey",panel.first=grid(),
         ylim=c(0,ymax),ylab=saunames[sau],xlab="")
    for (i in 1:reps) lines(yrs,dat[,i],lwd=1,col="grey")
    if (Q > 0) {
      lines(yrs,meds["50%",],lwd=2,col=4)
      if (Q == 90) {
        lines(yrs,meds["5%",],lwd=1,col=2)
        lines(yrs,meds["95%",],lwd=1,col=2)
      } else {
        lines(yrs,meds["2.5%",],lwd=1,col=2)
        lines(yrs,meds["97.5%",],lwd=1,col=2)
      }
    }
    if (!is.na(hline)) abline(h=hline,lwd=1,col="black",lty=2)
  }
  mtext(paste0(scene,"  ",label),side=2,line=-0.2,outer=TRUE,cex=1.1)
  if (nchar(filen) > 0)
    addplot(filen=filen,rundir=rundir,category="C_vs_MSY",caption)
  return(invisible(outmed))
} # end of plotsceneproj

#' @title plotzonechangerate plot the output from getzonechangerate function
#'
#' @description plotzonechangerate produces a plot of the percent changes within
#'     the whichvar by zone selected across all selected scenarios within
#'     do_comparison. The outcome is placed within the ScenarioPMs tab. The
#'     variables it can work with from zone from the output of do_comparison
#'     include: matureB, exploitB, midyrexpB, catch, acatch, harvestR, cpue,
#'     recruit, deplsB, and depleB.
#'
#' @param rundir the directory in which comparisons are being made. It is best
#'     to be a separate directory form any particular scenario.
#' @param res the output object from the getzonechangerate function
#' @param whichvar which variable within the dyn object to characterize
#' @param glb the globals object needed for the number and names of the sau,
#'     and the years and their names.
#' @param console should each plot go to the console or be saved to rundir.
#'     default = TRUE ie go to console
#'
#' @seealso \link{getzonechangerate}, \link{do_comparison}
#'
#' @returns invisibly the filename used without the path
#' @export
#'
#' @examples
#' # syntax
#' # pickvar <- "acatch"
#' # res <- getzonechangerate(zone=zone,whichvar=pickvar,glb=glb)
#' # filen <- plotrateofchange(rundir=rundir,res=res,whichvar=pickvar,glb=glb)
plotzonechangerate <- function(rundir,res,whichvar,glb,console=TRUE) {
  #  rundir=rundir;res=res;whichvar="catch";glb=glb;console=TRUE
  scenes <- colnames(res$pdiffer)
  nscen <- length(scenes)
  yrs <- glb$pyrnames
  nyrs <- glb$pyrs
  tmp <- res$pdiffer
  miny <- getmin(tmp)
  maxy <- getmax(tmp)
  fname <- paste0(whichvar,"_zonal_rate_of_change_across_scenarios")
  fname1 <- paste0(fname,".png")
  filen <- pathtopath(rundir,fname1)
  if (console) filen=""
  plotprep(width=8,height=4.5,newdev=FALSE,filename=filen,verbose=FALSE)
  parset(margin=c(0.5,0.5,0.1,0.1))

  label <- paste0("Percent Annual Change in zonal ",whichvar)
  plot(yrs,tmp[,1],type="l",lwd=2,col=1,ylab=label,ylim=c(miny,maxy),xlab="",
       panel.first=grid())
  abline(h=0,lwd=1,col=1)
  for (scen in 2:nscen) lines(yrs,tmp[,scen],lwd=2,col=scen)
  legend("topright",scenes,lwd=3,col=1:nscen,bty="n",cex=1.2)
  return(invisible(fname1))
} # plotzonechangerate

#' @title plotzonedevs plots projection deviates for a variable across scenerios
#'
#' @description plotzonedevs provides a plot of a scenario Performance Measure.
#'     In this case the PM are the deviations from a loess fit to whatever input
#'     variable has been selected from the zone objects in do_comparison(). One
#'     needs first generate a list of the selected variable from the the input
#'     zone objects (one from each scenario). A standard comparison would be on
#'     the variation of catches through the years of the projection. But this
#'     could be applied to any of the variables within the zone object, which
#'     includes matureB, exploitB, midyexpB, catch, acatch, harvestR, cpue,
#'     recruit, deplsB, depleB. The final catchN and Nt would need separate
#'     treatment.
#'
#' @param invar the sub-object from each scenario as a list, see examples
#' @param scenes a vector of names fr each scenario being compared
#' @param glb  a globals object from any scenario. This is to obtain the
#'     number of projection years, which MUST be the same across
#'     all scenarios (else the matrices will not line up).
#' @param filen the filename and path if the plot is to be saved, default = "".
#'     if a filename is given it should be a .png file
#'
#' @return a list of a matrix of the mean and sd devs per scenario, and a matrix
#'      of the actual deviates, all returned invisibly. It also plots a graph
#' @export
#'
#' @examples
#' \dontrun{
#' # example syntax, assumes zone has been extracted from all scenarios
#'    glb <- glbc[[1]]
#'    nscen <- length(scenes)
#'    whichyrs <- (glb$hyrs + 1):(glb$hyrs + glb$pyrs)
#'    x <- makelist(scenes)
#'    for (scen in 1:nscen) x[[scen]] <- zone[[scen]]$catch[whichyrs,]
#'    devout <- plotzonedevs(x,scenes,glb,filen="")
#' }
plotzonedevs <- function(invar,scenes,glb,filen=""){
  # invar=x;scenes=scenes;filen=""
  nscen <- length(scenes)
  yrs <- glb$pyrnames
  nyrs <- length(yrs)
  reps <- glb$reps
  summarydevs <- matrix(0,nrow=nscen,ncol=2,
                        dimnames=list(scenes,c("Mean","SD")))
  devs <- matrix(0,nrow=reps,ncol=nscen,dimnames=list(1:reps,scenes))
  for (scen in 1:nscen) {
    # scen = 1
    scenvar <- invar[[scen]]
    devs[,scen] <- apply(scenvar,2,getabdevs,years=yrs,smooth=TRUE)
    summarydevs[scen,"Mean"] <- mean(devs[,scen],na.rm=TRUE)
    summarydevs[scen,"SD"] <- sd(devs[,scen],na.rm=TRUE)
  }
  left <- getmin(devs)
  right <- getmax(devs)
  plotprep(width=9,height=6,newdev=FALSE,filename=filen,verbose=FALSE)
  parset(plots=pickbound(nscen),margin=c(0.5,0.5,0.1,0.1),byrow=TRUE)
  for (scen in 1:nscen) {
    hist(devs[,scen],breaks=20,main="",xlab="",ylab=scenes[scen],
         xlim=c(left,right))
    abline(v=summarydevs[scen,"Mean"],lwd=2,col=4)
    label <- paste0(scenes[scen],"-",round(summarydevs[scen,"Mean"],2))
    mtext(label,side=3,line=-1.1,cex=1.2)
  }
  return(invisible(list(summarydevs=summarydevs,devs=devs)))
} # end of plotzonedevs

#' @title plotzonedyn ribbonplot of zone catch, mature biomass, harvestR and CE
#'
#' @description plotzonedyn generates ribbonplots of the zone-wide total catches,
#'     the mature biomass, the harvest rate, and the CPUE across the years of
#'     projection.
#'
#' @param rundir the directory in which all results are held for a scenario or
#'     comparison of scenarios
#' @param scenes the names of the different scenarios being compared
#' @param zone a list of the out$outzone objects from each chosen scenario
#' @param glbc all of the global objects from out. They should all have the
#'     same reps and years but this will be checked.
#' @param console should the plot be sent to the console or saved for use in the
#'     web-page output? default=TRUE, set to FALSE to save it
#' @param q90 should the 90th quantiles be plotted? default=TRUE. If FALSE then
#'     the 95th quantiles are used.
#' @param polys should transparent polygons be plotted = TRUE, or lines = FALSE
#' @param intens if polys=TRUE then intens signifies the intensity of colour on
#'     a scale of 0 - 255. 127 is about 50 percent dense.
#' @param hlines should reference lines be drawn on each plot. default=NULL
#'     which = no lines if lines are required then hlines needs to be set to a
#'     list of h values for each scenario, for each plotted variable, see
#'     examples for syntax
#' @param scencol what sequence of colours should each scenario have. The default=NULL
#'     which implies the colour will reflect the sequence in which they are
#'     read in.
#'
#' @seealso {
#'   \link{poptozone}, \link{plotZone}
#' }
#'
#' @return nothing but it does generate a plot and saves it if console=FALSE
#' @export
#'
#' @examples
#' print("wait on data sets")
#' # syntax plotzonedyn(rundir,scenes,zone,glbc,console=FALSE,q90=TRUE,
#' #                    polyd=TRUE,intens=100,hlines=list(catch=outprod[,"MSY"],
#' #                    spawnB=outprod[,"Bmsy"],harvestR=0,
#' #                    cpue=outprod[,"CEmsy"]))
plotzonedyn <- function(rundir,scenes,zone,glbc,console=TRUE,q90=TRUE,
                        polys=TRUE,intens=127,hlines=NULL,scencol=NULL) {
  # rundir=rundir;scenes=scenes;zone=zone;glbc=glbc;console=TRUE;q90=TRUE;
  # polys=TRUE;intens=127; hlines=list(catch=outprod[,"MSY"],harvestR=0,
  #  spawnB=outprod[,"Bmsy"],cpue=outprod[,"CEmsy"], expB=outprod[,"Bexmsy"])
  nscen <- length(scenes)
  allreps <- unlist(lapply(glbc,"[[","reps"))
  if (min(allreps) != max(allreps)) {
    label <- paste0("WARNING. Cannot compare scenarios with different replicate",
                    " numbers so some plots of dynamics omitted \n")
    cat(label)
  } else {
    reps <- max(allreps)
    glb <- glbc[[1]]
    hyrs <- glb$hyrs
    pyrs <- glb$pyrs
    yrs <- hyrs + pyrs
    pyrnames <- glb$pyrnames
    yrnames <- as.numeric(c(tail(glb$hyrnames,1),pyrnames))
    varx <- array(0,dim=c((pyrs+1),reps,nscen),
                  dimnames=list(yrnames,1:reps,scenes))
    qlab <- c("2.5%","5%","50%","95%","97.5%")
    varq <- array(0,dim=c(5,(pyrs+1),nscen),
                  dimnames=list(qlab,yrnames,scenes))
    filen=""
    if (!console) {
      filen <- filenametopath(rundir,"zonedynamics.png")
      caption <- paste0("Outputs for each scenario summarized across the whole ",
                     " zone. Plots of Catch, MatureB, HarvestR, exploitB, and ",
                  "Recruits. Horizintal lines are MSY, Bmsy, CEmsy, and eBmsy.")
    }
    plotprep(width=8,height=9,newdev=FALSE,filename=filen,verbose=FALSE)
    parset(plots=c(3,2),margin=c(0.3,0.4,0.05,0.05),outmargin=c(0,1,0,0),
           byrow=TRUE)
    for (i in 1:nscen) {
      varx[,,i] <- zone[[i]]$catch[hyrs:yrs,]
      varq[,,i] <- apply(varx[,,i],1,quantile,probs=c(0.025,0.05,0.5,0.95,0.975),
                         na.rm=TRUE)
    }
    projcatch <- varx  # save for phase plot
    doquantplot(varq,varname="Catch (t)",yrnames,scenes,q90=q90,polys=polys,
                intens=intens,hlin=hlines[[1]],scencol=scencol)
    for (i in 1:nscen) {
      varx[,,i] <- zone[[i]]$matureB[hyrs:yrs,]
      varq[,,i] <- apply(varx[,,i],1,quantile,probs=c(0.025,0.05,0.5,0.95,0.975),
                         na.rm=TRUE)
    }
    projmatB <- varx  # save for phase plot
    doquantplot(varq,varname="Mature Biomass",yrnames,scenes,q90=q90,polys=polys,
                intens=intens,scencol=scencol)
    for (i in 1:nscen) {
      varx[,,i] <- zone[[i]]$harvestR[hyrs:yrs,]
      varq[,,i] <- apply(varx[,,i],1,quantile,probs=c(0.025,0.05,0.5,0.95,0.975),
                         na.rm=TRUE)
    }
    doquantplot(varq,varname="Harvest Rate",yrnames,scenes,q90=q90,polys=polys,
                intens=intens,scencol=scencol)
    for (i in 1:nscen) {
      varx[,,i] <- zone[[i]]$cpue[hyrs:yrs,]
      varq[,,i] <- apply(varx[,,i],1,quantile,probs=c(0.025,0.05,0.5,0.95,0.975),
                         na.rm=TRUE)
    }
    doquantplot(varq,varname="CPUE",yrnames,scenes,q90=q90,polys=polys,
                intens=intens,scencol=scencol)
    for (i in 1:nscen) {
      varx[,,i] <- zone[[i]]$exploitB[hyrs:yrs,]
      varq[,,i] <- apply(varx[,,i],1,quantile,probs=c(0.025,0.05,0.5,0.95,0.975),
                         na.rm=TRUE)
    }
    projexpB <- varx
    doquantplot(varq,varname="Exploitable Biomass",yrnames,scenes,q90=q90,
                polys=polys,intens=intens,scencol=scencol)
    for (i in 1:nscen) {
      varx[,,i] <- zone[[i]]$recruit[hyrs:yrs,]
      varq[,,i] <- apply(varx[,,i],1,quantile,probs=c(0.025,0.05,0.5,0.95,0.975),
                         na.rm=TRUE)
    }
    doquantplot(varq,varname="Recruitment",yrnames,scenes,q90=q90,
                polys=polys,intens=intens,scencol=scencol)
    label <- paste0(scenes,collapse=",")
    legend("bottomright",legend=scenes,lwd=3,col=1:nscen,bty="n",cex=1.1)
    mtext("Zone Wide Dynamics",side=2,line=-0.2,outer=TRUE,cex=1.1)
    if (!console) {
      addplot(filen=filen,rundir=rundir,category="zone",caption)
    }
    # phase plot of catchdivMSY vs MatB/Bmsy
    catbymsy <- projcatch
    matBbyBmsy <- projmatB
    for (i in 1:nscen) {
      catbymsy[,,i] <- projcatch[,,i]/hlines$catch[i]
      matBbyBmsy[,,i] <- projmatB[,,i]/hlines$spawnB[i]
    }
    yrs <- c(tail(glb$hyrnames,1),glb$pyrnames)
    nyrs <- length(yrs)
    yvar <- xvar <- array(0,c(nyrs,nscen),
                          dimnames=list(yrs,scenes))
    for (scen in 1:nscen) {
      xvar[,scen] <- apply(matBbyBmsy[,,scen],1,median)
      yvar[,scen] <- apply(catbymsy[,,scen],1,median)
    }
    xlabel="Median Mature Biomass div Bmsy"
    ylabel="Median Actual Catches div MSY"
    simplephaseplot(rundir=rundir,xvar=xvar,yvar=yvar,xlabel=xlabel,ylabel=ylabel,
                    console=console,width=6,height=6,hline=1.0,vline=1.0,
                    scenes=scenes)
    # phase plot of catchdivMSY vs matBdepletion
    pickyrs <- glb$hyrs:(glb$hyrs+glb$pyrs)
    for (scen in 1:nscen) {
      xvar[,scen] <- apply(zone[[scen]]$deplsB[pickyrs,],1,median)
      yvar[,scen] <- apply(catbymsy[,,scen],1,median)
    }
    xlabel="Median Mature Biomass Depletion"
    ylabel="Median Actual Catches div MSY"
    simplephaseplot(rundir=rundir,xvar=xvar,yvar=yvar,xlabel=xlabel,ylabel=ylabel,
                    console=console,width=6,height=6,hline=1.0,vline=1.0,
                    scenes=scenes)
    # phase plot of catchdivMSY vs expB/Bexmsy
    catbymsy <- projcatch
    expBbyBexmsy <- projexpB
    for (i in 1:nscen) {
      catbymsy[,,i] <- projcatch[,,i]/hlines$catch[i]
      expBbyBexmsy[,,i] <- projmatB[,,i]/hlines$expB[i]
    }
    yrs <- c(tail(glb$hyrnames,1),glb$pyrnames)
    nyrs <- length(yrs)
    yvar <- xvar <- array(0,c(nyrs,nscen),
                          dimnames=list(yrs,scenes))
    for (scen in 1:nscen) {
      xvar[,scen] <- apply(expBbyBexmsy[,,scen],1,median)
      yvar[,scen] <- apply(catbymsy[,,scen],1,median)
    }
    xlabel="Median Exploitable Biomass div Bexmsy"
    ylabel="Median Actual Catches div MSY"
    simplephaseplot(rundir=rundir,xvar=xvar,yvar=yvar,xlabel=xlabel,ylabel=ylabel,
                    console=FALSE,width=6,height=6,hline=1.0,vline=1.0,
                    scenes=scenes)
    # phase plot of catchdivMSY vs expBdepl
    pickyrs <- glb$hyrs:(glb$hyrs+glb$pyrs)
    for (scen in 1:nscen) {
      xvar[,scen] <- apply(zone[[scen]]$depleB[pickyrs,],1,median)
      yvar[,scen] <- apply(catbymsy[,,scen],1,median)
    }
    xlabel="Median Exploitable Biomass Depletion"
    ylabel="Median Actual Catches div MSY"
    simplephaseplot(rundir=rundir,xvar=xvar,yvar=yvar,xlabel=xlabel,ylabel=ylabel,
                    console=console,width=6,height=6,hline=1.0,vline=1.0,
                    scenes=scenes)
  } # end of else statement
} # end of plotzonedyn

#' @title plotzoneyrtovar plots the time to invar at a zone scale
#'
#' @description plotzoneyrtovar it is common to want to know how long each
#'     replicate in a scenario run takes to be >= to given value, eg when is the
#'     catch/msy ratio first = 1. The zonecatchvsmsy function generates that
#'     catch/msy information, this plots the histograms. Other varnames such
#'     as matureB would be of interest when compared to Bmsy.
#'
#' @param rundir the directory in which all results are held for a scenario or
#'     comparison of scenarios
#' @param invar a list of the zone-scale replicate runs for the varname from
#'     each scenario. the zone object is generated inside do_comparison
#' @param varname just the name of the variable being plotted, to make sure the
#'     figures all have the correct labelling
#' @param glbc the list of the global objects generated within do_comparison
#' @param category if a file is saved this argument identifies the webpage tab
#'     into which the plot will be placed.
#' @param console should the plot be sent to the console or saved for use in the
#'     web-page output? default=TRUE, set to FALSE to save it
#'
#' @seealso{
#'    \link{whenfirst},  \link{zonecatchvsmsy}, \link{do_comparison}
#' }
#'
#' @returns nothing but it does generate a plot of nscen histograms
#' @export
#'
#' @examples
#' print("wait on example data")
#' # syntax  plotzoneyrtovar(rundir=rundir,invar=zyrtomsy,varname="MSY",
#' #                         glbc=glbc,category="C_vs_MSY",console=TRUE)
plotzoneyrtovar <- function(rundir,invar,varname,glbc,category,console=TRUE) {
  #  rundir=rundir; invar=zyrtomsy;varname="MSY";glbc=glbc;category="C_vs_MSY";console=TRUE
  scenes <- colnames(invar)
  nscen <- length(scenes)
  nyrs <- glbc[[1]]$pyrs
  reps <- glbc[[1]]$reps
  filen=""
  if (!console) {
    filen <- pathtopath(rundir,paste0("Zone_scale_years_to_",varname,".png"))
    caption <- paste0("Histograms of years to ",varname," for Zones in ",
                      paste0(scenes," ",collapse=""))
  }
  numNA <- apply(invar,2,countNAs)

  plotprep(width=9,height=6,newdev=FALSE,filename=filen,verbose=FALSE)
  parset(plots=pickbound(nscen),cex=1.0,margin=c(0.35,0.5,0.2,0.05),
         outmargin=c(1,1,0,0))
  for (scen in 1:nscen) {  # scen = 1
    if (numNA[scen] == reps) {
      plotnull(msg="None at MSY",cex=1.2)
    } else {
      histdat <- invar[,scen]
      inthist(x=histdat,col="lightblue",border=1,width=0.8,xmin=1,xmax=nyrs,
              ylab=scenes[scen],panel.first=grid())
      mtext(paste0("NAs = ",numNA[scen]),side=3,outer=FALSE,cex=1.1,line=0)
    }
  }
  mtext(paste0("Years to First Reach ",varname),side=1,outer=TRUE,cex=1.2,
        line=-0.2)
  mtext("Count of each Interval",side=2,outer=TRUE,cex=1.2,line=-0.2)
  if (!console) {
    addplot(filen=filen,rundir=rundir,category=category,caption=caption)
  }
} # end of plotzoneyrtovar

#' @title sauquantbyscene calcualtes the quantiles for all sau and all scenarios
#'
#' @description sauquantbyscene takes the output from the scenebyvar function
#'     and for each sau extracts the quantiles for the input variable for
#'     each scenario. The output is a list nsau long each containing the five
#'     quantiles for the variable input for each scenario ready for use by
#'     doquantplot, or whatever other analysis is required.
#'
#' @param invar the output from the scenebyvar function which should be an array
#'     of nyrs x nsau x replicates for the variable selected, be it catch, cpue,
#'     matureB, etc.
#' @param glb the globals object. The comparisons invovled imply that all
#'     scenarios considered should have the same number of years and replicates.
#'
#' @seealso{
#'     \link{scenebyvar}, \link{doquantplot}, \link{do_comp_outputs}
#' }
#'
#' @return a list of length nsau containing the quantiles for each scenario in
#'     each sau, for the input variable. The quantiles include 0.025, 0.05, 0.5,
#'     0.95, and 0.975.
#' @export
#'
#' @examples
#' print("wait on data sets")
#' # invar=catch; glb=glbc
sauquantbyscene <- function(invar,glb) {
  # nyrs <- dim(invar[[1]])[1]
  # yrnames <- as.numeric(unlist(dimnames(invar[[1]])[1]))
  nscen <- length(invar)
  scenes <- names(invar)
  qs <- c("2.5%","5%","50%","95%","97.5%")
#  scenquant <- array(0,dim=c(5,nyrs,nscen),dimnames=list(qs,yrnames,scenes))
  scenquant <- makelist(scenes)
  for (i in 1:nscen) {  # i = 1
    scen <- invar[[i]]
    nsau <- glb[[i]]$nSAU
    saunames <- glb[[i]]$saunames
    sauquant <- makelist(saunames)
    for (sau in 1:nsau) {
      sauquant[[sau]] <- apply(scen[,sau,],1,quants)
    }
    scenquant[[i]] <- sauquant
  }
  return(scenquant)
} # end of sauquantbyscene


#' @title sauribbon a ribbon plot the quantiles x scenes for the input variable
#'
#' @description sauribbon is used to produce a ribbon plot for a given sau across
#'     a set of scenarios for a given variable (cpue, catch, harvest rate, etc)
#'     out of the dynamics object. After running do_comparison
#'
#' @param rundir the directory in which all results are held for a scenario or
#'     comparison of scenarios
#' @param scenes the names of the different scenarios being compared
#' @param varqnts the quantiles for each scenario as output by sauquantbyscene
#' @param glbc the list of the global objects. They should all have the same
#'     number of sau.
#' @param varname just the name of the variable being plotted, to make sure the
#'     figures all have the correct labelling
#' @param console should the plot be sent to the console or saved for use in the
#'     web-page output? default=TRUE, set to FALSE to save it
#' @param q90 should the 90th quantiles be plotted? default=TRUE. If FALSE then
#'     the 95th quantiles are used.
#' @param intens if polys=TRUE then intens signifies the intensity of colour on
#'     a scale of 0 - 255. 127 is about 50 percent dense.
#' @param addleg add a legend or not. If not use "", if yes use any of topleft,
#'     topright, bottomleft, or bottomright, default = bottomright.
#' @param addmedian should each polygon also have its median plotted as a line.
#' @param defpar should the formatting of the plot use the internal version
#'     (the default setting) or, if set = FALSE use an external definition.
#'
#' @seealso {
#'   \link{sauquantbyscene}, \link{doquantribbon}, \link{do_comp_outputs}
#' }
#'
#' @return nothing but it generates a plot and saves it if console=FALSE
#' @export
#'
#' @examples
#' \dontrun{
#'  # a possible example run where scenarios 1, 13, and 14 are compared.
#'    result <- do_comparison(rundir=rundir,postfixdir=postfixdir,outdir=outdir,
#'                            files=files,pickfiles=c(1,2,3),verbose=TRUE)
#'    out <- do_comp_outputs(result,projonly=TRUE)
#'    catch <- scenebyvar(dyn=out$dyn,byvar="catch",glbc=out$glbc)
#'    sauribbon("",scenes=out$scenes,sau=8,varqnts=catqnts,glb=out$glbc,
#'    varname="Catch",console=TRUE,q90=TRUE,intens=100,addleg="bottomright")
#' }
sauribbon <- function(rundir,scenes,varqnts,glbc,varname,
                      console=TRUE,q90=TRUE,intens=127,addleg="bottomright",
                      addmedian=3,defpar=TRUE) {
  # rundir=rundir;scenes=scenes;varqnts=catqnts;glb=glbc;varname="catch";
  # console=TRUE;q90=Q90;intens=intensity;addleg=ribbonleg;defpar=TRUE;addmedian=TRUE
  filen <- ""
  nscen <- length(scenes)
  numsau <- sapply(glbc,"[[","nSAU")
  maxsau <- max(numsau)
  if (all(numsau == maxsau)) {
    for (sau in 1:maxsau) { # sau = 2
      sauname <- glbc[[1]]$saunames[sau]
      if (!console) {
        filename <- paste0(sauname,"_",varname,"_ribbon.png")
        filen <- pathtopath(rundir,filename)
        caption <- paste0(varname," ribbon plot for ",glbc[[1]]$saunames[sau])
      }
      if (defpar) {
        oldpar <- par(no.readonly=TRUE)
        plotprep(width=8,height=4,newdev=TRUE,filename=filen,verbose=FALSE)
        parset(cex=1.1,margin=c(0.35,0.5,0.05,0.05))
      }
      rc <- dim(varqnts[[1]][[1]])
      rcnames <- dimnames(varqnts[[1]][[1]])
      rcnames[[3]] <- scenes
      qntvar <- array(0,dim=(c(rc,nscen)),dimnames=rcnames)
      for (i in 1:nscen) qntvar[,,i] <- varqnts[[i]][[sau]]
      doquantribbon(qntvar,varname=varname,sauname=sauname,yrnames=glbc[[1]]$pyrnames,
                    scenes=scenes,q90=q90,intens=intens,addleg=addleg,
                    addmedian=addmedian)
      if (!console) {
        addplot(filen=filen,rundir=rundir,category=varname,caption=caption)
      }
    }
  } else {
    warning("Differing number of SAU, not ribbon plts produced. \n")
  }
  if (defpar) return(invisible(oldpar))
} # end of sauribbon

#' @title scenarioproperties tabulates important properties of each scenario
#'
#' @description scenarioproperties is used with the do_comparison function to
#'     tabulate a set of important properties of each scenario so that the
#'     comparability of each scenario can be quickly assessed.
#'
#' @param scenes the runnames of each scenario
#' @param glbc the globals object or list for each scenario
#' @param ctrlc the ctrl object for each scenario
#' @param condCc the condC object for each scenario
#'
#' @seealso{
#'   \link{do_comparison}
#' }
#'
#' @return a list containing a matrix of a set of scenario properties to check
#'     comparability, and a vector identifying whether all columns are identical
#'     = 1 or some are different = 0
#' @export
#'
#' @examples
#' print("wait on data sets")
scenarioproperties <- function(scenes,glbc,ctrlc,condCc) {
  rows <- c("reps","randseed","sigR","sigB","sigCE","projyrs",
            "numpop","nSAU","Nclass","hyrs","pyrs","larvdisp","initdepl")
  numrow <- length(rows)
  numcol <- length(scenes)
  ansout <- matrix(0,nrow=numrow,ncol=numcol,dimnames=list(rows,scenes))
  for (i in 1:numcol) {
    initdepl <- min(condCc[[i]]$initdepl,na.rm=TRUE)
    ansout[,i] <- c(ctrlc[[i]]$reps,ctrlc[[i]]$randseed,ctrlc[[i]]$withsigR,
                    ctrlc[[i]]$withsigB,ctrlc[[i]]$withsigCE,
                    ctrlc[[i]]$projection,glbc[[i]]$numpop,glbc[[i]]$nSAU,
                    glbc[[i]]$Nclass,glbc[[i]]$hyrs,glbc[[i]]$pyrs,
                    glbc[[i]]$larvdisp,initdepl)
  }
  same <- numeric(numrow)
  for (i in 1:numrow) { # ansout <- scenprops; i=1
    same[i] <- ifelse(min(ansout[i,]) == max(ansout[i,]),1,0)
  }
  ansout <- cbind(ansout,same)
  return(ansout)
} # end of scenarioproperties

#' @title scenebyvar extracts as sau x variable x scenario from do_comparison output
#'
#' @description scenebyvar facilitates making comparisons between scenarios by
#'     at the sau level by extracting individual variables from the 'dyn'
#'     object that contains each scenarios dynamics by sau and is part of the
#'     output from the do_comparison function.
#'
#' @param dyn a list of the dynamics at an sau scale for each scenario in the
#'     comparison made within do_comparison
#' @param byvar the variable to be extracted from the dyn object. Valid values
#'     are matureB, exploitB, midyexpB, catch, acatch, harvestR, cpue, recruit,
#'     deplsB, and depleB. Anything will stop the function and shed an error
#'     message.
#' @param glb the globals object, which shuold be the same across scenarios
#' @param projonly should only the projection years be output, default=TRUE
#'
#' @seealso{
#'    \link{do_comparison}
#' }
#'
#' @return a list of the nscene projyrs x nsau x reps for the byvar values
#' @export
#'
#' @examples
#' print("wait on example data sets")
scenebyvar <- function(dyn,byvar,glb,projonly=TRUE) {
  knownvar <- c("matureB","exploitB","midyexpB","catch","acatch","harvestR",
                "cpue","recruit","deplsB","depleB")
  if (byvar %in% knownvar) {
    nscen <- length(dyn)
    scenes <- names(dyn)
    varout <- makelist(scenes)
    if (projonly) {
      for (i in 1:nscen) varout[[i]] <- projectiononly(dyn[[i]][[byvar]],glb[[i]])
    } else {
      for (i in 1:nscen) varout[[i]] <- dyn[[i]][[byvar]]
    }
  } else {
    cat("Unknown variable name in [scenebyvar] function \n")
    cat("Allowed terms: ",knownvar,"\n")
    stop()
  }
  return(invisible(varout))
} # end of scenbyvar

#' @title scenebyzone extracts as variable x scenario from do_comparison output
#'
#' @description scenebyzone facilitates making comparisons between scenarios
#'     at the zone level by extracting individual variables from the 'zone'
#'     object that contains each scenarios dynamics at the zone scale, where
#'     each variable has been amalgamated across populations and sau in an
#'     appropriate manner. cpue, for example, is catch weighted, whereas matureB
#'     is simply summed. The zone object is part of the
#'     output from the do_comparison function.
#'
#' @param zone a list of the dynamics at a zone scale for each scenario in the
#'     comparison made within do_comparison
#' @param byvar the variable to be extracted from the dyn object. Valid values
#'     are matureB, exploitB, midyexpB, catch, acatch, harvestR, cpue, recruit,
#'     deplsB, and depleB. Anything will stop the function and shed an error
#'     message.
#' @param glb the globals object, which should be the same across scenarios
#' @param projonly should only the projection years be output, default=TRUE
#'
#' @seealso{
#'    \link{do_comparison}
#' }
#'
#' @return a list of the nscene projyrs x reps for the byvar values
#' @export
#'
#' @examples
#' print("wait on example data sets")
scenebyzone <- function(zone,byvar,glb,projonly=TRUE) {
  knownvar <- c("matureB","exploitB","midyexpB","catch","acatch","TAC",
                "harvestR","cpue","recruit","deplsB","depleB")
  if (byvar %in% knownvar) {
    nscen <- length(zone)
    scenes <- names(zone)
    varout <- makelist(scenes)
    if (projonly) {
      projyr <- (glb$hyrs + 1):(glb$hyrs + glb$pyrs)
      for (i in 1:nscen) varout[[i]] <- zone[[i]][[byvar]][projyr,]
    } else {
      for (i in 1:nscen) varout[[i]] <- zone[[i]][[byvar]]
    }
  } else {
    cat("Unknown variable name in [scenebyvar] function \n")
    cat("Allowed terms: ",knownvar,"\n")
    stop()
  }
  return(invisible(varout))
} # end of scenbyzone

#' @title state_by_year applies a selected function to replicates in dynamics
#'
#' @description state_by_year applies a user selected or generated function to
#'     each replicate in each SAU within one of the variables inside the
#'     dynamics array. This can be used to determine the mean, median, stdev,
#'     or any other custom function one cares to define. In the dyn object from
#'     do_comparison the valid whichvar values are: matureB, exploitB,
#'     midyexpB, catch, acatch, harvestR, cpue, recruit, deplsB, and depleB.
#'     The year numbers refer to the total number of years in the simulation,
#'     both the historic catches and the projection years.
#'
#'
#' @param indyn usually the dynamics object from do_comparison
#' @param whichvar which dynamics variable to use as a character string
#' @param whichyrs which year numbers to use
#' @param whichfun any R function applicable to vectors or numbers or a custom
#'     user-defined function
#' @param glb the globals object needed for the number and names of the sau,
#'     and the years and their names.
#'
#' @return array of function values for each of the years, scenarios, and sau
#' @export
#'
#' @examples
#' \dontrun{
#' # assumes that do_comparison function has been run.
#'   res <- state_by_year(indyn=dyn,whichvar="catch",whichyrs=c(63,68),
#'                        whichfun=mean,glb=glb)
#' # Using the first 4 examples from aMSEGuide should give something like
#'   round(res[1,,],3)
#'            sau6   sau7   sau8   sau9  sau10  sau11   sau12  sau13
#'   EG      7.208 17.537  8.811 41.957 37.559  89.514 48.686 33.714
#'   EGMR12  7.015 18.477  9.357 43.826 40.893  90.421 52.538 34.061
#'   EGMR3   8.754 18.774 11.639 47.301 35.875  94.545 51.320 36.831
#'   EGMRall 7.962 19.887 12.278 53.922 41.175 101.453 57.521 35.796
#' }
state_by_year <- function(indyn,whichvar,whichyrs,whichfun,glb) {
  # indyn=dyn; whichvar="catch"; whichyrs=c(63,68); whichfun=mean; glb=glb
  nsau <- glb$nSAU
  saunames <- glb$saunames
  years <- c(glb$hyrnames,glb$pyrnames)
  nyrs <- length(whichyrs)
  yrnames <- years[whichyrs]
  scenes <- names(indyn)
  nscen <- length(scenes)
  indexvar <- match(whichvar,names(indyn[[1]]))
  res <- array(0,dim=c(nyrs,nscen,nsau),dimnames=list(yrnames,scenes,saunames))
  for (scen in 1:nscen) {# scen = 1
    pickvar <- indyn[[scen]][[indexvar]]
    for (yr in 1:nyrs) res[yr,scen,] <- apply(pickvar[whichyrs[yr],,],1,whichfun)
  }
  return(res)
} # end of state-by-year

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

#' @title tabulateproductivity produces tables of the MSY and related statistics
#'
#' @description tabulateproductivity writes out csv files for each scenario
#'     that contain the B0, Bmsy, MSY, Depmsy, and CEmsy for each sau
#'
#' @param rundir the full path to the directory holding the all result files
#' @param prods the sauprod objects from each scenario in a list
#' @param scenes a vector of the names for each scenario
#'
#' @return It adds a csv file to rundir (or nscenario csv files if unexpectedly
#'     the productivity tables differ) and returns either a single identical
#'     matrix of productivity stats, or the full list.
#' @export
#'
#' @examples
#' print("wait on data sets")
tabulateproductivity <- function(rundir,prods,scenes) {
  nscen <- length(scenes)
  dimsame <- matrix(0,nrow=(nscen-1),ncol=2)
  colnames(dimsame) <- c("rows","cols")
  for (i in 2:nscen) { # i = 2
    dimsame[i-1,"rows"] <- abs(dim(prods[[i-1]])[1] - dim(prods[[i]])[1])
    dimsame[i-1,"cols"] <- abs(dim(prods[[i-1]])[2] - dim(prods[[i]])[2])
  }
  if (!all(dimsame < 1)) {
    answer <- paste0("Some productivity matrices have different dimensions.",
                     " most likely that some scenarios were run using an ",
                     "old version of aMSE \n")
    warning(cat(answer))
    for (i in 1:nscen) {
      filen <- paste0("productivity_",scenes[i],".csv")
      caption <- paste0("SAU productivity statistics for ",scenes[i],".",
                        " Productivity matrices NOT identical.")
      addtable(prods[[i]],filen,rundir=rundir,category="productivity",caption)
      out <- prods
    }
  } else {
    same <- vector(mode="logical",length=(nscen-1))
    count <- 1
    for (i in 2:nscen) same[i-1] <- all(prods[[i-1]] == prods[[i]])
    if (all(same)) {
      filen <- paste0("productivity_",scenes[i],".csv")
      caption <- paste0("SAU productivity statistics for ",scenes[i],".",
                        " All productivity matrices identical.")
      addtable(prods[[i]],filen,rundir=rundir,category="productivity",caption)
      out <- prods[[1]]
    } else {
      for (i in 1:nscen) {
        filen <- paste0("productivity_",scenes[i],".csv")
        caption <- paste0("SAU productivity statistics for ",scenes[i],".",
                          " Productivity matrices NOT identical.")
        addtable(prods[[i]],filen,rundir=rundir,category="productivity",caption)
        out <- prods
      }
    }
  }
  return(invisible(out))
} # end of tabulateproductivity

#' @title tabulateproductivity produces tables of the MSY and related statistics
#'
#' @description tabulateproductivity writes out csv files for each scenario
#'     that contain the B0, Bmsy, MSY, Depmsy, and CEmsy for each sau
#'
#' @param rundir the full path to the directory holding the all result files
#' @param prods the sauprod objects from each scenario in a list
#' @param scenes a vector of the names for each scenario
#'
#' @return It adds a csv file to rundir (or nscenario csv files if unexpectedly
#'     the productivity tables differ) and returns either a single identical
#'     matrix of productivity stats, or the full list.
#' @export
#'
#' @examples
#' print("wait on data sets")
#' # syntax tabulatezoneprod(rundir=rundir,prods=prods,scenes=scenes)
tabulatezoneprod <- function(rundir,prods,scenes) {
  nscen <- length(scenes)
  columns <- c("B0","Bmsy","MSY","Dmsy","CEmsy","Bexmsy")
  out <- matrix(0,nrow=nscen,ncol=length(columns),dimnames=list(scenes,columns))
  for (i in 1:nscen) { # i=2
    scenprod <- prods[[i]]
    out[i,1:3] <- colSums(scenprod[,1:3])
    out[i,"Dmsy"] <- mean(scenprod[,"Dmsy"])
    out[i,"CEmsy"] <- sum((scenprod[,3]/out[i,3]) * scenprod[,5])
    out[i,"Bexmsy"] <- sum(scenprod[,"Bexmsy"])
  }
  filen <- paste0("zone_productivity.csv")
  caption <- paste0("Productivity statistics for the whole zone.")
  addtable(out,filen,rundir=rundir,category="zone",caption)
  return(invisible(out))
} # end of tabulatezoneprod

#' @title whenfirst returns the first occurrence of a value in a vector
#'
#' @description whenfirst is a utility function to be used primarily within an
#'     apply function call. It could then search across all replicates for the
#'     first occurrence of a value >= a given target value.
#'
#' @param x a vector of yearly values within a single replicate
#' @param value the value to be searched for within x as in x >= value
#'
#' @returns a vector replicates long of the first year in which the value occurs
#' @export
#'
#' @examples
#' x <- seq(1,100,1)
#' whenfirst(x,75)
whenfirst <- function(x,value) {
  return(which(x >= value)[1])
}

#' @title zonecatchvsmsy calc catch/msy and year-to-msy by zone by scenario
#'
#' @description zonecatchvsmsy uses the prods, zone, and glbc lists to
#'     estimate the ratio of catch/MSY for each replicate of the zone catches,
#'     plus it determines the first year in each replicate in which MSY is
#'     achieved, If MSY is not met in a replicate an NA is returned. To compare
#'     scenarios they must all have the same number of replicate projections
#'     and must all have the same number of projection years
#'
#' @param prods a list of the production statistics matrix for each scenario,
#'     this is the out$sauprod matrix stored for each scenario.
#' @param zone the list of zone outputs for each scenario
#' @param glbc the globals object for each scenario (needed in case the hyrs
#'     differ between scenarios)
#'
#' @return an invisible list consisting of a list of the catch/msy by zone for
#'     each scenario and a matrix of the year to first meeting the MSY in each
#'     replicate zone catch projection in each scenario
#' @export
#'
#' @examples
#' print("wait on data sets")
#' # syntax:  zonecatchvsmsy(prods,zone,glbc)
zonecatchvsmsy <- function(prods,zone,glbc) {
  scenes <- names(prods)
  nscen <- length(scenes)
  columns <- c("B0","Bmsy","MSY","Bexmsy")
  zprod <- matrix(0,nrow=nscen,ncol=length(columns),
                  dimnames=list(scenes,columns))
  zcatch <- makelist(scenes)
  glb <- glbc[[1]]
  yrs <- glb$pyrnames
  nyrs <- glb$pyrs    # all scenarios need same number of projection years
  reps <- glb$reps    # and same number of replicates
  for (scen in 1:nscen) { #  scen = 1
    hyrs <- glbc[[scen]]$hyrs
    pickprod <- prods[[scen]]
    zprod[scen,] <- colSums(pickprod[,columns])
    zcatch[[scen]] <- zone[[scen]]$catch[(hyrs+1):(hyrs+nyrs),]
  }
  zcvsmsy <- makelist(scenes)
  zyrtomsy <- matrix(0,nrow=reps,ncol=nscen,dimnames=list(1:reps,scenes))
  for (scen in 1:nscen) { # scen = 1
    catches <- zcatch[[scen]]
    msy <- zprod[scen,"MSY"]
    zcvsmsy[[scen]] <- catches/msy
    zyrtomsy[,scen] <- apply(catches,2,whenfirst,value=msy)
  }
  return(invisible(list(zcvsmsy=zcvsmsy,zyrtomsy=zyrtomsy)))
} # end of zonecatchvsmsy

#' @title zoneribbon a ribbon plot the quantiles x scenes for the input variable
#'
#' @description zoneribbon is used to produce a ribbon plot for the zone across
#'     a set of scenarios for a given variable (cpue, catch, harvest rate, etc)
#'     out of the zone object. After running do_comparison
#'
#' @param rundir the directory in which all results are held for a scenario or
#'     comparison of scenarios
#' @param scenes the names of the different scenarios being compared
#' @param invar a list of the zone-scale replicate runs for the varname from
#'     each scenario. the zone object is generated inside docomparison
#' @param glbc the list of the global objects.
#' @param varname just the name of the variable being plotted, to make sure the
#'     figures all have the correct labelling
#' @param category if a file is saved this argument identifies the webpage tab
#'     into which the plot will be placed.
#' @param console should the plot be sent to the console or saved for use in the
#'     web-page output? default=TRUE, set to FALSE to save it
#' @param q90 should the 90th quantiles be plotted? default=TRUE. If FALSE then
#'     the 95th quantiles are used.
#' @param intens if polys=TRUE then intens signifies the intensity of colour on
#'     a scale of 0 - 255. 127 is about 50 percent dense.
#' @param addleg add a legend or not. If not use "", if yes use any of topleft,
#'     topright, bottomleft, or bottomright, default = bottomright.
#' @param addmedian should each polygon also have its median plotted as a line.
#' @param add1line default = TRUE, adds a thicker dashed black line at 1.0,
#'     where the catch/msy ratio = 1.0
#'
#' @seealso {
#'   \link{doquantribbon}, \link{do_comp_outputs}
#' }
#'
#' @return a list containing the quantiles for the varname for each of the
#'     scenarios - invisibly
#' @export
#'
#' @examples
#' \dontrun{
#'  # a possible run where scenarios 1, 13, and 14 are compared.
#'    result <- do_comparison(rundir=rundir,postfixdir=postfixdir,outdir=outdir,
#'                            files=files,pickfiles=c(1,2,3),verbose=TRUE)
#' }
zoneribbon <- function(rundir,scenes,invar,glbc,varname,category,
                       console=TRUE,q90=TRUE,intens=127,addleg="topleft",
                       addmedian=3,add1line=TRUE) {
# rundir=""; scenes=scenes;invar=zcpue;glbc=glbc;varname="CPUE";category="cpue"
# console=TRUE; q90=TRUE;intens=127;addleg="topleft";addmedian=0
  zonequant <- makelist(scenes)
  glb <- glbc[[1]]
  reps <- glb$reps
  yrs <- glb$pyrnames
  nyrs <- length(yrs)
  nscen=length(scenes)
  maxq <- 0
  for (scen in 1:nscen) {  # scen=1
    zpick <- invar[[scen]]
    if (nrow(zpick) != nyrs) {
      zpick <- invar[[scen]][(glb$hyrs+1):(glb$hyrs + glb$pyrs),]
    }
    maxy <- getmax(zpick)
    zquant <- apply(zpick,1,quants)
    zonequant[[scen]] <- zquant
    if (q90) {
      zquant <- zquant[c(2,3,4),]
    } else {
      zquant <- zquant[c(1,3,5)]
    }
    maxy <- getmax(zquant)
    maxq <- ifelse((maxq < maxy),maxy,maxq)
  }
  filen=""
  if (!console) {
    filename <- paste0(varname,"_by_zone_ribbon.png")
    filen <- pathtopath(rundir,filename)
    caption <- paste0(varname," ribbon plot for Zones in ",
                      paste0(scenes," ",collapse=""))
  }
  plotprep(width=8,height=4.5,newdev=FALSE,filename=filen,verbose=FALSE)
  parset(cex=1.0,margin=c(0.35,0.5,0.05,0.05))
  plot(yrs,zonequant[[1]][3,],type="l",lwd=0,col=0,ylim=c(0,maxq),xlab="",
       ylab=paste0("Zonewide ",varname),panel.first=grid(),cex.lab=1.25,yaxs="i")
  for (i in 1:nscen) { # i =2
    if (q90) {
      poldat <- makepolygon(zonequant[[i]][2,],zonequant[[i]][4,],yrs)
      polygon(poldat,col=RGB(i,alpha=intens))
    } else {
      poldat <- makepolygon(zonequant[[i]][1,],zonequant[[i]][5,],yrs)
      polygon(poldat,col=RGB(i,alpha=intens))
    }
    if (addmedian > 0) lines(yrs,zonequant[[i]][3,],lwd=addmedian,col=i)
  }
  if (add1line) abline(h=1,lwd=1,col=1,lty=2)
  if (nchar(addleg) > 0) {
    if (addleg %in% c("topleft","topright","bottomleft","bottomright")) {
      legend(addleg,legend=scenes,col=1:nscen,lwd=3,bty="n",cex=1.2)
    } else {
      cat("Inappropriate legend location given, see ?doquantribbon \n")
    }
  }
  if (!console) {
    addplot(filen=filen,rundir=rundir,category=category,caption=caption)
  }
  return(invisible(zonequant))
} # end of zoneribbon

