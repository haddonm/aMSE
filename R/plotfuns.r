

#' @title compareCPUE plots the historical cpue against conditioned cpue
#'
#' @description compareCPUE generates a plot of the historical cpue and compares
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
compareCPUE <- function(histCE,saucpue,glb,rundir,filen="") {
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
  doplots=getparplots(nsau)
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
  return(ssq)
} # end of compareCPUE

#' @title diagnosticsproj plots a series of diagnostics to DiagProj
#'
#' @description diagnosticsproj provides a series of plots and results that
#'     illustrate the properties of the projections. These include the
#'     residuals between SAU actual catches and their predicted catches. But
#'     also a series of plots of the projections for only nrep trajectories to
#'     illustrate that the dynamic variables are changing through time in a
#'     plausible or 'realistic' manner.
#'
#' @param zonePsau the SAU scale object containing the dynamics results
#' @param glb the global constants object
#' @param rundir the rundir for the scenario
#' @param nrep the number of replicate trajectories to plot; default=3
#'
#' @return Nothing but it does add some plots to rundir
#' @export
#'
#' @examples
#' print("wait on suitable data sets")
diagnosticsproj <- function(zonePsau,glb,rundir,nrep=3) {
  # zonePsau=sauout$zonePsau;glb=glb;rundir=rundir;nrep=3
  reps <- glb$reps
  nsau <- glb$nSAU
  hyrs <- glb$hyrs
  pickyr <- (hyrs+1):(hyrs+glb$pyrs)
  yrnames <- glb$pyrnames
  plotvar <- function(var1,nsau,saunames,nrep,filen,caption,label,var2=NULL) {
    # a helper function that plots nrep trajectories from invar
    #var1=catch;nsau=nsau;saunames=saunames;nrep=3;filen="";caption="";label="_Actual_Catch";var2=acatch
    plotprep(width=7,height=7,filename=filen,cex=1.0,verbose=FALSE)
    parset(plots=getparplots(nsau))
    for (sau in 1:nsau) { # sau =1
      pickrep <- sample(1:reps,nrep,replace=FALSE)
      ymax <- getmax(var1[pickyr,sau,pickrep])
      ylabel <- paste0(saunames[sau],label)
      plot(yrnames,var1[pickyr,sau,pickrep[1]],type="l",lwd=2,ylim=c(0,ymax),
           ylab=ylabel,xlab="Years",panel.first=grid())
      for (i in 2:nrep) lines(yrnames,var1[pickyr,sau,pickrep[i]],lwd=2,col=i)
      if (is.numeric(var2))
        for (i in 1:nrep) lines(yrnames,var2[pickyr,sau,pickrep[i]],lwd=2,col=i,lty=3)
    } # end of actual catches
    addplot(filen,rundir=rundir,category="DiagProj",caption)
  } # end of internal function plotvar
  #  zonePsau=sauout$zonePsau; filen=""
#  nsau <- glb$nSAU
  saunames <- glb$saunames
  # plot catch residuals
  catch <- zonePsau$catch
  acatch <- zonePsau$acatch
  filen <- filenametopath(rundir,"catch-acatch_residuals.png") # filen=""
  caption <- paste0("The residuals between the actual SAU catches and the",
                    "predicted catches from the HCR.")
  plotprep(width=7,height=7,filename=filen,cex=1.0,verbose=FALSE)
  resid <- catch[pickyr,,] - acatch[pickyr,,]
  parset(plots=getparplots(nsau),margin=c(0.3,0.35,0.05,0.05),
         outmargin = c(1,1,0,0))
  for (sau in 1:nsau)
    hist(resid[,sau,],breaks=25,main="",ylab=saunames[sau],xlab="",
         panel.first=grid())
  mtext("Difference between Observed and Expected SAU Catches",
        side=1,outer=TRUE,cex=1.1,line=-0.1)
  mtext("Frequency",side=2,outer=TRUE,cex=1.1,line=-0.1)
  addplot(filen,rundir=rundir,category="DiagProj",caption)
  # plot the limited trajectory plots # filen=""
  filen <- filenametopath(rundir,paste0("actual_catch_projections_",nrep,".png"))
  caption <- paste0(nrep," projections of actual catches, dotted lines are",
                    " the aspirational catches.")
  plotvar(catch,nsau,saunames,nrep,filen,caption,"_Actual_Catch",var2=acatch)

  cpue <- zonePsau$cpue
  filen <- filenametopath(rundir,paste0("cpue_projections_",nrep,".png"))
  caption <- paste0(nrep," cpue projections")
  plotvar(cpue,nsau,saunames,nrep,filen,caption,"_CPUE")

} # end of diagnosticsproj

#' @title dosau plots the conditioning history for the dynamics
#'
#' @description dosau plots the deplsB, cpue, matureB, catch, harvestR, and
#'     recruits for a given SAU during the years of conditioning. This aims to
#'     assist the conditioning process by illustrating the state of the sau
#'     during and at the end of conditioning on the fishery. If extra is TRUE
#'     then it also plots the exploitable biomass and the depletion of the
#'     exploitable biomass. A loess fit is also plotted onto the recruitments
#'     plot to give some insight into the recruitment deviates. Unfortunately,
#'     the 'true' predicted recruitment levels (without variation) cannot be
#'     easily estimated because the implementation of larval movement disturbs
#'     how many recruits are present each year in each population and hence
#'     each SAU.
#'
#' @param inzone the conditioned zone (zoneDD) after it has been converted to
#'     sau scale by using getsauzone
#' @param glb the global constants object
#' @param picksau which sau should be plotted
#' @param histCE the historical cpue series from condC
#' @param yrnames the years of historical catches eg 1973:2019
#' @param extra should exploitable B and depleB be plotted? default=FALSE
#'
#' @return nothing but it does generate a plot
#' @export
#'
#' @examples
#' print("wait on suitable data-sets")
dosau <- function(inzone,glb,picksau,histCE,yrnames,extra=FALSE) {
  # inzone=sauZone;glb=glb;picksau=sau;histCE=histCE;yrnames=glb$hyrnames; extra=FALSE
  histceyr <- as.numeric(rownames(histCE))
  saunames <- glb$saunames #as.numeric(glb$saunames)
  saunum <- 1:glb$nSAU
  indexsau <- which(saunum == picksau)
  label=c("deplsB","cpue","matB","catch","harvestR","recruit")
  if (extra) {
    label=c("deplsB","cpue","matB","catch","harvestR","recruit","expB","depleB")
  }
  doplots=getparplots(length(label))
  nvar <- length(label)
  matB <- inzone$matB
  nyrs <- dim(matB)[1]
  parset(plots=doplots,margin=c(0.3,0.4,0.05,0.05),outmargin=c(1,0,1,0),
         byrow=FALSE)
  for (invar in 1:nvar) {  #  invar=1
    dvar <- inzone[[label[invar]]][,indexsau]
    ymax <- getmax(dvar)
    plot(yrnames,dvar,type="l",lwd=2,col=4,panel.first=grid(),yaxs="i",
         ylab=paste0(label[invar]),ylim=c(0,ymax),xlab="",cex=1.0)
    if (label[invar] == "cpue") lines(histceyr,histCE[,indexsau],lwd=2,col=3)
    if (label[invar] == "deplsB") {
      abline(h=0.2,col=2,lwd=1)
      pickyr <-round(nyrs/2)
      text(x=yrnames[pickyr],y=0.9*ymax,round(tail(dvar,1),4),cex=1.1,pos=4)
    }
    if (label[invar] == "harvestR") abline(h=c(0.4,0.75),col=2,lwd=1,lty=2)
    if (label[invar] == "recruit") {
      loessfit <- loess(dvar ~ glb$hyrnames,span=0.625)
      lines(loessfit$x,loessfit$fitted,lwd=2,col=2)
    }
  }
  mtext(paste0("Years ",glb$saunames[picksau]),side=1,outer=TRUE,cex=1.1,line=-0.1)
} # end of dosau

#' @title dosauplot generates the plot of the chosen variable for each sau
#'
#' @description dosauplot takes the results of the prepareproject and the
#'     projection function, it combines the trajectories across the years of
#'     conditioning and the projection years, for all replicates, and plots all
#'     replicates for the years from 'startyr' to the final year of the
#'     projections. It includes the median and inner 90% CI, and, for the CPUE
#'     plot, if the historical conditioning CPUE is included in histCE, it also
#'     includes the original cpue series. Beware of trying to plot the
#'     historical cpue in years where there are no historical cpue data,
#'     although there are error captures at work to limit the pot to available.
#'
#' @param ylabel the variable name to be plotted, used as a y-axis label
#' @param postrep the replicate values from the replicated projection zoneDP for
#'     the variable selected taken from the output of 'zonetosau'.
#' @param glb the global constants object
#' @param startyr the index for the first year of the conditioning dynamics to
#'     include
#' @param addCI should quantile confidence bounds be included on the plots,
#'     default=FALSE
#' @param CIprobs the quantiles used to generate the confidence intervals. The
#'     default = c(0.05,0.5,0.95)
#' @param histCE the historical cpue data if included, default = NULL
#' @param addCE should historical cpue be included on the plots, default=FALSE
#'
#' @return invisibly, a list of CI and median for each SAU
#' @export
#'
#' @examples
#' print("wait on suitable built in data sets")
dosauplot <- function(ylabel,postrep,glb,startyr,addCI=FALSE,
                      CIprobs=c(0.05,0.5,0.95),histCE=NULL,addCE=FALSE) {
  # ylabel="Catch";postrep=zonePsau[["catch"]];glb=glb; addCE=FALSE
  # startyr=startyr;histCE=histCE; CIprobs=c(0.05,0.5,0.95); addCI=TRUE
  label <- glb$saunames
  nsau <- glb$nSAU
  sauCI <- vector("list",nsau)
  names(sauCI) <- glb$saunames
  hyrs <- glb$hyrs
  allyrs <- hyrs + glb$pyrs
  yrnames <- c(glb$hyrnames,glb$pyrnames)
  reps <- glb$reps
  nplot <- getparplots(nsau)
  pyrnames <- glb$pyrnames
  yrnames <- c(glb$hyrnames,pyrnames)
  projyrs <- (hyrs + 1):allyrs
  parset(plots=nplot,byrow=FALSE)
  xvar <- yrnames[startyr:allyrs]
  for (sau in 1:nsau) {  # sau = 2
    ymax <- getmax(postrep[startyr:allyrs,sau,])
    plot(xvar,postrep[startyr:allyrs,sau,1],type="l",lwd=1,col="grey",panel.first=grid(),
         ylab=paste0(ylabel,"  ",label[sau]),ylim=c(0,ymax),xlab="Years",yaxs="i")
    for (i in 2:reps)
      lines(xvar,postrep[startyr:allyrs,sau,i],lwd=1,col="grey")
    CI <- apply(postrep[projyrs,sau,],1,quantile,probs=CIprobs)
    sauCI[[sau]] <- CI
    lines(yrnames[projyrs],CI[2L,],lwd=2L,col=4L)
    if (addCI) {
      lines(yrnames[projyrs],CI[1L,],lwd=1L,col=2L)
      lines(yrnames[projyrs],CI[3L,],lwd=1L,col=2L)
    }
    if (addCE) {  # add the historical cpue line
      histyr <- as.numeric(rownames(histCE))
      pickce <- which(histyr %in% xvar)
      lines(histyr[pickce],histCE[pickce,sau],lwd=2L,col=3L)
    }
    abline(v=(yrnames[hyrs]+0.5),col=2L,lty=2L)
    if (addCE) abline(h=150,col=2,lty=2)
  }
  return(invisible(sauCI))
} # end of dosauplot


#' @title onesau plots the dynamics for a single SAU
#'
#' @description onesau plots the details of matureB, exploitB, catch, acatch,
#'     harvestR, cpue, recruit, deplsB, depleB for a single SAU on one plot
#'
#' @param prerep the zoneDsau object that represents the SAU data
#' @param postrep the zonePsau object that represents the SAU data
#' @param glb the global constants object
#' @param startyr the year from which to begin the conditioned dynamics
#' @param picksau which sau should be plotted
#' @param addCI should the 90 percent CI be added, default=FALSE
#' @param CIprobs what CI should be fitted, default=c(0.05,0.5,0.95)
#' @param histCE should the historical CPUE be added to the cpue plot
#'
#' @return a list of the CI for each variable
#' @export
#'
#' @examples
#' print("Wait on suitable data sets")
onesau <- function(prerep,postrep,glb,startyr,picksau,addCI=FALSE,
                   CIprobs=c(0.05,0.5,0.95),histCE=NULL) {
  saunames <- glb$saunames
  sauindex <- which(saunames == picksau)
  label=c("matureB","exploitB","catch","acatch", "harvestR","cpue",
          "recruit","deplsB","depleB")
  nvar <- length(label)
  prematB <- prerep$matureB
  postmatB <- postrep$matureB
  preyrs <- dim(prematB)[1]
  projyrs <- dim(postmatB)[1]
  allyrs <- preyrs + projyrs
  reps <- dim(prematB)[3]
  varCI <- vector("list",nvar)
  names(varCI) <- label
  if (is.numeric(histCE)) ceyr <- startyr:preyrs
  parset(plots=getparplots(length(label)),margin=c(0.3,0.4,0.05,0.05),outmargin=c(1,0,1,0),
         byrow=FALSE)
  for (invar in 1:nvar) {  #  invar=1
    premat <- prerep[[label[invar]]][,sauindex,]
    postmat <- postrep[[label[invar]]][,sauindex,]
    ymax <- getmax(rbind(premat[startyr:preyrs,],postmat[,]))
    traj <- c(premat[startyr:preyrs,1],postmat[,1])
    plot(startyr:allyrs,traj,type="l",lwd=1,col="grey",panel.first=grid(),
         ylab=paste0(label[invar]),ylim=c(0,ymax),xlab="",cex=1.0)
    for (i in 2:reps)
      lines(startyr:allyrs,c(premat[startyr:preyrs,i],postmat[,i]),
            lwd=1,col="grey")
    CI <- apply(postmat[,],1,quantile,probs=CIprobs)
    varCI[[invar]] <- CI
    lines((preyrs+1):allyrs,CI[2,],lwd=2,col=4)
    if (addCI) {
      lines((preyrs+1):allyrs,CI[1,],lwd=1,col=2)
      lines((preyrs+1):allyrs,CI[3,],lwd=1,col=2)
    }
    abline(v=(preyrs+0.5),col=2)
    if ((is.numeric(histCE)) & (label[invar] == "cpue")) {
      nhistce <- dim(histCE)[1]
      if (length(ceyr) > nhistce) ceyr <- (preyrs - nhistce + 1):preyrs
      oldce <- tail(histCE[,sauindex],length(ceyr))
      lines(ceyr,oldce,lwd=2,col=3)
    }
  }
  mtext(paste0("SAU ",picksau),side=3,outer=TRUE,cex=1.0)
  mtext("Years",side=1,outer=TRUE,cex=1.1,line=-0.1)
  return(invisible(varCI))
} # end of onesau

#' @title onezoneplot plots out one variable from the zone
#'
#' @description once the projected zone has been summarized by poptozone one
#'     can add plots to the results using plotZone. This in turn uses
#'     onezoneplot to plot each variable in turn.
#'
#' @param invar which zone variable to plot
#' @param rundir the run directory for the results
#' @param glb the global constants object
#' @param CIprobs the quantiles to use for the CI.
#' @param varname the name of the variable for use in labels
#' @param startyr in what year number should the plot begin.
#' @param addfile should a png file be added to the results, default = TRUE
#'
#' @return it generates a png file of the plot if addfile remains TRUE
#' @export
#'
#' @examples
#' print("wait on suitable data")
onezoneplot <- function(invar,rundir,glb,CIprobs,varname,startyr,addfile=TRUE) {
  #  invar=inzone$catch;rundir=rundir;glb=glb;CIprobs=c(0.05,0.5,0.95);addfile=TRUE;startyr=(glb$hyrs + 1)
  # varname="Catch"
  yrnames <- c(glb$hyrnames,glb$pyrnames)
  allyr <- glb$hyrs + glb$pyrs
  pickyr <- startyr:allyr
  reps <- glb$reps
  namefile <- paste0("Projected_Zonal_",varname,".png")
  filen <- filenametopath(rundir,namefile)
  if (!addfile) filen <- ""
  plotprep(width=7,height=4,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  parset(plots=c(1,1))
  ymax <- getmax(invar[pickyr,])
  CI <- apply(invar,1,quantile,probs=CIprobs)
  plot(yrnames[pickyr],invar[pickyr,1],type="l",lwd=1,col="grey",ylab=varname,
       xlab="Years",panel.first=grid(),ylim=c(0,ymax),yaxs="i")
  for (i in 1:reps) lines(yrnames[pickyr],invar[pickyr,i],col="grey",lwd=1)
  lines(yrnames[pickyr],CI[1,pickyr],lwd=1,col=2)
  lines(yrnames[pickyr],CI[3,pickyr],lwd=1,col=2)
  lines(yrnames[pickyr],CI[2,pickyr],lwd=2,col=4)
  caption <- paste0("Projected ",varname," across all replicates.")
  if (addfile) addplot(filen,rundir=rundir,category="zonescale",caption)
} # end of onezoneplot


#' @title plotCNt plots the historical conditioning Nt or catchN for all years
#'
#' @description plotCNt the historical conditioning on the fisheries data leads
#'     to a three dimensional array of numbers-at-size x years x numpop. We use
#'     'prepareDDNt' to convert this to a 3D array by SAU. PlotCNt then plots,
#'     for each SAU, numbers-at-size for all years of conditioning with the
#'     first year emphasized in black and the last year in red.
#'
#' @param Nt the 3D array of numbers-at-size x years x SAU from prepareDDNt
#' @param glb the global constants object
#' @param vline add a vertical dashed line, perhaps at the LML, default=NULL
#' @param start from what size class should the maximum Y value be measured.
#'     The default = 3 (= 6mm)
#'
#' @return nothing but it does plot a graph
#' @export
#'
#' @examples
#'  print("wait on suitable built in data-sets")
plotCNt <- function(Nt,glb,vline=NULL,start=3) {
  #   Nt <- lcomp$Nt; year=1;start=3; glb=glb;
  Nclass <- glb$Nclass
  midpts <- glb$midpts
  nsau <- glb$nSAU
  saunames <- glb$saunames
  endyr <- dim(Nt)[2]
  parset(plots=getparplots(nsau),margin = c(0.25, 0.45, 0.05, 0.05),byrow=FALSE,
         outmargin = c(1.1, 0, 0, 0))
  for (sau in 1:nsau) { #  sau=1
    label <- paste0("N '000s ",saunames[sau])
    ymax <- getmax(Nt[start:Nclass,,sau])/1000.0
    plot(midpts,Nt[,1,sau]/1000.0,type="l",lwd=1,panel.first=grid(),
         ylim=c(0,ymax),ylab=label,xlab="")
    for (i in 2:(endyr-1)) lines(midpts,Nt[,i,sau]/1000.0,lwd=1,col="grey")
    lines(midpts,Nt[,1,sau]/1000.0,lwd=2,col=1)
    lines(midpts,Nt[,endyr,sau]/1000.0,lwd=2,col=2)
    if (!is.null(vline)) abline(v=vline,lwd=1,lty=2)
  }
  mtext("Shell Length mm",side=1,outer=TRUE,cex=1.1)
} # end of plotCNt


#' @title plotconditioning illustrate the zones dynamics after conditioning
#'
#' @description plotconditioning converts the conditioned dynamic zone object
#'     from population scale to sau scale and then plots the various components
#'     to depict their state after conditioning. This is designed to place such
#'     plots into the rundir to assist in the conditioning process.
#'
#' @param zoneDD the conditioned dynamic zone object
#' @param glb the global constants object
#' @param zoneC the constant zone object
#' @param histCE the historical CPUE matrix
#' @param rundir the rundir for the given scenario
#'
#' @return the sau scale zone dynamics, invisibly
#' @export
#'
#' @examples
#' print("wait on suitable data sets")
plotconditioning <- function(zoneDD,glb,zoneC,histCE,rundir) {
  # zoneDD=zoneDD;glb=glb;zoneC=zoneC;histCE=condC$histCE;rundir=rundir
  sauindex <- glb$sauindex
  popB0 <- getlistvar(zoneC,"B0")
  B0 <- tapply(popB0,sauindex,sum)
  popExB0 <- getlistvar(zoneC,"ExB0")
  ExB0 <- tapply(popExB0,sauindex,sum)
  sauZone <- getsauzone(zoneDD,glb,B0=B0,ExB0=ExB0)
  allsau <- 1:glb$nSAU #
  snames <- glb$saunames
  yrnames <- glb$hyrnames
  ssq <- compareCPUE(histCE,sauZone$cpue,glb,rundir,filen="compareCPUE.png")
  for (sau in allsau) {  #  sau=allsau[1]; filen=""
    filen <- filenametopath(rundir,paste0(snames[sau],"_conditioned.png"))
    plotprep(width=8,height=8,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
    dosau(sauZone,glb,picksau=sau,histCE=histCE,yrnames=glb$hyrnames)
    caption <- "Dynamics during conditioning on the fishery."
    addplot(filen,rundir=rundir,category="condition",caption)
  }
  return(invisible(list(sauZone=sauZone,ssq=ssq)))
} # end plotconditioning

#' @title plothsstats plots some HS performance statistics
#'
#' @description plothstats generates plots for the sum of catches for the
#'     first 5 and 10 years of the projections. This is really only
#'     interesting when compared with the same statistics from a different
#'     operating model definition or different harvest strategy. The input
#'     data is generated during the running of do_MSE.
#'
#' @param rundir the scenario's results directory
#' @param hsstats the statistics about the harvest strategies performance
#' @param glb the global constants object
#'
#' @return nothing but it does generate two plots into rundir
#' @export
#'
#' @examples
#' print("wait on suitable in ternal data sets")
plothsstats <- function(rundir,hsstats,glb) {
  # hsstats=out$HSstats; rundir=rundir; glb=out$glb
  sum5 <- hsstats$sum5
  sum10 <- hsstats$sum10
  labelnames <- colnames(sum10)
  nsau <- glb$nSAU
  filen <- filenametopath(rundir,"sum_5yr_proj_catches.png")
 # filen=""
  plotprep(width=8,height=7,filename=filen,cex=1.0,verbose=FALSE)
  parset(plots=getparplots(ncol(sum5)))
  for (i in 1:(nsau+1)) {
    label <- paste0("Tonnes     ",labelnames[i])
    hist(sum5[,i],breaks=20,main="",xlab=label)
  }
  caption <- paste0("The sum of the first 5 years of he projections for each ",
                    "SAU and the complete zone")
  addplot(filen,rundir=rundir,category="HSperf",caption)
  #now 10 years
  filen <- filenametopath(rundir,"sum_10yr_proj_catches.png")
  #filen=""
  plotprep(width=8,height=7,filename=filen,cex=1.0,verbose=FALSE)
  parset(plots=getparplots(ncol(sum10)))
  for (i in 1:(nsau+1)) {
    label <- paste0("Tonnes     ",labelnames[i])
    hist(sum10[,i],breaks=20,main="",xlab=label)
  }
  caption <- paste0("The sum of the first 10 years of he projections for ",
                    "each SAU and the complete zone")
  addplot(filen,rundir=rundir,category="HSperf",caption)
} # end of plothsstats

#' @title plotNt plots the size-composition for each SAU
#'
#' @description plotNt accepts size-composition data (either the population Nt,
#'     of the catchN), and plots the replicate size-distributions for each SAU.
#'     It then puts the first year's size-distribution on top and if the medcol
#'     that defines the colour of hte mdeian of the replicates is different
#'     from 0 (not plotted) then it plots the median as well. A contrasting
#'     colour choice is 4 = Blue.
#'
#'
#' @param Nt the 4 dimensional array of numbers-at-size Nclass,years,nsau,reps.
#' @param year which year from dimension 2 should be plotted?
#' @param glb the global constants object
#' @param start from what size class should the maximum Y value be measured.
#'     The default = 3 (= 6mm)
#' @param medcol colour of the median line; default=0 = not plotted
#'
#' @return invisibly the matrix of median values, Nclass x nsau
#' @export
#'
#' @examples
#' print("wait on suitable built in data-sets")
plotNt <- function(Nt,year,glb,start=3,medcol=0) {
  #   Nt <- zonePsau$Nt; origNt=zoneD$Nt;year=25;start=3; glb=glb; newdev=FALSE
  Nclass <- glb$Nclass
  midpts <- glb$midpts
  reps <- dim(Nt)[4]
  nsau <- glb$nSAU
  saunames <- glb$saunames
  saumedians <- matrix(0,nrow=Nclass,ncol=nsau,dimnames=list(midpts,saunames))
  parset(plots=getparplots(nsau),byrow=FALSE)
  for (sau in 1:nsau) { #  sau=1
    sdat <- Nt[,year,sau,]
    saumedians[,sau] <- apply(sdat,1,median)
    ymax <- getmax(sdat[start:Nclass,])/1000.0
    label <- paste0("N '000s  ",saunames[sau])
    plot(midpts,sdat[,1]/1000.0,type="l",lwd=1,col="grey",panel.first=grid(),
         ylim=c(0,ymax),ylab=label)
    for (i in 1:reps) lines(midpts,sdat[,i]/1000.0,col="grey")
   # lines(midpts,sdat[,1]/1000.0,lwd=2,col=2)
    lines(midpts,saumedians[,sau]/1000.0,lwd=2,col=medcol)
  }
  return(invisible(saumedians))
} # end of plotNt

#' @title plotprod graphs a selected productivity variable
#'
#' @description plotprod graphs a selected productivity variable from
#'     the choice of ExB exploitable biomass, MatB mature of spawning
#'     biomass = Bmsy, AnnH the actual annual harvest rate Hmsy, Catch
#'     the yield at Bmsy and Hmsy = MSY, Deplet the mature biomass
#'     depletion level, and RelCE the relative cpue at MSY
#'
#' @param product the output from doproduction
#' @param xname the name of the productivity variable for the x-axis,
#'     defaults to MatB
#' @param yname the name of the y-xis variable, default=Catch=Yield
#' @param xlimit default=NA, enables the range of the x-axis to be
#'     constrained, for example using xlimit=c(0.15,0.4)
#' @param xlab the x-axis label, default=""
#' @param ylab the y-axis label, default=""
#' @param font the type of font used in the plot. default = 7, bold
#'     Times, 6 is not bold, 1 is sans-serif, 2 bold sans-serif
#' @param filename the complete path and filename of where to save the
#'     png plot file. default is empty, meaning no file is produced.
#' @param devoff a boolean to allow the plot device to remain open for
#'     the user to add more components. Of course, if this is set to
#'     FALSE then it is up to the user to use dev.off
#'
#' @return invisibly a list of the x and y matrices plotted
#' @export
#'
#' @examples
#' data(zone)
#' product <- zone$product
#' plotprod(product)
#' stat <- findmsy(product)
#' abline(h=stat[,"Catch"],col=1:6,lwd=2)
#' abline(v=stat[,"MatB"],col=1:6,lwd=2)
#' print(stat)
plotprod <- function(product,xname="MatB",yname="Catch",xlimit=NA,
                     xlab="Mature Biomass t",ylab="Production t",
                     font=7,filename="",devoff=FALSE) {
  x <- product[,xname,]
  y <- product[,yname,]
  numpop <- ncol(x)
  maxy <- getmax(y)
  if (length(xlimit) == 1) xlimit <- c(0,getmax(x))
  plotprep(width=7,height=4,newdev=FALSE,filename=filename,cex=0.9,
           verbose=FALSE)
  parset(font=font)
  plot(x[,1],y[,1],type="l",lwd=2,col=1,ylim=c(0,maxy),xlim=xlimit,
       xlab=xlab,ylab=ylab,panel.first=grid(),yaxs="i")
  if (numpop > 1) {
    for (i in 2:numpop) lines(x[,i],y[,i],lwd=2,col=i)
    legend("topright",paste0("P",1:numpop),lwd=3,col=c(1:numpop),
           bty="n",cex=1.2)
  }
  if ((nchar(filename) > 0) & (devoff)) dev.off()
  return(invisible(list(xname=x,yname=y)))
} # end of plotprod


#' @title plotZone plots the projected TAC across all replicates
#'
#' @description plotZone plots out the projected TAC for all replicates. The
#'     first value is the sum of the final catches in the conditioning period.
#'
#' @param inzone the dynamic zone projection object
#' @param rundir the rundir for the scenario
#' @param glb the global constants object
#' @param startyr which year number to begin the plot
#' @param CIprobs the quantiles used to generate the confidence intervals. The
#'     default = c(0.05,0.5,0.95)
#' @param addfile should the plots be added to the rundir =TRUE or only sent to
#'     the console = FALSE. Default = TRUE
#'
#' @return nothing but it adds a plot to the zonescale tab in the rundir
#' @export
#'
#' @examples
#' print("wait on suitable datasets")
#' # inzone=outzone;rundir=rundir;glb=glb;CIprobs=c(0.05,0.5,0.95);addfile=TRUE; startyr=40
plotZone <- function(inzone,rundir,glb,startyr,CIprobs=c(0.05,0.5,0.95),addfile=TRUE) {
  onezoneplot(inzone$catch,rundir,glb,CIprobs=CIprobs,"Catch",startyr=startyr,
              addfile=addfile)
  onezoneplot(inzone$cpue,rundir,glb,CIprobs=CIprobs,"CPUE",startyr=startyr,
              addfile=addfile)
  onezoneplot(inzone$deplsB,rundir,glb,CIprobs=CIprobs,"Mature_Biomass_Depletion",
              startyr=startyr,addfile=addfile)
  onezoneplot(inzone$matureB,rundir,glb,CIprobs=CIprobs,"Mature_Biomass",
              startyr=startyr,addfile=addfile)
  onezoneplot(inzone$harvestR,rundir,glb,CIprobs=CIprobs,"Harvest_Rate",
              startyr=startyr,addfile=addfile)
  onezoneplot(inzone$recruit,rundir,glb,CIprobs=CIprobs,"Recruitment",
              startyr=startyr,addfile=addfile)
  onezoneplot(inzone$TAC,rundir,glb,CIprobs=CIprobs,"TAC",startyr=startyr,
              addfile=addfile)
} # end of plotZone

#' @title plotzonesau generates a plot of the zone scale and sau scale
#'
#' @description plotsonesau enables the zone catch or catch weighted zone cpue
#'    (from catchweightCE) to be compared with the sau catch or the sau cpue
#'    within a single plot. This has value for gaining insights into the
#'    fishery dynamics of a zone and how it responds to fishing pressure.
#'    NOTE WELL: Currently, this function is setup for Tasmania's set of 8 SAU.
#'    Assuming other jurisdictions have more or less then the layout of the
#'    plots would need to be altered. Until these are known this function will
#'    likely only work as intended with Tasmanian data!
#'
#' @param zonetot a vector of either the total catch or zone summary cpue for a
#'     zone. Ideally, the vector should be named using the year in which each
#'     value was taken.
#' @param saudat either the catch or the cpue for each SAU (year x sau), again
#'     this matrix should have both row and column names
#' @param saunames a vector of nSAU names. Ideally very short ones to ensure
#'     they fit on each small plot
#' @param label the y-axis label for the zone plot
#' @param labelsau the y-axis label for each row of the SAU plots
#' @param side should the SAU label be at the top of bottom of each plot.
#'     default = 3 (= top of plot). The main alternative would be 1 (the bottom)
#'     but 2 (left side) and 4 (right side) also work.
#' @param sauscale should all SAU plots have the own y-axis scale = TRUE, the
#'     default. If you want each plot to have the same scale sauscale=FALSE.
#'
#' @return nothing but this does generate a plot.
#' @export
#'
#' @seealso catchweightCE
#'
#' @examples
#' print("wait on suitable internal dats sets")
plotzonesau <- function(zonetot,saudat,saunames,label,labelsau,side=3,
                        sauscale=TRUE) {
  yrs <- as.numeric(names(zonetot))
  nyrs <- length(yrs)
  addn <- nyrs - nrow(saudat)
  nsau <- ncol(saudat)
  layout(rbind(c(1,1,1,1),
               c(2,3,4,5),
               c(6,7,8,9)),
         heights=c(4,1.5,1.5))
  par(mai=c(0.3,0.45,0.05,0.05),oma=c(0,1,0,0.4),cex=1.0)
  #layout.show(9)
  ymax <- getmax(zonetot)
  plot(yrs,zonetot,type="l",lwd=2,xlab="",ylab=label,ylim=c(0,ymax),
       panel.first=grid(),yaxs="i")
  par(mai=c(0.25,0.3,0.05,0.05),cex=0.7)
  ymax <- getmax(saudat)
  for (i in 1:nsau) {
    if (sauscale) ymax <- getmax(saudat[,i])
    plot(yrs,putNA(saudat[,i],addn,0),type="l",lwd=2,xlab="",ylab="",ylim=c(0,ymax),
         yaxs="i",panel.first=grid())
    if ((i == 1) | (i == 5)) mtext(labelsau,side=2,outer=FALSE,line=1.5,cex=1.0)
    mtext(saunames[i],side=side,outer=FALSE,line = -1.1,cex=1.0)
  }
} # end of plotzonesau


#' @title sauplots generates the dosauplots for the dynamic variables
#'
#' @description sauplots generates the dosauplots for the dynamic variables
#'    including cpue, catch, acatch (aspirational catch), matureB, exploitB,
#'    recruit, and harvestR
#'
#' @param zoneDP the dynamic object produced by the projections
#' @param NAS the numbers-at-size 4D arrays from doprojection
#' @param glb the object containing the global variables
#' @param rundir the results directory
#' @param B0 the B0 values by population use getvar(zoneC,"B0")
#' @param ExB0 the ExB0 values by population use getvar(zoneC,"ExB0")
#' @param startyr the index for the first year of the conditioning dynamics to
#'     include in the trajectory, to give startyr:indexoflastyear,eg startyr=40
#' @param addCI should confidence intervals be added to envelope plots,
#'     default=TRUE
#' @param histCE historical CPUE data used in CPUE plots, default=NULL
#' @param tabcat the name of the results website tab for the plots
#'
#' @return a list of lists of CI for each SAU and variable as well as the
#'     zoneDsau and zonePsau
#' @export
#'
#' @examples
#' print("wait on suitable internal data-sets")
sauplots <- function(zoneDP,NAS,glb,rundir,B0,ExB0,startyr,addCI=TRUE,
                     histCE=NULL,tabcat="projSAU") {
  # zoneDP=zoneDP;NAS=NAS;glb=glb;rundir=rundir;B0=B0;ExB0=ExB0;
  # startyr=48; addCI=TRUE;histCE=histCE; tabcat="projSAU"
  zonePsau <- zonetosau(zoneDP,NAS,glb,B0,ExB0)
  label <-  c("cpue","catch","acatch","matureB","exploitB","recruit","harvestR")
  out <- vector("list",length(label))
  names(out) <- label
  #CPUE
  filen <- filenametopath(rundir,"proj_cpue_SAU.png")  # filen=""
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  CI <- dosauplot("cpue",zonePsau[["cpue"]],glb,addCE=TRUE,
                  startyr=startyr,addCI=TRUE,histCE=histCE)
  caption <- "The CPUE projections for each SAU."
  addplot(filen,rundir=rundir,category=tabcat,caption)
  out[["cpue"]] <- CI
  #Catches
  filen <- filenametopath(rundir,"proj_catch_SAU.png")
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  CI <- dosauplot("catch",zonePsau[["catch"]],glb,
                  startyr=startyr,addCI=TRUE,histCE=NULL)
  caption <- "The catch projections for each SAU."
  addplot(filen,rundir=rundir,category=tabcat,caption)
  out[["catch"]] <- CI
  #Aspirational catches
  filen <- filenametopath(rundir,"proj_aspcatch_SAU.png")
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  CI <- dosauplot("acatch",zonePsau[["acatch"]],glb,
                  startyr=startyr,addCI=TRUE,histCE=NULL)
  caption <- paste0("The Aspirational catch projections for each SAU. Catches ",
                    "prior to HS are actual catches.")
  addplot(filen,rundir=rundir,category=tabcat,caption)
  out[["acatch"]] <- CI
  #MatureBiomass
  filen <- filenametopath(rundir,"proj_matureB_SAU.png")
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  CI <- dosauplot("matureB",zonePsau[["matureB"]],glb,
                  startyr=startyr,addCI=TRUE,histCE=NULL)
  caption <- "The mature biomass projections for each SAU."
  addplot(filen,rundir=rundir,category=tabcat,caption)
  out[["matureB"]] <- CI
  #exploitable biomass
  filen <- filenametopath(rundir,"proj_exploitB_SAU.png")
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  CI <- dosauplot("exploitB",zonePsau[["exploitB"]],glb,
                  startyr=startyr,addCI=TRUE,histCE=NULL)
  caption <- "The exploitable biomass projections for each SAU."
  addplot(filen,rundir=rundir,category=tabcat,caption)
  out[["exploitB"]] <- CI
  #recruitment
  filen <- filenametopath(rundir,"proj_recruit_SAU.png")
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  CI <- dosauplot("recruit",zonePsau[["recruit"]],glb,
                  startyr=startyr,addCI=TRUE,histCE=NULL)
  caption <- "The recruitment projections for each SAU."
  addplot(filen,rundir=rundir,category=tabcat,caption)
  out[["recruit"]] <- CI
  filen <- filenametopath(rundir,"proj_harvestR_SAU.png")
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  CI <- dosauplot("harvestR",zonePsau[["harvestR"]],glb,
                  startyr=startyr,addCI=TRUE,histCE=NULL)
  caption <- "The Harvest Rate projections for each SAU."
  addplot(filen,rundir=rundir,category=tabcat,caption)
  out[["harvestR"]] <- CI
  return(invisible(list(outCI=out,zonePsau=zonePsau)))
}# end of sauplots

#' @title saurecdevs plots the recruitment deviates for each SAU
#'
#' @description saurecdevs enables a visual comparison of the 'fitted'
#'     recruitment deviates for each SAU.
#'
#' @param recdevs the matrix of recdevs (out$condC$recdevs)
#' @param glb the globals object
#' @param rundir the scenario directory for results storage
#' @param filen the filename (no path needed) to be used in the rundir
#'
#' @return nothing, but it does generate a plot and save a png file if filen
#'     is populated
#' @export
#'
#' @examples
#' print("wait on suitable data sets")
saurecdevs <- function(recdevs,glb,rundir,filen="") {
  if (nchar(filen) > 0) filen <- filenametopath(rundir,filen)
  years <- as.numeric(rownames(recdevs))
  arecdevs <- abs(recdevs)
  nsau <- glb$nSAU
  label <- glb$saunames
  doplots=getparplots(nsau)
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,verbose=FALSE)
  parset(plots=doplots,margin=c(0.3,0.3,0.05,0.05),outmargin=c(0,1,0,0),
         byrow=FALSE)
  for (sau in 1:nsau) {
    ymax <- getmax(arecdevs[,sau])
    plot(years,arecdevs[,sau],type="l",lwd=2,col=1,xlab="",ylab="",
         panel.first=grid(),ylim=c(0,ymax),yaxs="i")
    #   lines(years,histCE[,sau],lwd=2,col=3)
    mtext(label[sau],side=1,line=-1.3,cex=1.25)
    abline(h=1.0,lwd=1,col=2)
  }
  mtext("Recruitment",side=2,line=-0.3,outer=TRUE,cex=1.25)
  if (nchar(filen) > 0) {
    caption <- "Historical Recruitment Deviations."
    addplot(filen,rundir=rundir,category="condition",caption)
  }
} # end of compareCPUE


