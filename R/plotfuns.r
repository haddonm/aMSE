

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
#' @param obscol the colour of the observed CPUE on the plots, default=2=red
#'
#' @return a vector of length nsau containing the ssq for each SAU
#' @export
#'
#' @examples
#' print("wait on suitable internal data sets")
#' # histCE=histCE;saucpue=sauZone$cpue;glb=glb;rundir=rundir;filen=filen
#' # sauZone=out$condout$sauZone; saucpue=sauZone$cpue; filen="";
#' # histCE=out$condC$histCE;  glb=out$glb;obscol=2
compareCPUE <- function(histCE,saucpue,glb,rundir,filen="",obscol=2) {
  if (nchar(filen) > 0) filen <- filenametopath(rundir,filen)
  nsau <- glb$nSAU
  hce <- TRUE
  if (is.null(histCE)) {
    hyrs <- glb$hyrs
    years <- 1:hyrs
    cpue <- saucpue[,1:nsau]
    hce <- FALSE
   } else {
    years <- as.numeric(rownames(histCE))
    hyrs <- glb$hyrnames
    pick <- match(years,hyrs)
    cpue <- as.matrix(saucpue[pick,1:nsau])
    rownames(cpue) <- years
  }
  label <- glb$saunames
  colnames(cpue) <- label
  ssq <- numeric(nsau)
  doplots=pickbound(nsau)
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,verbose=FALSE)
  parset(plots=doplots,margin=c(0.3,0.3,0.05,0.05),outmargin=c(0,1,0,0),
         byrow=FALSE)
  for (sau in 1:nsau) {  #  sau = 1
    ymax <- getmax(cpue[,sau])
    if (hce) {
      ssq[sau] <- sum((histCE[,sau] - cpue[,sau])^2,na.rm=TRUE)
      ymax <- getmax(c(cpue[,sau],histCE[,sau]))
    }
    plot(years,cpue[,sau],type="l",lwd=2,col=1,xlab="",ylab="",
         panel.first=grid(),ylim=c(0,ymax),yaxs="i")
   if (hce) {
     lines(years,histCE[,sau],lwd=2,col=2)
     lab2 <- paste0(label[sau],"  ",round(ssq[sau],1))
     mtext(lab2,side=1,line=-1.3,cex=1.25)
   }
  }
  legend("topright",c("Observed CE","Predicted CE"),col=c(2,1),lwd=3,bty="n",
         cex=1.0)
  mtext("CPUE",side=2,line=-0.3,outer=TRUE,cex=1.25)
  if (nchar(filen) > 0) {
    if (hce) {
      caption <- paste0("Comparison of Observed CPUE (red line) with that ",
                        "predicted by operating model (black line).")
      } else {
      caption  <- "Predicted CPUE from Operating Model."
    }
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
  # zonePsau=sauout;glb=glb;rundir=rundir;nrep=3
  reps <- glb$reps
  nsau <- glb$nSAU
  hyrs <- glb$hyrs
  pickyr <- (hyrs+1):(hyrs+glb$pyrs)
  yrnames <- glb$pyrnames
  plotvar <- function(var1,nsau,saunames,nrep,filen,caption,label,var2=NULL) {
    # a helper function that plots nrep trajectories from invar
    #var1=catch;nsau=nsau;saunames=saunames;nrep=3;filen="";caption="";label="_Actual_Catch";var2=acatch
    plotprep(width=7,height=7,filename=filen,cex=1.0,verbose=FALSE)
    parset(plots=pickbound(nsau))
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
  saunames <- glb$saunames
  # plot catch residuals
  catch <- zonePsau$catch
  acatch <- zonePsau$acatch
  filen <- filenametopath(rundir,"catch-acatch_residuals.png") # filen=""
  caption <- paste0("The residuals between the actual SAU catches and the",
                    "predicted catches from the HCR.")
  plotprep(width=7,height=7,filename=filen,cex=1.0,verbose=FALSE)
  resid <- catch[pickyr,,] - acatch[pickyr,,]
  parset(plots=pickbound(nsau),margin=c(0.3,0.35,0.05,0.05),
         outmargin = c(1,1,0,0))
  for (sau in 1:nsau) {
    if (length(dim(resid)) == 2) {
      hist(resid,breaks=25,main="",ylab=saunames[sau],xlab="")
    } else {
      hist(resid[,sau,],breaks=25,main="",ylab=saunames[sau],xlab="")
    }
  }
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
#' @param histCatch the historical catch series from condC
#' @param yrnames the years of historical catches eg 1973:2019
#' @param recdev the recdevs for a single SAU
#'
#' @return nothing but it does generate a plot
#' @export
#'
#' @examples
#' print("wait on suitable data-sets")
dosau <- function(inzone,glb,picksau,histCE,histCatch,yrnames,recdev) {
  # inzone=sauZone;glb=glb;picksau=sau;histCE=histCE;yrnames=glb$hyrnames; recdev=recdevs[,sau]
  hce <- TRUE
  if (is.null(histCE)) {
    hce <- FALSE
    yrnames <- 1:glb$hyrs
  } else {
    histceyr <- as.numeric(rownames(histCE))
  }
  saunames <- glb$saunames #as.numeric(glb$saunames)
  saunum <- 1:glb$nSAU
  indexsau <- which(saunum == picksau)
  label=c("deplsB","cpue","matB","catch","harvestR","depleB",
          "recruit")
  nvar <- length(label) + 2 # add an extra 1 for the recdevs
  doplots=pickbound(nvar)
  matB <- inzone$matB
  nyrs <- dim(matB)[1]
  parset(plots=doplots,margin=c(0.3,0.4,0.05,0.05),outmargin=c(1,0,1,0),
         byrow=FALSE)
  for (invar in 1:(nvar-2)) {  #  invar=2
    dvar <- inzone[[label[invar]]][,indexsau]
    ymax <- getmax(dvar)
    plot(yrnames,dvar,type="l",lwd=2,col=4,panel.first=grid(),yaxs="i",
         ylab=paste0(label[invar]),ylim=c(0,ymax),xlab="",cex=1.0)
    if ((label[invar] == "cpue") & (hce))
        lines(histceyr,histCE[,indexsau],lwd=2,col=3)
    if (label[invar] %in% c("deplsB","depleB")) {
      abline(h=0.2,col=2,lwd=1)
      pickyr <-round(nyrs/2)
      text(x=yrnames[pickyr],y=0.9*ymax,round(tail(dvar,1),4),cex=1.1,pos=4)
    }
    if (label[invar] == "harvestR") abline(h=c(0.4,0.75),col=2,lwd=1,lty=2)
    if (label[invar] == "recruit") {
      loessfit <- loess(dvar ~ yrnames,span=0.625)
      lines(loessfit$x,loessfit$fitted,lwd=2,col=2)
    }
  }
  arecdev <- abs(recdev)
  ymax <- getmax(arecdev)
  plot(yrnames,arecdev,type="l",lwd=2,col=4,panel.first=grid(),yaxs="i",
       ylab="recdevs",ylim=c(0,ymax),xlab="",cex=1.0)
  abline(h=1.0,col=2,lwd=1)
  mtext(paste0("Years ",glb$saunames[picksau]),side=1,outer=TRUE,cex=1.1,line=-0.1)
  # now plot the CPUE vs Catch
  histyr <- as.numeric(rownames(histCatch))
  label <- "Correlation"
  cedat <- histCE[,indexsau]
  cpuedat <- cedat[which(cedat > 0)]
  pickC <- which(histyr %in% as.numeric(names(cpuedat)))
  catdat <- histCatch[pickC,indexsau]
  ccfout <- ccf(x=catdat,y=cpuedat,type="correlation",
                ylab=label,plot=TRUE,xlab="Lag Years")
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
#' @param hlin a vector of horizontal reference lines, one for each sau, or NULL
#'     if a vector then a horizontal line will be added to each plot.
#'
#' @return invisibly, a list of CI and median for each SAU
#' @export
#'
#' @examples
#' print("wait on suitable built in data sets")
dosauplot <- function(ylabel,postrep,glb,startyr,addCI=FALSE,
                      CIprobs=c(0.05,0.5,0.95),histCE=NULL,addCE=FALSE,
                      hlin=NULL) {
  # ylabel="Catch";postrep=zonePsau[["harvestR"]];glb=glb; addCE=FALSE
  # startyr=startyr;histCE=NULL; CIprobs=c(0.05,0.5,0.95); addCI=NULL
  label <- glb$saunames
  nsau <- glb$nSAU
  sauCI <- vector("list",nsau)
  names(sauCI) <- glb$saunames
  hyrs <- glb$hyrs
  allyrs <- hyrs + glb$pyrs
  yrnames <- c(glb$hyrnames,glb$pyrnames)
  reps <- glb$reps
  nplot <- pickbound(nsau)
  pyrnames <- glb$pyrnames
  yrnames <- c(glb$hyrnames,pyrnames)
  projyrs <- (hyrs + 1):allyrs
  parset(plots=nplot,byrow=FALSE)
  xvar <- yrnames[startyr:allyrs]
  for (sau in 1:nsau) {  # sau = 1
    ymax <- getmax(postrep[startyr:allyrs,sau,])
    plot(xvar,postrep[startyr:allyrs,sau,1],type="l",lwd=1,col="grey",panel.first=grid(),
         ylab=paste0(ylabel,"  ",label[sau]),ylim=c(0,ymax),xlab="Years",yaxs="i")
    for (i in 2:reps)
      lines(xvar,postrep[startyr:allyrs,sau,i],lwd=1,col="grey")
    CI <- apply(postrep[projyrs,sau,],1,quantile,probs=CIprobs,na.rm=TRUE)
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
    if (!is.null(hlin)) abline(h=hlin[sau],lwd=2,lty=2,col=1)
  }
  return(invisible(sauCI))
} # end of dosauplot

#' @title finalcondyeardepletion plots and tabulates final depletion after conditioning
#'
#' @description finalcondyeardepletion generates histograms and tabulations of
#'     either the spawning biomass or exploitable biomass depletion levels in the
#'     final year of conditioning. This includes the varyrs of recruitment
#'     variation prior to projections. It plots a histogram for each sau and
#'     generates a table of quantiles at probs of 0, 0.05, 0.5, 0.95, and 1.0.
#'
#' @param rundir the scenario's directory
#' @param sauzone the output from zonetosau, found in sauplots. This is
#'     the projection results summarized into their respective sau rather than
#'     the individual populations
#' @param glb the globals object
#' @param deplvar which depletion variable to use, either 'sB' or 'eB' are the
#'     only options. Any other input will use 'sB'
#' @param console should the plot go to the console or be saved into rundir?
#'     default=TRUE
#'
#' @seealso{
#'    \link{zonetosau}, \link{sauplots}
#' }
#'
#' @return invisibly returns the matrix of quantiles for each sau
#' @export
#'
#' @examples
#' print("wait on data sets")
finalcondyeardepletion <- function(rundir,sauzone,glb,deplvar="sB",console=TRUE) {
  nsau <- glb$nSAU
  hyrs <- glb$hyrs
  saunames <- glb$saunames
  if (!(deplvar %in% c("sB","eB"))) deplvar <- "sB"
  depl <- switch(deplvar,"sB" = sauzone$deplsB,
                 "eB" = sauzone$depleB)
  deplquant <- matrix(0,nrow=nsau,ncol=5,
                      dimnames=list(saunames,c("0%","5%","50%","95%","100%")))
  filen <- ""
  if (!console) {
    nfile <- paste0("finalcondyear_",deplvar,"_depletion.png")
    filen <- filenametopath(rundir,nfile)
    caption <- "Histograms of the depletion in the final year of conditioning."
  }
  plotprep(width=8,height=9,newdev=TRUE,filename=filen,cex=0.9,verbose=FALSE)
  parset(plots=pickbound(nsau),byrow=FALSE,margin=c(0.3,0.4,0.05,0.05),
         outmargin=c(1,1,0,0))
  for (i in 1:nsau) {
    deplquant[i,] <- quantile(depl[hyrs,i,],probs=c(0,0.05,0.5,0.95,1),na.rm=TRUE)
    hist(depl[hyrs,i,],main="",xlab="",ylab=saunames[i])
    abline(v=deplquant[i,3],lwd=2,col=2)
  }
  label <- paste0("Final Conditioning Year ",deplvar," Depletion")
  mtext(label,side=1,outer=TRUE,cex=1.1,line=-0.2)
  mtext("Frequency",side=2,outer=TRUE,cex=1.1,line=-0.2)
  if (!console) {
    addplot(filen,rundir=rundir,category="condition",caption)
    tabname <- paste0("finalcondyear_",deplvar,"_depletion.csv")
    addtable(round(deplquant,4),tabname,rundir,category="condition",caption=
               "Quantiles on depletion after conditioning on historical catches.")
  }
  return(invisible(deplquant))
} # end of finalcondyeardepletion

#' @title historicalplots a wrapper function to call plots of fishery plots
#'
#' @description historicalplots is a wrapper used by do_MSE and do_condition,
#'     to plot out details of the historical fishery data used in the
#'     conditioning. It currently plots out the historical catches and CPUE.
#'
#' @param rundir the directory where all results are stored
#' @param condC the R object contaiining the historical fishery data generated
#'     by the readctrlfile function
#' @param glb the global variables object
#' @param commonscale should the SAU plots share a common y-axis scale,
#'     default =  FALSE
#' @param proportion should the plot be of proportion relative to the maximum
#'     catch through each time-series or as raw tonnes x year, default=FALSE
#'
#' @return nothing but it does add plot files to the rundir
#' @export
#'
#' @examples
#' print("wait on suitable data sets")
historicalplots <- function(rundir,condC,glb,commonscale=FALSE,proportion=FALSE) {
  # condC=out$condC; glb=out$glb; commonscale=FALSE; proportion=TRUE
  hyrnames <- glb$hyrnames
  hyrs <- glb$hyrs
  nsau <- glb$nSAU
  saunames <- glb$sauname
  histC <- condC$histCatch
  filen <- filenametopath(rundir,"Historical_Catches.png")  #
  plothistcatch(yrs=hyrnames,histC,nsau,saunames,commonscale=commonscale,
                proportion=proportion,filen=filen)
  caption <- "Historical catches for each SAU, and across the zone."
  addplot(filen,rundir=rundir,category="Fishery",caption)
  filen <- filenametopath(rundir,"Historical_CPUE.png")
  plothistce(rundir,condC,glb,filen=filen)
  caption <- "Historical CPUE for each SAU."
  addplot(filen,rundir=rundir,category="Fishery",caption)
} # end of historical plots

#' @title HSphaseplot makes a phase plot of predicted mature depletion vs H
#'
#' @description HSphaseplot generates, for a single sau, a Kobe or phase plot
#'     of the predicted harvest rate against the predicted mature biomass
#'     depletion. It highlights the startyr of cpue conditioning (1992 in
#'     Tasmania) the end of conditioning (usually 2020 in Tasmania), and the
#'     final year (2050 if projected for 30 years). The default limit on harvest
#'     rate is set at 0.175, this may not be precautionary for all sau.
#'
#' @param dyn a list of the dynamics for every sau. This function only requires
#'     the mature biomass depletion and the harvest rates
#' @param glb the globals object
#' @param sau  which sau index (1 - 8 in Tasmania) to use?
#' @param rundir the directory for all files and results, default=""
#' @param startyr default = 1992. Otherwise could use condC$yearCE[1]
#' @param console should the plot go to the console or be saved? Default=TRUE
#' @param setpar if this is a stand-alone plot setting this TRUE will define a
#'     plot size and set the par values. If the plot is to be part of a
#'     multiple plot then setpar should be FALSE. default=TRUE
#' @param targdepl what target depletion is wanted, default=0.4
#' @param limdepl what is the limit depletion reference point, default = 0.2
#' @param limH what is the maximum H of limit fishing mortality, default = 0.2
#' #'
#' @return a list of the median mature biomass depletion and median harvest rate
#'     by year invisibly, it also plots these vectors.
#' @export
#'
#' @examples
#' \dontrun{
#'  print("wait on a data set")
#'  for (plotsau in 1:glb$nSAU) {
#'    HSphaseplot(dyn=sauout,glb=glb,sau=plotsau,
#'                rundir=rundir,startyr=condC$yearCE[1],
#'                console=FALSE,targdepl=kobeRP[1],
#'                limdepl=kobeRP[2],limH=kobeRP[3])
#'  }
#' }
HSphaseplot <- function(dyn,glb,sau,rundir="",startyr=1992,console=TRUE,
                        setpar=TRUE,targdepl=0.4,limdepl=0.2,limH=0.175) {
  # dyn <- sauout; maxdepl=0.6; startyr=condC$yearCE[1]; glb=glb; console=TRUE
  # targdepl <- kobeRP[1]; limdepl=kobeRP[2]; limH <-kobeRP[3]; sau=6; outhcr=outhcr
  yrs <- c(glb$hyrnames,glb$pyrnames)
  picky <- which(yrs >= startyr)
  condy <- which(yrs[picky] == glb$hyrnames[length(glb$hyrnames)])
  nyrs <- length(picky)
  depl <- dyn$deplsB[picky,sau,]
  H <- dyn$harvestR[picky,sau,]
  filen <- ""
  if ((!console) & (setpar)) {
    nfile <- paste0("PhasePlot_for_",glb$saunames[sau],".png")
    filen <- filenametopath(rundir,nfile)
    caption <- paste0("PhasePlot_of Harvest Rate vs Mature Biomass for ",
                      glb$saunames[sau])
  }
  medH <- apply(H,1,median)
  meddepl <- apply(depl,1,median)
  if (setpar) {
    plotprep(width=7,height=6,filename=filen,cex=0.9,verbose=FALSE)
    parset(cex.lab=1.25)
  }
  maxH <- getmax(medH)
  maxdepl <- getmax(meddepl,mult=1.25)
  plot(meddepl,medH,type="p",pch=16,xlim=c(0,maxdepl),ylim=c(0,maxH),xaxs="i",yaxs="i",
       xlab=paste0("Mature Biomass Depletion ",glb$saunames[sau]),ylab="Harvest Rate")
  phaseplotpolygons(maxy=maxH,maxx=maxdepl,targx=targdepl,limx=limdepl,limy=limH)
  suppressWarnings(arrows(x0=meddepl[1:(nyrs-1)],y0=medH[1:(nyrs-1)],
                          x1=meddepl[2:nyrs],y1=medH[2:nyrs],lwd=2,
                          length=0.075,code=2))
  # when H values repeatedly zero one gets zero length arrows and a warning
  # this is not a problem so it is suppressed.
  points(meddepl[c(1,condy,nyrs)],medH[c(1,condy,nyrs)],pch=16,cex=3,col=c(1,6,4))
  legend("topright",legend=c(startyr,yrs[picky[condy]],yrs[picky[length(picky)]]),
         col=c(1,6,4),lwd=5,cex=1.5,bty="n")
  if ((!console) & (setpar)) {
    addplot(filen,rundir=rundir,category="phaseplot",caption)
  }
  return(invisible(c(mediandepl=meddepl[nyrs],medianH=medH[nyrs])))
} # end of HSphaseplot

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
  parset(plots=pickbound(length(label)),margin=c(0.3,0.4,0.05,0.05),outmargin=c(1,0,1,0),
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
    CI <- apply(postmat[,],1,quantile,probs=CIprobs,na.rm=TRUE)
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
  CI <- apply(invar,1,quantile,probs=CIprobs,na.rm=TRUE)
  plot(yrnames[pickyr],invar[pickyr,1],type="l",lwd=1,col="grey",ylab=varname,
       xlab="Years",panel.first=grid(),ylim=c(0,ymax),yaxs="i")
  for (i in 1:reps) lines(yrnames[pickyr],invar[pickyr,i],col="grey",lwd=1)
  lines(yrnames[pickyr],CI[1,pickyr],lwd=1,col=2)
  lines(yrnames[pickyr],CI[3,pickyr],lwd=1,col=2)
  lines(yrnames[pickyr],CI[2,pickyr],lwd=2,col=4)
  caption <- paste0("Projected ",varname," across all replicates.")
  if (addfile) addplot(filen,rundir=rundir,category="zonescale",caption)
  return(invisible(CI))
} # end of onezoneplot

#' @title phaseplotpolygons sets up the colour blocks for phase plots
#'
#' @description phaseplotpolygons is used to put together coloured phase plots
#'     a common example of which would be kobe plots. This function assumes that
#'     a plot has already been defined as well as the maximum x and y values.
#'     In a standard kobe plot y would relate to fishing mortality or harvest
#'     rate, and x would relate to biomass or depletion. If plotting a kobe plot
#'     then the target and limit depletion are important as is the limit harvest
#'     rate. This function only generates the coloured blocks it is for use
#'     within other functions
#'
#' @param maxy the upper bound on the y-axis
#' @param maxx the upper bound on the x-axis
#' @param targx the target reference point depletion or biomass
#' @param limx the limit reference point depletion or biomass
#' @param limy the limit reference point for fishing mortality or harvest rate
#'
#' @return nothing but it does add to a plot
#' @export
#'
#' @examples
#' \dontrun{
#' # illustrating a form of use
#'  medH <- apply(H,1,median)
#'  meddepl <- apply(depl,1,median)
#'  plotprep(width=7,height=6,filename=filen,cex=0.9,verbose=FALSE)
#'  parset(cex.lab=1.25)
#'  maxH <- getmax(medH)
#'  maxdepl <- getmax(meddepl,mult=1.25)
#'  plot(meddepl,medH,type="p",pch=16,xlim=c(0,maxdepl),ylim=c(0,maxH),xaxs="i",
#'       yaxs="i",xlab=paste0("Mature Biomass Depletion ",glb$saunames[sau]),
#'       ylab="Harvest Rate")
#'  phaseplotpolygons(maxH=maxH,maxdepl=maxdepl,targdepl=targdepl,
#'                    limdepl=limdepl,limH=limH)
#' }
phaseplotpolygons <- function(maxy,maxx=1.0,targx=0.4,limx=0.2,limy=0.175) {
  polygon(x=c(0,limx,limx,0,0),y=c(0,0,maxy,maxy,0),
          col=rgb(255,0,0,175,maxColorValue=255))
  polygon(x=c(limx,maxx,maxx,limx,limx),
          y=c(maxy,maxy,limy,limy,maxy),
          col=rgb(255,255,0,140,maxColorValue=255))
  polygon(x=c(limx,targx,targx,limx,limx),y=c(limy,limy,0,0,limy),
          col=rgb(190,255,0,130,maxColorValue=255))
  polygon(x=c(targx,maxx,maxx,targx,targx),y=c(0,0,limy,limy,0),
          col=rgb(0,255,0,120,maxColorValue=255))
  abline(h=limy,lwd=2,col=1)
  abline(v=targx,lwd=2,col=1)
} # end of phaseplotpolygons

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
  parset(plots=pickbound(nsau),margin = c(0.25, 0.45, 0.05, 0.05),byrow=FALSE,
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
#' @param histCE the historical CPUE matrix for all sau
#' @param histCatch the historical catch matrix for all sau
#' @param rundir the rundir for the given scenario
#' @param recdevs the recruitment deviates from the control file, condC
#' @param console should the plot go to the console or a file, default=TRUE
#'
#' @return the sau scale zone dynamics, invisibly
#' @export
#'
#' @examples
#' print("wait on suitable data sets")
plotconditioning <- function(zoneDD,glb,zoneC,histCE,histCatch,rundir,
                             recdevs,console=TRUE) {
  # zoneDD=zoneDD;glb=glb;zoneC=zoneC;histCE=condC$histCE;rundir=rundir;recdevs=condC$recdevs
  sauindex <- glb$sauindex
  popB0 <- getlistvar(zoneC,"B0")
  B0 <- tapply(popB0,sauindex,sum)
  popExB0 <- getlistvar(zoneC,"ExB0")
  ExB0 <- tapply(popExB0,sauindex,sum)
  sauZone <- getsauzone(zoneDD,glb,B0=B0,ExB0=ExB0)
  allsau <- 1:glb$nSAU #
  snames <- glb$saunames
  yrnames <- glb$hyrnames
  filen <- ""
  if (!console) filen <- "compareCPUE.png"
  ssq <- compareCPUE(histCE,sauZone$cpue,glb,rundir,filen=filen)
  for (sau in allsau) {  #  sau=allsau[1]; filen="";  sau=allsau[1]
    filen="" #filen <- filenametopath(rundir,paste0(snames[sau],"_conditioned.png"))
    if (!console) {
      filen <- filenametopath(rundir,paste0(snames[sau],"_conditioned.png"))
      caption <- "Dynamics during conditioning on the fishery."
    }
    plotprep(width=9,height=9,newdev=TRUE,filename=filen,cex=0.9,verbose=FALSE)
    dosau(sauZone,glb,picksau=sau,histCE=histCE,histCatch=histCatch,
          yrnames=glb$hyrnames,recdevs[,sau])
    if (!console) addplot(filen,rundir=rundir,category="condition",caption)
  }
  return(invisible(list(sauZone=sauZone,ssq=ssq)))
} # end plotconditioning.

#' @title plothistcatch generates a plot of the historical catches by SAU
#'
#' @description plothistcatch generates a plot of histircal catches by SAU
#'     and by Zone. It can plot can tonnes x year or proportion. All SAU
#'     plots can have independent or common y-axis scales.
#'
#' @param yrs a vector of the years of catch history
#' @param histC the matrix of catches by years
#' @param nsau the number of SAU
#' @param saunames the SAU names as a vector
#' @param commonscale should the SAU plots share a common y-axis scale,
#'     default =  FALSE
#' @param proportion should the plot be of proportion relative to the maximum
#'     catch through each time-series or as raw tonnes x year, default=FALSE
#' @param filen the name of the file used to save the plot, default="", which
#'     sends the plot to a separate window
#'
#' @return nothing but it does generate a plot either to a window or a file
#' @export
#'
#' @examples
#' print("wait on suitable data sets")
plothistcatch <- function(yrs,histC,nsau,saunames,commonscale=FALSE,
                          proportion=FALSE,filen="") {
  totC <- rowSums(histC,na.rm=TRUE)
  if (proportion) {
    for (i in 1:nsau) {
      maxy <- max(histC[,i])
      histC[,i] <- histC[,i]/maxy
    }
    commonscale=TRUE
    totC <- totC/max(totC)
  }
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  parset(plots=pickbound(nsau+1),margin=c(0.3,0.35,0.05,0.05),
         outmargin=c(1,1,0,0),byrow=FALSE)
  if (commonscale) ymax <- getmax(histC)
  for (i in 1:nsau) {
    if (!commonscale) ymax <- getmax(histC[,i])
    plot(yrs,histC[,i],type="l",lwd=2,xlab="",ylab=saunames[i],yaxs="i",
         ylim=c(0,ymax),panel.first=grid())
  }
  ymax <- getmax(totC)
  plot(yrs,totC,type="l",lwd=2,xlab="",ylab="Zone",yaxs="i",
       ylim=c(0,ymax),panel.first=grid())
  mtext("Years",side=1,line=-0.1,outer=TRUE,cex=1.1)
  label <- "Annual Catch (t)"
  if (!commonscale) label <- paste0(label," (note different y-axis scales)")
  mtext(label,side=2,line=-0.1,outer=TRUE,cex=1.1)
} # end of plothistcatch

#' @title plothistce generates a plot of the historical cpue by SAU
#'
#' @description plothistce generates a plot of historical cpue by SAU. It can
#'     plot as cpue or proportion of maximum cpue.
#'
#' @param rundir the direcotry where all results are stored
#' @param condC the R object contiaining the historical fishery data generated
#'     by the readctrlfile function
#' @param glb the globol variables object
#' @param proportion should the plot be of proportion relative to the maximum
#'     catch through each time-series or as raw tonnes x year, default=FALSE
#' @param filen the name of the file used to save the plot, default="", which
#'     sends the plot to a separate window
#'
#' @return nothing but it does generate a plot either to a window or a file
#' @export
#'
#' @examples
#' print("wait on suitable data sets")
plothistce <- function(rundir,condC,glb,proportion=FALSE,filen="") {
  nsau <- glb$nSAU
  saunames <- glb$saunames
  yrs <- condC$yearCE
  histCE <- condC$histCE
  if (proportion) {
    for (sau in 1:nsau) {
      maxy <- max(histCE[,sau],na.rm=TRUE)
      histCE[,sau] <- histCE[,sau]/maxy
    }
  }
  plotprep(width=7,height=8,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  parset(plots=pickbound(nsau),margin=c(0.3,0.35,0.05,0.05),
         outmargin=c(1,1,0,0),byrow=FALSE)
  for (sau in 1:nsau) {
    ymax <- getmax(histCE[,sau])
    plot(yrs,histCE[,sau],type="l",lwd=2,ylim=c(0,ymax),yaxs="i",xlab="",
         ylab=saunames[sau],panel.first=grid())
  }
  mtext("Years",side=1,line=-0.1,outer=TRUE,cex=1.1)
  label <- "Annual CPUE"
  if (proportion) label <- paste0(label," as proportion of maximum")
  mtext(label,side=2,line=-0.1,outer=TRUE,cex=1.1)
} # end of plothistce

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
#' @param average should the histograms be of the average catch per year,
#'     default=FALSE, so it plots the actual sum of catches
#'
#' @return nothing but it does generate two plots into rundir
#' @export
#'
#' @examples
#' print("wait on suitable in ternal data sets")
plothsstats <- function(rundir,hsstats,glb,average=FALSE) {
  # hsstats=out$HSstats; rundir=rundir; glb=out$glb
  sum5 <- hsstats$sum5
  sum10 <- hsstats$sum10
  labelnames <- colnames(sum10)
  nsau <- glb$nSAU
  if (average) {
    sum5 <- sum5/5.0
    sum10 <- sum10/10.0
    pngname5 <- "Mean_Annual_Catch_5yr.png"
    pngname10 <- "Mean_Annual_Catch_10yr.png"
    captionadd <- "mean catch across "
  } else {
    pngname5 <- "sum_5yr_proj_catches.png"
    pngname10 <- "sum_10yr_proj_catches.png"
    captionadd <- "sum of catches across "
  }
  filen <- filenametopath(rundir,pngname5)
 # filen=""
  plotprep(width=8,height=7,filename=filen,cex=1.0,verbose=FALSE)
  parset(plots=pickbound(ncol(sum5)))
  for (i in 1:(nsau+1)) {
    label <- paste0("Tonnes     ",labelnames[i])
    hist(sum5[,i],breaks=20,main="",xlab=label)
  }
  caption <- paste0("The ",captionadd,"the first 5 years of he projections",
                    " for each SAU and the complete zone")
  addplot(filen,rundir=rundir,category="HSperf",caption)
  #now 10 years
  filen <- filenametopath(rundir,pngname10)
  #filen=""
  plotprep(width=8,height=7,filename=filen,cex=1.0,verbose=FALSE)
  parset(plots=pickbound(ncol(sum10)))
  for (i in 1:(nsau+1)) {
    label <- paste0("Tonnes     ",labelnames[i])
    hist(sum10[,i],breaks=20,main="",xlab=label)
  }
  caption <- paste0("The ",captionadd,"the first 10 years of he projections",
                    " for each SAU and the complete zone")
  addplot(filen,rundir=rundir,category="HSperf",caption)
} # end of plothsstats

#' @title plotNt plots the size-composition for each SAU
#'
#' @description plotNt accepts size-composition data (either the population Nt,
#'     of the catchN), and plots the replicate size-distributions for each SAU.
#'     It then puts the first year's size-distribution on top and if the medcol
#'     that defines the colour of the median of the replicates is different
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
  parset(plots=pickbound(nsau),byrow=FALSE)
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

#' @title plotpopprops plots out some properties of input populations
#'
#' @description once the depleted zone (zoneCP and zoneDD) has been generated
#'     the component populations have their properties summarized in propD. This
#'     is saved as propertyDD.csv and output to the zoneDD tab. To provide a
#'     visual consideration of these properties this function generated
#'     histograms for an array of input array column names.
#'
#' @param x the array containing columns to plot
#' @param rundir the run directory for the results, can be set to "" if plotting
#'     to the console
#' @param glb the global constants object
#' @param varnames the names of each variable to select the columns and to use
#'     in labels
#' @param startyr in what year number should the plot refer to?
#' @param console  should plot go to the console or a png be saved? Default=TRUE
#' @param bins the number of breaks to plot in each histogram, default=25
#'
#' @return it generates a png file of the plot if console = FALSE, otherwise it
#'     is drawn to the console
#' @export
#'
#' @examples
#' print("wait on suitable data")
#' # x=propD;rundir=rundir;glb=glb;varnames=columns;startyr=hyrs;bins=21;console=FALSE
plotpopprops <- function(x,rundir,glb,varnames,startyr,console=TRUE,bins=25) {
  npop <- glb$numpop
  n <- length(varnames)
  if (console) { filen="" } else {
    filen <- filenametopath(rundir,"population_properties.png")
  }
  plotprep(width=8,height=7,newdev=FALSE,filename=filen,verbose=FALSE)
  parset(pickbound(n),byrow=FALSE,outmargin=c(0,1,0,0),
         margin=c(0.4,0.3,0.05,0.05))
  for (i in 1:n) {
    xcol <- x[1:npop,varnames[i]]
    hist(xcol,main="",xlab=varnames[i],ylab="",breaks=bins)
  }
  mtext("Frequency",side=2,outer=TRUE,line=-0.2,cex=1.2)
  if (!console) {
    caption <- paste0(" Selected Biological properties of all popualtions",
                      " some are for year ",startyr,"; from propertyDD.csv.")
    addplot(filen,rundir=rundir,category="zoneDD",caption)
  }
} # end of plotpopprops



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

#' @title plotsizecomp generates a plot of available size composition data
#'
#' @description plotsizecomp generates a plot of available size composition
#'     data with the option of plotting the predicted size-composition of
#'     catch on top of it; this latter option is only possible when plotting
#'     the distributions as proportions.
#'
#' @param rundir the directory for all ctrl, data, and output files.
#' @param incomp filtered size-composition data from the fishery
#' @param SAU the name of the SAU used to label the y-axis
#' @param lml a vector of the lml in each year represented by incomp
#' @param catchN the predicted size-composition of the catch
#' @param start which size-class to start from, default=NA, which means all
#'     size-classes in the samples will be used.
#' @param proportion should the plots be as proportions or counts,default=TRUE
#' @param console should the plot be sent ot the console or saved to rundir.
#'     default=TRUE
#' @param width the width of the plot, default = 10
#' @param height the height of the plot, default = 9
#' @param fnt what font to use, default font = 7, bold times
#' @param tabcategory what name to give to the webage tab when the plots are
#'     saved. Default="predictedcatchN"
#'
#' @seealso{
#'  \link{preparesizecomp}, \link{popNAStosau}
#' }
#'
#' @return nothing but it does generate a plot
#' @export
#'
#' @examples
#' print("wait on data sets")
plotsizecomp <- function(rundir,incomp,SAU="",lml=NULL,catchN=NULL,start=NA,
                         proportion=TRUE,console=TRUE,width=10,height=9,
                         fnt=7,tabcategory="predictedcatchN") {
  # rundir=rundir;incomp=lfs;SAU=glb$saunames[plotsau];lml=LML;width=10;height=9
  # catchN=sauCt[,,plotsau];start=NA;proportion=TRUE;console=TRUE;fnt=7
  if (console) { filen <- "" } else {
    filen <- filenametopath(rundir,paste0(SAU,"sizecomp_nas_by_year.png"))
  }
  mids <- as.numeric(rownames(incomp))
  nmids <- length(mids)
  if (is.numeric(start)) {
    mids <- mids[start:nmids]
    nmids <- length(mids)
  }
  yrsize <- as.numeric(colnames(incomp))
  nyrs <- length(yrsize)
  nobs <- colSums(incomp,na.rm=TRUE)
  if (proportion) {
    incomp <- prop.table(incomp,margin=2)
    if (!is.null(catchN)) {
      years <- as.numeric(colnames(catchN))
      midpts <- as.numeric(rownames(catchN))
      picky <- match(yrsize,years)
      picks <- match(mids,midpts)
      pNt <- prop.table(catchN[picks,picky],margin=2)
    }
  }
  plotprep(width=width,height=height,newdev=FALSE,filename=filen,
           verbose=FALSE,usefont=fnt)
  parset(plots=pickbound(nyrs),margin=c(0.225,0.225,0.05,0.05),
         outmargin=c(1,1,0,0),byrow=FALSE)
  for (yr in 1:nyrs) { # yr = 3
    y <- incomp[,yr]
    ymax <- getmax(y)
    if ((proportion) & (!is.null(catchN))) ymax <- getmax(c(y,pNt[,yr]))
    plot(mids,y,type="l",lwd=2,ylab="",ylim=c(0,ymax),yaxs="i",
         panel.first = grid(),xlab="")
    mtext(text=yrsize[yr],side=4,line=-1.1,cex=1.0,font=7)
    text(mids[trunc(nmids*0.65)],0.9*ymax,nobs[yr],cex=1.25,pos=4)
    if ((proportion) & (!is.null(catchN))) {
      lines(mids,pNt[,yr],lwd=2,col=2)
    }
    if (!is.null(lml)) {
      if (ncol(lml) > 2) {
        pickLML <- grep(SAU,colnames(lml))
        abline(v=lml[yr,pickLML],lwd=2,col="blue")
      } else {
        abline(v=lml[yr,2],lwd=2,col="blue")
      }
    }
  }
  mtext("Shell Length (mm)",side=1,line=-0.2,outer=TRUE,cex=1.2)
  mtext(paste0("Proportion for ",SAU),side=2,line=-0.2,outer=TRUE,cex=1.2)
  if (!console) {
    caption <- paste0("Catch size-composition data to which the SBM is being ",
                      "fitted. The vertical blue lines are the LML by year.")
    addplot(filen,rundir,category=tabcategory,caption=caption)
  }
} # end of plotsizecomp



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
#' @return the median and 90th quantiles (invisibly) and it adds a plot to the
#'     zonescale tab in the rundir
#' @export
#'
#' @examples
#' print("wait on suitable datasets")
#' # inzone=outzone;rundir=rundir;glb=glb;CIprobs=c(0.05,0.5,0.95);addfile=TRUE;
#' # startyr=startyr; addfile=TRUE
plotZone <- function(inzone,rundir,glb,startyr,CIprobs=c(0.05,0.5,0.95),addfile=TRUE) {
  catchCI <- onezoneplot(inzone$catch,rundir,glb,CIprobs=CIprobs,"Catch",
                         startyr=startyr,addfile=addfile)
  cpueCI <- onezoneplot(inzone$cpue,rundir,glb,CIprobs=CIprobs,"CPUE",
                        startyr=startyr,addfile=addfile)
  deplCI <- onezoneplot(inzone$deplsB,rundir,glb,CIprobs=CIprobs,
                        "Mature_Biomass_Depletion",startyr=startyr,addfile=addfile)
  matbCI <- onezoneplot(inzone$matureB,rundir,glb,CIprobs=CIprobs,"Mature_Biomass",
                        startyr=startyr,addfile=addfile)
  HCI <- onezoneplot(inzone$harvestR,rundir,glb,CIprobs=CIprobs,"Harvest_Rate",
                     startyr=startyr,addfile=addfile)
  recCI <- onezoneplot(inzone$recruit,rundir,glb,CIprobs=CIprobs,"Recruitment",
                       startyr=startyr,addfile=addfile)
  tacCI <- onezoneplot(inzone$TAC,rundir,glb,CIprobs=CIprobs,"TAC",startyr=startyr,
                       addfile=addfile)
  meds <- cbind(catchCI[2,],cpueCI[2,],deplCI[2,],matbCI[2,],HCI[2,],recCI[2,],
                tacCI[2,])
  colnames(meds) <- c("Catch","cpue","depl","matureB","HarvestR","Recruit","TAC")
  rownames(meds) <- c(glb$hyrnames,glb$pyrnames)
  CI <- cbind(catchCI[1,],catchCI[3,],cpueCI[1,],cpueCI[3,],deplCI[1,],deplCI[3,],
              matbCI[1,],matbCI[3,],HCI[1,],HCI[3,],recCI[1,],recCI[3,],
              tacCI[1,],tacCI[3,])
  rownames(CI) <- c(glb$hyrnames,glb$pyrnames)
  colnames(CI) <- c("CatchL90","CatchU90","cpueL90","cpueU90","deplL90","deplU90",
                    "matBL90","matBU90","HL90","HU90","recL90","recU90",
                    "tacL90","tacU90")
  addtable(meds,"zone-medians.csv",rundir,category="zonescale",
           caption="Medians of zonescale dynamics.")
  return(invisible(list(zonemeds=meds,zoneCI=CI)))
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
#' @seealso{
#'   \link{catchweightCE}
#' }
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

#' @title poplevelplot generates population x replicate trajectory plots
#'
#' @description poplevelplot is used to produce population level plots of an
#'     input variable derived from zoneDP. It produces a single plot for a given
#'     sau which depicts each population's replicate trajectories, all in grey.
#'     However, medians of the variable through each year are calculated and
#'     plotted in distinct colours for each population up to a maximum of 12
#'     populations per sau. If more than 12 then the population index 1:npop is
#'     used. These plots should be considered as a partial diagnostic for how
#'     well the implied fleet dynamics is distributing catch among an sau's
#'     populations. The trajectories are plotted from the last year of
#'     conditioning to the last year of projections.
#'
#' @param rundir the scenario's directory for datafiles and results
#' @param sau which sau is to be worked on. In Tas = 1:8 in the west
#' @param popvar which 3D array population variable is to be used. This will be
#'     generally derived from zoneDP, eg zoneDP$catch or zoneDP$depleB
#' @param glb the globals object used for sauindex, saunames, hyrs, and reps
#' @param label the sau name will be added to this label for the y-axis, this
#'     should be simple, such as 'catch' or have words connected by hyphens or
#'     underscores as it will also be added to the filename.
#' @param console should the plot go to the console or be saved into rundir.
#'     default=TRUE = to the console
#'
#' @return invisibly the population medians
#' @export
#'
#' @examples
#' print("poplevelplot(rundir=rundir,sau=7,popvar=zoneDP$catch,glb=glb,")
#' print("label='Catch',console=TRUE")
poplevelplot <- function(rundir,sau,popvar,glb,label="",console=TRUE) {
  popcol <- c("red","#FF7790","orange","yellow","#80FE90","green","cyan",
              "#0080FF","blue","#8000FF","magenta","#FF0080")
  sauindex <- glb$sauindex
  saunames <- glb$saunames
  nyrs <- length(dimnames(popvar)[[1]])
  saudepl <- popvar[glb$hyrs:nyrs,(sauindex == sau),]
  popnums <-  as.numeric(dimnames(saudepl)[[2]])
  npop <- length(popnums)
  reps <- glb$reps
  if (console) {
    filen <- ""
  } else {
    nfile <- paste0("Pop_Trajs_and_Medians_for_",label,"_",saunames[sau],".png")
    filen <- pathtopath(rundir,nfile)
    caption <- paste0(saunames[sau]," trajectories for each population ",
                      "x replicates of ",label)
  }
  yrs <- as.numeric(dimnames(saudepl)[[1]])
  numyrs <- length(yrs)
  maxy <- getmax(saudepl)
  label <- paste0(saunames[sau],"  ",label)
  plotprep(width=9,height=4.5,newdev=TRUE,filename=filen,cex=0.9,verbose=FALSE)
  parset()
  plot(yrs,saudepl[,1,1],type="l",lwd=1,col="grey",ylim=c(0,maxy),
       panel.first=grid(),xlab="",ylab=label)
  for (iter in 1:reps) {
    for (pop in 1:npop) {
      lines(yrs,saudepl[,pop,iter],lwd=1,col="grey")
    }
  }
  meds <- matrix(0,nrow=numyrs,ncol=npop,dimnames=list(yrs,popnums))
  if (npop > 12) {
    outcol <- 1:npop
  } else {
    outcol <- popcol
  }
  for (pop in 1:npop) {
    meds[,pop] <- apply(saudepl[,pop,],1,median)

    lines(yrs,meds[,pop],lwd=3,col=outcol[pop])
  }
  legend("topright",legend=popnums,col=c(popcol[1:npop]),lwd=4,bty="n",cex=1)
  if (!console) {
    addplot(filen,rundir=rundir,category="poplevelplots",caption)
  }
  return(invisible(meds))
} # end of poplevelplot

#' @title predsizecomp plots size-composition for each replicate in a subset of years
#'
#' @description predsizecomp generates plots the size composition data for a
#'     subset of the projection years (determined by the interval argument).
#'     All plots omit the first five size classes so that the post-settlement
#'     spikes do not obscure the rest of the plots. The minSL (minimum shell
#'     length argument) determines the left-most size in the plots to allow for
#'     details to be made clearer.
#'
#' @param sau which sau to plot, the index number not its name, in TAS western
#'     zone this would be between 1 : 8
#' @param NSC the numbers-at-size matrix from sauout, either Nt or catchN
#' @param glb the globals object
#' @param minSL the minimum shell length to be plotted. This constrains each
#'     plot so that details can be clearer. default=10mm which implies all
#'     sizes above the post-larval sizes.
#' @param interval default=5 which selects the year just before projections and
#'     then the data for every interval years along the projection period. The
#'     years selected always includes the last year of projection and the year
#'     before projections started.
#' @param prop default=TRUE. Should proportions be plotted or actual numbers
#' @param console should the plot go to the console or will it be saved
#' @param rundir the directory into which all results are placed. default="" but
#'     generally this would = rundir
#' @param omit default = 5. how many size classes should be omitted from the
#'     start of each compositions plot. The default = 5 which in Tasmania,
#'     which has 2mm size classes starting at 2mm, this omits all animals less
#'     than 10mm.
#'
#' @return the outcome which is OK if all replicates OK, but names the last year
#'     group to have problems in its replicates if any hit zero catch, and it
#'     plots a graph.
#' @export
#'
#' @examples
#' print("wait on suitable data")
#' # sau=1;NSC=sauNAS$catchN; glb=glb; minSL=minsizecomp[2];omit=5
#' # interval=nasInterval;prop=TRUE;console=TRUE;rundir=rundir
#' # sau=sau; NSC=sauNAS$Nt; glb=glb; minSL=minsizecomp[1];omit=5
#' # interval=nasInterval;prop=TRUE;console=FALSE;rundir=rundir;omit=5
predsizecomp <- function(sau,NSC,glb,minSL=10,interval=5,prop=TRUE,
                         console=TRUE,rundir="",omit=5) {
  outcome <- "OK"
  Nclass <- glb$Nclass
  pickSC <- (omit+1):Nclass
  if (minSL > 0) {
    pickS <- which(glb$midpts >= minSL)
    pickSC <- pickS[1]:Nclass  # which size classses
  }
  pickyrs <- seq(glb$hyrs,(glb$hyrs+glb$pyrs),interval)  # which years
  yrnames <- c(glb$hyrnames,glb$pyrnames)
  picknames <- yrnames[pickyrs]
  numyrs <- length(pickyrs)
  medyrs <- matrix(0,nrow=Nclass,ncol=numyrs,dimnames=list(glb$midpts,picknames))
  sauNt <- NSC[,pickyrs,sau,]
  reps <- glb$reps
  shlength <- as.numeric(dimnames(sauNt)[[1]])
  filen <- ""
  if (!console) {
    if (sum(sauNt[1:40,1,1]) > 1) label="Stock_" else label = "Catch_"
    nfile <- paste0(label,"Size_Composition_for_",glb$saunames[sau],".png")
    filen <- pathtopath(rundir,nfile)
    caption <- paste0(label,"size-composition for ",glb$saunames[sau])
  }
  plotprep(width=8,height=9,filename=filen,cex=0.9,verbose=FALSE)
  parset(plots=pickbound(numyrs+1),margin=c(0.25,0.45,0.1,0.1),byrow=FALSE,
         outmargin=c(1.5,1,0,0))
  for (i in 1:numyrs) {  # i=3
      yrsauNt <- psauNt <- sauNt[,i,]
      if (prop) {
        reptotal <- colSums(yrsauNt,na.rm=TRUE); reptotal
        pickrep <- which(reptotal > 0)
        if (length(pickrep) < length(reptotal)) {
          reptotal[-pickrep] <- 1
          nreps <- length(reptotal) - length(pickrep)
          issue <- paste0("In predsizecomp ", nreps," Replicates in ",
                            picknames[i]," for sau ",sau," were zero \n")
          addwarning(glb$warnfile,issue)
        }
        for (rp in 1:reps) psauNt[,rp] <- yrsauNt[,rp]/reptotal[rp]
        reptotal <- colSums(psauNt,na.rm=TRUE)
        pickpos <- which(reptotal > 0)
        lenpickpos <- length(pickpos)
      }
      if (lenpickpos > 1) {
         qts <- apply(psauNt[,pickpos],1,quants)
         medyrs[,i] <- qts[3,]
      } else {
        if (lenpickpos == 1) {
           medyrs[,i] <- psauNt[,pickpos]
        } else {
          medyrs[,i] <- 0
        }
      }
      if (lenpickpos > 0) {
        maxy <- max(psauNt[pickSC,pickpos])*1.025
        plot(shlength[pickSC],psauNt[pickSC,pickpos[1]],type="l",lwd=1,col="grey",
             xlab="",ylim=c(0,maxy),ylab=picknames[i],panel.first=grid())
      }
      if (lenpickpos > 1) {
        for (j in 2:lenpickpos) lines(shlength[pickSC],psauNt[pickSC,pickpos[j]],
                                    lwd=1,col="grey")
        lines(shlength[pickSC],qts[3,pickSC],lwd=2,col="red")
      }
  }
  maxy <- getmax(medyrs[pickSC,])
  plot(shlength[pickSC],medyrs[pickSC,1],type="l",lwd=2,col=1,xlab="",
       ylim=c(0,maxy),ylab="Medians",panel.first=grid())
  for (j in 2:numyrs) lines(shlength[pickSC],medyrs[pickSC,j],lwd=2,col=j)
  legend("topright",legend=picknames,col=1:numyrs,lwd=3,bty="n",cex=0.9)
  mtext(paste0("Shell Length (mm)      ",glb$saunames[sau]),side=1,line=-0.1,
        outer=TRUE,cex=1.2)
  label <- "Frequency"
  if (prop) label <- "Proportion"
  mtext(label,side=2,line=-0.2,outer=TRUE,cex=1.2)
  if (!console) {
    addplot(filen,rundir=rundir,category="NumSize",caption)
  }
 return(invisible(outcome))
} # end of predsizecomp

#' @title preparesizecomp strips out empty columns and identifies samples
#'
#' @description preparesizecomp takes in the size composition of catch data
#'     and removes empty columns and those with totals of less than the
#'     defined mincount. The minimum used depends on sampling routines. In
#'     Tasmania, a sample of 100 is the usual minimum taken from individual
#'     landings. Single landings may not be representative of much.
#'
#' @param sizecomp the size-composition of catch data from the dat object
#'     made by readLBMdata.
#' @param mincount the minimum number of observations within a year for
#'     inclusion in the analysis, default = 100.
#' @param deleteyears default = 0, which means all years will be retained.
#'     Otherwise, this can be a vector of years, all of which will be removed
#'     before the mincount is applied.
#'
#' @return a cleaned version of the sizecomp matrix
#' @export
#'
#' @examples
#' print("wait on data sets")
#' # sizecomp=compdat[,,sau]; mincount=100; deleteyears=0
preparesizecomp <- function(sizecomp,mincount=100,deleteyears=0) {
  if (deleteyears[1] > 0) {
    pickcols <- which(as.numeric(colnames(sizecomp)) %in% deleteyears)
    if (length(pickcols) > 0) sizecomp <- sizecomp[,-pickcols]
  }
  mids <- as.numeric(rownames(sizecomp))   #sizecomp[,"length"]
  totals <- colSums(sizecomp,na.rm=TRUE)
  pick <- which(totals > mincount)
  sizecomp <- sizecomp[,pick]
  # numcol <- ncol(sizecomp)
  # sizecomp <- sizecomp[,2:numcol] # ignore the length column
  return(sizecomp)
} # end of preparesizecomp

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
#' @param hlines default = NULL, otherwise a dataframe of vectors of sau
#'     properties. eg sauprod containing Bmsy, MSY, Dmsy, CEmsy for adding to
#'     each projSAU plot
#'
#' @return a list of lists of CI for each SAU and variable as well as the
#'     zoneDsau and zonePsau
#' @export
#'
#' @examples
#' print("wait on suitable internal data-sets")
sauplots <- function(zoneDP,NAS,glb,rundir,B0,ExB0,startyr,addCI=TRUE,
                     histCE=NULL,tabcat="projSAU",hlines=NULL) {
  # zoneDP=zoneDP;NAS=NAS;glb=glb;rundir=rundir;B0=B0;ExB0=ExB0;
  # startyr=48; addCI=TRUE;histCE=condC$histCE; tabcat="projSAU"; hlines=NULL
  zonePsau <- zonetosau(zoneDP,NAS,glb,B0,ExB0)
  label <-  c("cpue","catch","acatch","matureB","exploitB","recruit","harvestR")
  finalcondyeardepletion(rundir,sauzone=zonePsau,glb,deplvar="sB",console=FALSE)
  finalcondyeardepletion(rundir,sauzone=zonePsau,glb,deplvar="eB",console=FALSE)
  #CPUE
  filen <- filenametopath(rundir,"proj_cpue_SAU.png")  # filen=""
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  CI <- dosauplot("cpue",zonePsau[["cpue"]],glb,addCE=TRUE,startyr=startyr,
                  addCI=TRUE,histCE=histCE,hlin=hlines["CEmsy",])
  caption <- "The CPUE projections for each SAU."
  addplot(filen,rundir=rundir,category=tabcat,caption)
  #Catches
  filen <- filenametopath(rundir,"proj_catch_SAU.png")
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  CI <- dosauplot("catch",zonePsau[["catch"]],glb,
                  startyr=startyr,addCI=TRUE,histCE=NULL,hlin=hlines["MSY",])
  caption <- "The catch projections for each SAU."
  addplot(filen,rundir=rundir,category=tabcat,caption)
  #Aspirational catches
  filen <- filenametopath(rundir,"proj_aspcatch_SAU.png")
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  CI <- dosauplot("acatch",zonePsau[["acatch"]],glb,
                  startyr=startyr,addCI=TRUE,histCE=NULL,hlin=NULL)
  caption <- paste0("The Aspirational catch projections for each SAU. Catches ",
                    "prior to HS are actual catches.")
  addplot(filen,rundir=rundir,category=tabcat,caption)
  #MatureBiomass
  filen <- filenametopath(rundir,"proj_matureB_SAU.png")
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  CI <- dosauplot("matureB",zonePsau[["matureB"]],glb,startyr=startyr,
                  addCI=TRUE,histCE=NULL,hlin=hlines["Bmsy",])
  caption <- "The mature biomass projections for each SAU."
  addplot(filen,rundir=rundir,category=tabcat,caption)
  #exploitable biomass
  filen <- filenametopath(rundir,"proj_exploitB_SAU.png")
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  CI <- dosauplot("exploitB",zonePsau[["exploitB"]],glb,
                  startyr=startyr,addCI=TRUE,histCE=NULL,hlin=NULL)
  caption <- "The exploitable biomass projections for each SAU."
  addplot(filen,rundir=rundir,category=tabcat,caption)
  #recruitment
  filen <- filenametopath(rundir,"proj_recruit_SAU.png")
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  CI <- dosauplot("recruit",zonePsau[["recruit"]],glb,
                  startyr=startyr,addCI=TRUE,histCE=NULL,hlin=NULL)
  caption <- "The recruitment projections for each SAU."
  addplot(filen,rundir=rundir,category=tabcat,caption)
 # harvest rate
  filen <- filenametopath(rundir,"proj_harvestR_SAU.png")
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  CI <- dosauplot("harvestR",zonePsau[["harvestR"]],glb,
                  startyr=startyr,addCI=TRUE,histCE=NULL,hlin=hlines["Hmsy",])
  caption <- "The Harvest Rate projections for each SAU."
  addplot(filen,rundir=rundir,category=tabcat,caption)
  return(invisible(zonePsau))
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
#' # recdevs=condC$recdevs;glb=glb;rundir=rundir; filen=""
saurecdevs <- function(recdevs,glb,rundir,filen="") {
  if (nchar(filen) > 0) filen <- filenametopath(rundir,filen)
  years <- as.numeric(rownames(recdevs))
  arecdevs <- abs(recdevs)
  nsau <- glb$nSAU
  label <- glb$saunames
  doplots=pickbound(nsau)
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

#' @title tasphaseplot generates a Tasmanian proxy phase (kobe) plot
#'
#' @description tasphaseplot generates a Tasmanian proxy phase (kobe) plot. This
#'     means it uses the gradient 4 scores as a proxy for the harvest rate and
#'     the target CPUE score as a proxy for the biomass depletion level. The
#'     scores rank from 0 - 10. The target CPUE score behaves appropriately in
#'     that high scores = things are good low scores things are bad, but the
#'     gradient 4 score is the reverse. So, the scores are reversed 10 - score
#'     and the axis labels adjusted accordingly. Thus, original high g4 scores
#'     indicate lower harvest rates while low g4 scores indicate high harvest
#'     rates. 10 - g4s reverses this trend and to indicate the change we use
#'     a scale from 5 - minus 5 with a limit refernece point at 0.
#'
#' @param proxyB the target CPUE score from outhcr. This is a 3D matrix so the
#'     median across replicates is used.
#' @param proxyH the gradient 4 score from outhcr. This is a 3D matrix so the
#'     median across replicates is used. The plot is of 10 - median(proxyH)
#' @param glb the globals object - just for the saunames
#' @param sau which sau is being plotted. For the western zone this is 1 - 8
#' @param rundir the rundir for the scenario where the plot will be stored if
#'     console = FALSE
#' @param console plot to the console or save the plot. default=TRUE
#' @param fnt what font to use, default = 7 bold times roman. This option is here
#'     because, for some reason, many people prefer to use fnt=2, which is some
#'     pale sans-serif difficult to read horror. But we must recognize that
#'     tastes differ, though really we should encourage people to be more
#'     considerate.
#' @param setpar if this is a stand-alone plot setting this TRUE will define a
#'     plot size and set the par values. If the plot is to be part of a
#'     multiple plot then setpar should be FALSE. default=TRUE
#'
#' @return a list of the median mature biomass depletion proxy and median
#'     harvest rate proxy by year, invisibly, it also plots these vectors.
#' @export
#'
#' @examples
#' \dontrun{
#' # eample syntax
#' tasphaseplot(proxyB=outhcr$g4s,proxyH=outhcr$targsc,glb=glb,sau=6,
#'              rundir="",console=TRUE,fnt=7)
#' }
tasphaseplot <- function(proxyB,proxyH,glb,sau,rundir="",console=TRUE,fnt=7,
                         setpar=TRUE){
  # proxyH=outhcr$g4s; proxyB=outhcr$targsc; glb=glb; sau=sau;console=console; fnt=7
  # setpar=FALSE;
  targdepl <- 5
  limdepl <- 1
  limH <- 5
  yrs <- as.numeric(rownames(as.matrix(proxyH[,,1])))
  nyrs <- length(yrs)
  filen <- ""
  if ((!console) & (setpar)) {
    nfile <- paste0("Tasmanian_PhasePlot_for_",glb$saunames[sau],".png")
    filen <- pathtopath(rundir,nfile)
    caption <- paste0("Tasmanian PhasePlot_of Harvest Rate vs Mature Biomass for ",
                      glb$saunames[sau])
  }
  medH <- as.matrix(apply(proxyH[,sau,],1,median))
  usemedH <- 10 - medH # because a high score = low H so re3versal is needed
  meddepl <- as.matrix(apply(proxyB[,sau,],1,median))
  if (setpar) {
    plotprep(width=7,height=6,filename=filen,cex=0.9,verbose=FALSE)
    parset(cex.lab=1.25,margin=c(0.45,0.45,0.1,0.1),font=fnt)
  }
  maxH <- maxdepl <- 10
  labelx <- "eHS mean CPUE Target score ~ Abundance proxy for "
  labely <- "eHS (10 - CPUE_Gradient_4_Score) ~ Fishing Mortality proxy "
  plot(meddepl,usemedH,type="p",pch=16,xlim=c(0,10),ylim=c(0,10),xaxs="i",
       yaxs="i",xlab=paste0(labelx,glb$saunames[sau]),ylab=labely,yaxt="n")
  axis(side=2,at=seq(0,10,1),labels=seq(5,-5,-1)) # to match Craig's version
  phaseplotpolygons(maxy=maxH,maxx=maxdepl,targx=targdepl,limx=limdepl,limy=limH)
  suppressWarnings(arrows(x0=meddepl[1:(nyrs-1)],y0=(usemedH[1:(nyrs-1)]),
                          x1=meddepl[2:nyrs],y1=(usemedH[2:nyrs]),lwd=2,
                          length=0.075,code=2))
  #  when H values repeatedly zero one gets zero length arrows and a warning
  # this is no9t a problem so it is suppressed.
  points(meddepl[c(1,nyrs)],usemedH[c(1,nyrs)],pch=16,cex=3,col=c(6,4))
  legend("topright",legend=c(yrs[1],yrs[nyrs]),
         col=c(6,4),lwd=5,cex=1.5,bty="n")
  if ((!console) & (setpar)) {
    addplot(filen,rundir=rundir,category="phaseplot",caption)
  }
  return(invisible(c(mediandepl=meddepl[nyrs],medianH=medH[nyrs])))
} # end of tasphaseplot

#' @title twophaseplots generates a formal kobe plot plus a Tasmanian proxy plot
#'
#' @description twophaseplot generates a combination of the formal kobe plot of
#'     the predicted biomass depletion level against the predicted harvest rate,
#'     along with a Tasmanian proxy plot of the target CPUE score (proxy for
#'     biomass depletion) against the gradient 4 CPUE score (proxy for harvest
#'     rate).
#'
#' @param dyn a list of the dynamics for every sau. This function only requires
#'     the mature biomass depletion and the harvest rates
#' @param glb the globals object
#' @param outhcr the scores object from the hcr containing the targsc and g4s
#' @param sau which sau is to be plotted
#' @param kobeRP reference points for the eHS pseudo kobe plot. A vector of the
#'     targdepl, default=0.4, limdepl, default=0.2, and the limH, default=0.15
#' @param rundir the directory for all files and results, default=""
#' @param startyr default = 1992. Otherwise could use condC$yearCE[1]
#' @param console should the plot go to the console or be saved? Default=TRUE
#' @param fnt what font to use, default = 7 bold times roman. This option is here
#'     because, for some reason, many people prefer to use fnt=2, which is some
#'     pale sans-serif difficult to read horror. But we must recognize that
#'     tastes differ, though really we should encourage people to be more
#'     considerate.
#'
#' @return a vector of median depletion and H, and median Grad4 and TargCE, it
#'     also plots a graph containing two kobe plots if in Tas.
#' @export
#'
#' @examples
#' print("wait on data")
twophaseplots <- function(dyn,glb,outhcr,sau,kobeRP,rundir="",startyr=1992,
                          console=TRUE,fnt=7) {
# dyn=sauout;glb=glb;outhcr=outhcr;sau=1;kobeRP=kobeRP;rundir=rundir;startyr=1992
# console=TRUE;fnt=7
  filen <- ""
  if (!console) {
    nfile <- paste0("PhasePlots_for_",glb$saunames[sau],".png")
    filen <- pathtopath(rundir,nfile)
    caption <- paste0("PhasePlot_of Harvest Rate vs Mature Biomass for ",
                      glb$saunames[sau]," with a Tasmanian proxy phaseplot")
  }
  if (is.null(outhcr$targsc)) {
    plotprep(width=7,height=6,filename=filen,cex=0.9,verbose=FALSE)
    parset(plots=c(1,1),cex.lab=1.25,margin=c(0.45,0.45,0.1,0.1))
  } else {
    plotprep(width=12,height=6,filename=filen,cex=0.9,verbose=FALSE)
    parset(plots=c(1,2),cex.lab=1.25,margin=c(0.45,0.45,0.1,0.1))
  }
  outkobe <- HSphaseplot(dyn=dyn,glb=glb,sau=sau,rundir=rundir,
                         startyr=startyr,console=console,setpar=FALSE,
                         targdepl=kobeRP[1],limdepl=kobeRP[2],limH=kobeRP[3])
  outtas <- c(0,0)
  if (!is.null(outhcr$targsc)) {
     outtas <- tasphaseplot(proxyB=outhcr$targsc,proxyH=outhcr$g4s,glb=glb,
                            sau=sau,rundir=rundir,console=console,
                            setpar=FALSE,fnt=fnt)
  }
  if (!console) {
    addplot(filen,rundir=rundir,category="phaseplot",caption)
  }
  return(c(outkobe,outtas))
} # end of twophaseplots

