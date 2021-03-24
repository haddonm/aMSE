
#' @title dosauplot generates the plot of the chosen variable for each sau
#'
#' @description dosauplot takes the results of the prepareproject and the
#'     projection function, it combines the trajectories across the years of
#'     conditioning and the projection years, for all replicates, and plots all
#'     replicates for the years from 'startyr' to the final year of the
#'     projections. It includes the median and inner 90% CI, and, for the CPUE
#'     plot, if the historical conditioning CPUE is included in histCE, it also
#'     includes the original cpue series.
#'
#' @param ylabel the variable name to be plotted, used as a y-axis label
#' @param prerep the replicate values from the replicated initial zoneDDR after
#'     recruitment variability has been added to a given number of years for a
#'     selected variable taken from the output of 'zonetosau'.
#' @param postrep the replicate values from the replicated projection zoneDP for
#'     the variable selected
#' @param glb the global constants object
#' @param startyr the index for the first year of the conditioning dynamics to
#'     include
#' @param addCI should quantile confidence bounds be included on the plots,
#'     default=FALSE
#' @param CIprobs the quantiles used to generate the confidence intervals. The
#'     default = c(0.05,0.5,0.95)
#' @param histCE the historical cpue data if included, default = NULL
#'
#' @return invisibly, a list of CI and median for each SAU
#' @export
#'
#' @examples
#' print("wait on suitable built in data sets")
dosauplot <- function(ylabel,prerep,postrep,glb,startyr,addCI=FALSE,
                      CIprobs=c(0.05,0.5,0.95),histCE=NULL) {
  label <- glb$saunames
  nsau <- glb$nSAU
  sauCI <- vector("list",nsau)
  names(sauCI) <- glb$saunames
  preyrs <- dim(prerep)[1]
  projyrs <- dim(postrep)[1]
  allyrs <- preyrs + projyrs
  reps <- dim(prerep)[3]
  nplot <- getparplots(nsau)
  if (is.numeric(histCE)) ceyr <- startyr:preyrs
  parset(plots=nplot,byrow=FALSE)
  for (sau in 1:nsau) {
    ymax <- getmax(rbind(prerep[startyr:preyrs,sau,],postrep[,sau,]))
    traj <- c(prerep[startyr:preyrs,sau,1],postrep[,sau,1])
    plot(startyr:allyrs,traj,type="l",lwd=1,col="grey",panel.first=grid(),
         ylab=paste0(ylabel,"  ",label[sau]),ylim=c(0,ymax),xlab="Years")
    for (i in 2:reps)
      lines(startyr:allyrs,c(prerep[startyr:preyrs,sau,i],postrep[,sau,i]),
            lwd=1,col="grey")
    CI <- apply(postrep[,sau,],1,quantile,probs=CIprobs)
    sauCI[[sau]] <- CI
    lines((preyrs+1):allyrs,CI[2,],lwd=2,col=4)
    if (addCI) {
      lines((preyrs+1):allyrs,CI[1,],lwd=1,col=2)
      lines((preyrs+1):allyrs,CI[3,],lwd=1,col=2)
    }
    if (is.numeric(histCE)) {
      oldce <- tail(histCE[,sau],length(ceyr))
      lines(ceyr,oldce,lwd=2,col=3)
    }
    abline(v=(preyrs+0.5),col=2)
    if (ylabel == "cpue") abline(h=150,col=2)
  }
  return(invisible(sauCI))
} # end of dosauplot


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
  if (length(xlimit) ==1 ) {
    xlimit <- c(0,getmax(x))
  }
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


#' @title plotproj aids the plotting of a projected variable
#'
#' @description plotproj plots out the projections from teh mse for a selected
#'     variable. iters is included so that details of a few trajectories can
#'     be viewed as well as seeing all replicates.
#'
#' @param invar the array containing the year x SAU x reps values for a selected
#'     variable out of sauzoneDP
#' @param varlabel the label to place along the outer Y-axis
#' @param plotconst a list containing nsau, saunames, reps, projyrs, and plts
#'     a vector of two describing the layout of plots, default=c(4,2)
#' @param vline default=NULL, which means nothing extra is plotted. If given a
#'     year then a vertical line in red will be added to the plot
#' @param iters default=0, which means all iterations will be plotted. If iters
#'     has a value then only that many trajectories will be plotted
#' @param addqnts will calculate and add the median and 90 percent quantiles
#' @param miny sets the lower limit of the y-axis, default=0
#'
#' @return nothing but it does generate a plot
#' @export
#'
#' @examples
#' print("wait on data files")
plotproj <- function(invar,varlabel,plotconst,miny=0,
                     vline=NULL,iters=0,addqnts=FALSE) {
  nsau <- plotconst$nsau
  yrs <- 1:plotconst$projyrs
  saunames <- plotconst$saunames
  reps <- plotconst$reps
  parset(plots=plotconst$plts,byrow=FALSE,margin=c(0.25,0.4,0.1,0.05),
         outmargin=c(1,1,0,0))
  for (sau in 1:nsau) { # sau=1
    maxy <- getmax(invar[,sau,])
    plot(yrs,invar[,sau,1],lwd=1,type="l",col="grey",panel.first=grid(),
         ylim=c(miny,maxy),xlab="",ylab=saunames[sau])
    trajs <- reps
    if (iters > 0) trajs <- iters
    for (iter in 1:trajs) lines(yrs,invar[,sau,iter],lwd=1,col="grey")
    if (!is.null(vline)) abline(v=vline,lwd=1,col=2)
    if (addqnts) {
        CI <- apply(invar[,sau,],1,quantile,probs=c(0.05,0.5,0.95))
        lines(yrs,CI[1,],lwd=1,col=4)
        lines(yrs,CI[2,],lwd=2,col=4)
        lines(yrs,CI[3,],lwd=1,col=4)
    }
  }
  mtext("Years",side=1,line=-0.1,outer=TRUE,cex=1.0,font=7)
  mtext(varlabel,side=2,line=-0.1,outer=TRUE,cex=1.0,font=7)
} # end of plotproj


#' @title plotzoneproj plots the replicates of a single zone-wide variable
#'
#' @description plotzoneproj takes the results from the function aszone, which
#'     contains a set of variables zonesB, zoneeB, zoneC, zoneH, zoneR, zonece,
#'     zonedeplsB,and zonedepleB, and plots whichever is selected. It contains
#'     the option of also plotting the median and the inner 90 percent quantiles
#'     of each year's spread of values across replicates
#'
#' @param zoneV the variable from within the output of aszone
#' @param reps the number of replicates to plot. Generaly one would plot all of
#'     them. If not then it might be best to turn addqnts to FALSE
#' @param yrs the vector of years as in 1:(inityrs + projyrs)
#' @param label the y-axis label, which should obviously reflect which variable
#'     is chosen from the zone summary list.
#' @param addqnts should the median and inner 90 percent quantiles be added to
#'     the plot; default = TRUE
#' @param miny sets the lower limit of the y-axis, default=0
#'
#' @return if addqnts=TRUE the quantiles are returned invisibly
#' @export
#'
#' @examples
#' print("wait on more time")
#' # zoneV=zoneproj$zoneC;reps=reps;yrs=1:50;miny=0;label="Catches"; addqnts=TRUE
plotzoneproj <- function(zoneV,reps,yrs,label="",addqnts=TRUE,miny=0) {
  maxy <- getmax(zoneV)
  ylabel <- "Variable"
  if (nchar(label) > 0) ylabel <- label
  plot(yrs,zoneV[,1],type="n",panel.first=grid(),ylim=c(miny,maxy),yaxs="i",
       ylab=ylabel,xlab="Years")
  for (iter in 1:reps) lines(yrs,zoneV[,iter],lwd=1,col="grey")
  if (addqnts) {
    CI <- apply(zoneV,1,quantile,probs=c(0.05,0.5,0.95))
    lines(yrs,CI[1,],lwd=2,col=4)
    lines(yrs,CI[2,],lwd=2,col=2)
    lines(yrs,CI[3,],lwd=2,col=4)
  }
  if (addqnts) return(invisible(t(CI)))
} # end of plotzoneproj


#' @title sauplots generates the dosauplots for the dynamic variables
#'
#' @description sauplots generates the dosauplots for the dynamic variables
#'    including cpue, catch, acatch (aspirational catch), matureB, exploitB,
#'    recruit, and harvestR
#'
#' @param zoneDP the dynamic object produced by the projections
#' @param zoneDDR the dynamic conditioning object after preparing for replicates
#'     and the addition of recruitment variation in the final years.
#' @param glb the object containing the global variables
#' @param rundir the results directory
#' @param B0 the B0 values by population use getvar(zoneC,"B0")
#' @param ExB0 the ExB0 values by population use getvar(zoneC,"ExB0")
#' @param startyr the index for the first year of the conditioning dynamics to
#'     include in the trajectory, to give startyr:indexoflastyear,eg startyr=40
#' @param addCI should confidence intervals be added to envelope plots,
#'     default=TRUE
#' @param histCE historical CPUE data used in CPUE plots, default=NULL
#'
#' @return a list of lists of CI for each SAU and variable
#' @export
#'
#' @examples
#' print("wait on suitable internal data-sets")
sauplots <- function(zoneDP,zoneDDR,glb,rundir,B0,ExB0,startyr,addCI=TRUE,
                     histCE=NULL) {
  # zoneDP=zoneDP;zoneDDR=zoneDDR;glb=glb;rundir=rundir;B0=B0;ExB0=ExB0;addCI=TRUE;histCE=histCE
  zoneDsau <- zonetosau(zoneDDR,glb,B0,ExB0)
  zonePsau <- zonetosau(zoneDP,glb,B0,ExB0)
  label <-  c("cpue","catch","acatch","matureB","exploitB","recruit","harvestR")
  out <- vector("list",length(label))
  names(out) <- label
  nplot <- getparplots(glb$nSAU)
  #CPUE
  filen <- filenametopath(rundir,"proj_cpue_SAU.png")  # filen=""
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  CI <- dosauplot("cpue",zoneDsau[["cpue"]],zonePsau[["cpue"]],glb,
                  startyr=startyr,addCI=TRUE,histCE=histCE)
  caption <- "The CPUE projections for each SAU."
  addplot(filen,rundir=rundir,category="ProjSAU",caption)
  out[["cpue"]] <- CI
  #Catches
  filen <- filenametopath(rundir,"proj_catch_SAU.png")
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  CI <- dosauplot("catch",zoneDsau[["catch"]],zonePsau[["catch"]],glb,
                  startyr=startyr,addCI=TRUE,histCE=NULL)
  caption <- "The catch projections for each SAU."
  addplot(filen,rundir=rundir,category="ProjSAU",caption)
  out[["catch"]] <- CI
  #Aspirational catches
  filen <- filenametopath(rundir,"proj_aspcatch_SAU.png")
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  CI <- dosauplot("acatch",zoneDsau[["acatch"]],zonePsau[["acatch"]],glb,
                  startyr=startyr,addCI=TRUE,histCE=NULL)
  caption <- "The Aspirational catch projections for each SAU."
  addplot(filen,rundir=rundir,category="ProjSAU",caption)
  out[["acatch"]] <- CI
  #MatureBiomass
  filen <- filenametopath(rundir,"proj_matureB_SAU.png")
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  CI <- dosauplot("matureB",zoneDsau[["matureB"]],zonePsau[["matureB"]],glb,
                  startyr=startyr,addCI=TRUE,histCE=NULL)
  caption <- "The mature biomass projections for each SAU."
  addplot(filen,rundir=rundir,category="ProjSAU",caption)
  out[["matureB"]] <- CI
  #exploitable biomass
  filen <- filenametopath(rundir,"proj_exploitB_SAU.png")
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  CI <- dosauplot("exploitB",zoneDsau[["exploitB"]],zonePsau[["exploitB"]],glb,
                  startyr=startyr,addCI=TRUE,histCE=NULL)
  caption <- "The exploitable biomass projections for each SAU."
  addplot(filen,rundir=rundir,category="ProjSAU",caption)
  out[["exploitB"]] <- CI
  #recruitment
  filen <- filenametopath(rundir,"proj_recruit_SAU.png")
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  CI <- dosauplot("recruit",zoneDsau[["recruit"]],zonePsau[["recruit"]],glb,
                  startyr=startyr,addCI=TRUE,histCE=NULL)
  caption <- "The recruitment projections for each SAU."
  addplot(filen,rundir=rundir,category="ProjSAU",caption)
  out[["recruit"]] <- CI
  filen <- filenametopath(rundir,"proj_harvestR_SAU.png")
  plotprep(width=8,height=8,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  CI <- dosauplot("harvestR",zoneDsau[["harvestR"]],zonePsau[["harvestR"]],glb,
                  startyr=startyr,addCI=TRUE,histCE=NULL)
  caption <- "The Harvest Rate projections for each SAU."
  addplot(filen,rundir=rundir,category="ProjSAU",caption)
  out[["harvestR"]] <- CI
  return(invisible(out))
}# end of sauplots



