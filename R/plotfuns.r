

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

#' @title plotsau generates nSAU plots in one figure for a given variable
#'
#' @description plotsau summarizes a single variable from the dynamic object
#'     produced by the replicate projections by generating nSAU plots in a
#'     single figure. These are all lineplots in grey (change using col), with
#'     the median in red (change using medcol).
#'
#' @param invar the three dimensional array containing the variable to plot
#' @param glb the global object
#' @param plots the number of plots.eg for 8 plots one might put c(4,2)
#' @param ylab the prefix for the Y labels, which will have the SAU number
#'     added to it. Default="Catch"
#' @param xlab The single X label for the bottom of the plot, default="Year"
#' @param col the colour for the reps replicate projection lines, default="grey"
#' @param medcol the colour for the median line, default="red", set to 0 or NA
#'     to avoid adding to the plots.
#' @param addCI should quantile confidence bounds be included on the plots,
#'     default=FALSE
#' @param CIcol if CI are to be included let their colour be this, default=blue
#' @param CIprobs the CI quantiles. defaults are c(0.05,0.95), the inner 90 perc
#'
#' @return invisibly the median line for each SAU
#' @export
#'
#' @examples
#' print("wait on appropriate internal data sets")
plotsau <- function(invar,glb,plots,ylab="Catch",xlab="Year",col="grey",
                    medcol="red",addCI=FALSE,CIcol="blue",CIprobs=c(0.05,0.95)) {
  parset(plots=plots,byrow=FALSE,margin=c(0.25,0.5,0.1,0.1),outmargin=c(1,0,0,0))
  nsau <- glb$nSAU
  label <- glb$saunames
  projyrs <- dim(invar)[1]
  yrs <- 1:projyrs
  reps <- dim(invar)[3]
  saumedian <- matrix(0,nrow=projyrs,ncol=nsau,dimnames=list(yrs),label)
  for (sau in 1:nsau) {  # sau=1
    ymax <- getmax(invar[,sau,])
    plot(yrs,invar[,sau,1],type="l",lwd=1,col=col,panel.first = grid(),
         ylim=c(0,ymax),yaxs="i",ylab=paste0(ylab,"    ",label[sau]),xlab="")
    for (i in 1:ctrl$reps) lines(1:projyrs,invar[,sau,i],lwd=1,col=col)
    for (yr in 1:projyrs) saumedian[yr,sau] <- median(invar[yr,sau,])
    lines(yrs,saumedian[,sau],lwd=2,col=medcol)
    if (addCI) {
      CI <- apply(invar[,sau,],1,quantile,probs=CIprobs)
      lines(yrs,CI[1,],lwd=1,col=CIcol)
      lines(yrs,CI[2,],lwd=1,col=CIcol)
    }
  }
  mtext(text=xlab,side=1,line=-0.1,outer=TRUE,cex=1.1)
  return(invisible(saumedian))
} # end of plotsau


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



