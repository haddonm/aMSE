

#' @title product is the productivity curve matrix from doproduction
#'
#' @description product is the productivity curve matrix from
#'     doproduction when the example zone is generated using the
#'     inbuilt datasets ctrl, zone1, and constants. The slowest
#'     part of building the whole is to use the modregC function
#'     to adjust the zoneC and generate the production array. To
#'     save that time in the examples (to avoid time limits on
#'     examples should this package go to CRAN), then this dataset can
#'     be used instead. This is a three dimensional array of
#'     productivity variables.
#'
#' @name product
#'
#' @docType data
#'
#' @section contents:
#' \itemize{
#'   \item harvestrate the initial harvest rates applied
#'   \item productivity variables ExB, MatB, AnnH, Catch, Deplet, RelCE
#'   \item population the index of each population
#' }
#'
#' @examples
#'  data(product)
#'  product[1:20,,1]
NULL

#' @title testzoneC is a zone list made up of 6 equilibrium populations
#'
#' @description testzoneC is a zone list made up of 6 equilibrium
#'     populations. These have been run with a laral dispersal rate of
#'     0.03 so the change from B0 to effB0 is not great, but still
#'     required for an initial equilibrium. This is here to simplify
#'     the internal testing of funcitons that require a completed
#'     zone starting at equilibrium. Its name is to avoid conflict
#'     with any actual use of zoneC. use str(testzoneC, max.level=1)
#'     to see its format. It can be expected to be used with testzoneD
#'
#' @name testzoneC
#'
#' @docType data
#'
#' @section Subjects:
#'  \itemize{
#'    \item testing of functions that require a full zone
#'    \item initial equilibrium
#'  }
#'  @export
#'
#' @examples
#'  data(testzoneC)
#'  data(testzoneD)
#'  data(zone1)
#'  glb <- zone1$globals
#'  r0 <- getvar(testzoneC,"R0")
#'  move <- makemove(glb$numpop,r0,glb$larvdisp)
#'  glb$move <- move
#'  ans <- testequil(testzoneC, testzoneD, glb)
#'  str(testzoneC[[1]])
NULL

#' @title testzoneD is a list of 8 matrices and 2 arrays defining the dynamics of a zone
#'
#' @description testzoneD is a list of 8 matrices and 2 arrays defining
#'     the dynamics of a zone. These have been run with a larval
#'     dispersal rate of 0.03 to achieve an initial equilibrium. This
#'     is here to simplify the internal testing of functions that
#'     require a completed zone starting at equilibrium. Its name is
#'     to avoid conflict with any actual use of zoneD. use
#'     str(testzoneD, max.level=1) to see its format. It can be
#'     expected to be used with testzoneC.
#'
#' @name testzoneD
#'
#' @docType data
#'
#' @section Subjects:
#'  \itemize{
#'    \item testing of functions that require a full zone
#'    \item initial equilibrium
#'  }
#'  @export
#'
#' @examples
#'  data(testzoneC)
#'  data(testzoneD)
#'  data(zone1)
#'  glb <- zone1$globals
#'  r0 <- getvar(testzoneC,"R0")
#'  move <- makemove(glb$numpop,r0,glb$larvdisp)
#'  glb$move <- move
#'  ans <- testequil(testzoneC, testzoneD, glb)
#'  str(testzoneD)
NULL





#' @title plotzoneproj plots the replicates of a single zone-wide variable
#'
#' @description plotzoneproj takes the results from the function aszone, which
#'     contains a set of variables zonesB, zoneeB, zoneC, zoneH, zoneR, zonece,
#'     zonedeplsB,and zonedepleB, and plots whichever is selected. It contains
#'     the option of also plotting the median and the inner 90 percent quantiles
#'     of each year's spread of values across replicates. DEPRECATED
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

#' @title plotproj aids the plotting of a projected variable
#'
#' @description plotproj plots out the projections from the MSE for a selected
#'     variable. iters is included so that details of a few trajectories can
#'     be viewed as well as seeing all replicates. DEPRECATED
#'
#' @param invar the array containing the year x SAU x reps values for a selected
#'     variable out of sauzoneDP
#' @param varlabel the label to place along the outer Y-axis
#' @param plotconst a list containing nsau, saunames, reps, projyrs, and plts
#'     a vector of two describing the layout of plots
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
