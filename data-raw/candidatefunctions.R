
#' @title dohistoricC imposes the historical catches on an unfished zone
#'
#' @description dohistoricC is used during the conditioning of the zone/region
#'     and imposes the historical catches, held in the zone data object
#'     obtained using readzonefile, onto an unfished initial zone, and it does
#'     this by imposing the time-series of catches to each SAU/block. The
#'     operation is through the use of oneyearC, which imposes one year's
#'     catches, which are a vector of SAU catches for each year of the series.
#'
#' @param zoneDD The input unfished dynamic zone, zoneD, object
#' @param zoneC the zone constants object, zoneC
#' @param zone1 the zone data object obtained from readzonefile
#'
#' @return a zoneD object
#' @export
#'
#' @examples
#' print("wait on some data sets")
dohistoricC <- function(zoneDD,zoneC,zone1) {
  glb <- zone1$globals
  histC <- zone1$histCatch
  yrs <- zone1$histyr[,"year"]
  nyrs <- length(yrs)
  for (yr in 2:nyrs) {  # yr=2 # ignores the initial unfished year
    year <- yrs[yr]
    catchpop <- histC[yr,]
    zoneDD <- oneyearC(zoneC=zoneC,zoneD=zoneDD,Ncl=glb$Nclass,
                       catchp=catchpop,year=yr,sigmar=1e-08,
                       npop=glb$nSAU,movem=glb$move)
  }
  return(zoneDD)
} # end of dohistoricC

#' @title plotCPUE plots the scaled historic cpue against the predicted CPUe
#'
#' @description plotCPUE rescales the predicted CPUE from the OM to match the
#'     duration of each time-series of observed historical CPUE, currently it
#'     assumes the same number of years of CPUE for each observed series, which
#'     is wrong for blocks 6 and 13W!#'
#'
#' @param predCE the predicted cpue for each SAU from zoneD object
#' @param histCE the historical cpue for each SAU from the fishery
#' @param pickyr the indices of the years of overlap between the historical and
#'     predicted time-series of cpue
#' @param sau the SAU to be plotted. default=1:8 to suit the west coast of TAS
#'
#' @return invisibly the rescaled predCE matrix.
#' @export
#'
#' @examples
#' print("wait on some data sets and this function being adopted")
plotCPUE <- function(predCE,histCE,pickyr,sau=1:8) {
  plotprep(width=6,height=7,newdev=FALSE)
  parset(plots=c(4,2),margin=c(0.25,0.25,0.05,0.01),outmargin = c(1,1,0,0.4))
  for (i in sau) { # i = 1
    pCE <- predCE[,i]/mean(predCE[pickyr,i])
    ymax <- getmax(pCE)
    plot(yrs,pCE,type="l",lwd=2,ylim=c(0,ymax),xlab="",ylab="",col=4,
         panel.first=grid())
    lines(yrs[pickyr],scaletoOne(histCE[,i]),lwd=2,col="darkorange")
    mtext(text=zone1$SAUnames[i],side=3,line=-1.1,cex=1.0)
  }
  mtext("Predicted vs Observed CPUE",side=2,cex=1.0,outer=TRUE,line=-0.35)
  return(invisible(pCE))
} # end of plotCPUE

#' @title plotSAUdepl plots the depletion time-series for all SAU
#'
#' @description plotSAUdepl is a convenient plotting routine designed for the
#'     west coast of Tasmania. It is currently not a general SAU depletion
#'     plotting routine, but can potentially be generalized.
#'
#' @param outB the matrix of depletion values from the zoneD object
#'
#' @return nothing but it does plot a set of graphs
#' @export
#'
#' @examples
#'  print("wait on some data sets")
plotSAUdepl <- function(outB) {
  plotprep(width=7,height=6,newdev = FALSE)
  parset(plots=c(3,1))
  ymax <- getmax(outB[,c(1,2,3,8)])
  plot(yrs,outB[,1],type="l",lwd=2,ylim=c(0,ymax),panel.first=grid(),ylab="6,7,8,13")
  for (i in c(2,3,8)) lines(yrs,outB[,i],lwd=2,col=i)
  legend("topright",legend=c(6,7,8,13),col=c(1,2,3,8),lwd=3,bty="n",cex=1.0)
  ymax <- getmax(outB[,c(4,5,7)])
  plot(yrs,outB[,4],type="l",lwd=2,ylim=c(0,ymax),col=4,panel.first=grid(),
       ylab=c("9,10,12"))
  for (i in c(5,7)) lines(yrs,outB[,i],lwd=2,col=i)
  legend("topright",legend=c(9,10,12),col=c(4,5,7),lwd=3,bty="n",cex=1.0)
  ymax <- getmax(outB[,6])
  plot(yrs,outB[,6],type="l",lwd=2,ylim=c(0,ymax),col=6,panel.first=grid(),
       ylab="11")
  legend("topright",legend=c(11),col=c(6),lwd=3,bty="n",cex=1.0)
} # end of plotSAUdepl

