

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
#' data(product)
#' plotprod(product)
#' stat <- findmsy(product)
#' abline(h=stat[,"Catch"],col=1:6,lwd=2)
#' abline(v=stat[,"MatB"],col=1:6,lwd=2)
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



