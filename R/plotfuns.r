

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
#' @param xlab the x-axis label, default=""
#' @param ylab the y-axis label, default=""
#'
#' @return invisibly a list of the x and y matrices plotted
#' @export
#'
#' @examples
#' print("wait for a data set")
plotprod <- function(product,xname="MatB",yname="Catch",xlab="",
                     ylab="Production") {
  x <- product[,xname,]
  y <- product[,yname,]
  numpop <- ncol(x)
  maxy <- getmax(y)
  maxx <- getmax(x)
  plot(x[,1],y[,1],type="l",lwd=2,col=1,ylim=c(0,maxy),xlim=c(0,maxx),
       xlab=xlab,ylab=ylab,panel.first=grid(),yaxs="i")
  for (i in 2:numpop) lines(x[,i],y[,i],lwd=2,col=i)
  return(invisible(list(xname=x,yname=y)))
} # end of plotprod

