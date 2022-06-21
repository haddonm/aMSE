

#' @title comparevar generates the quantiles for eahc of a set of input scenarios
#'
#' @description comparevar generates quantiles for each of a set of input
#'     scenarios. It takes the full timeline of dynamics and outputs just the
#'     projections and the quantiles of those projections
#'
#' @param dyn this is a list of the out$sauout$zonePsau produced by each
#'     scenario. This has to be be generated from teh saved RData files from
#'     each scenario
#' @param glbc a list of the global objects from each scenario being compared
#' @param scenes a  list of the out$ctrl$runlabel from each scenario
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
#'
#' @return nothing but it does generate a plot with nsau panels
#' @export
#'
#' @examples
#'  print("wait on internal datasets")
plotscene <- function(scenquant,glbc,var="cpue",ymin=0,filen="",
                      legloc="topleft",legplot=1) {
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
    lines(yrs,scenquant[[1]][[j]][2,],lwd=1,col=1)
    lines(yrs,scenquant[[1]][[j]][4,],lwd=1,col=1)
    if (legplot == j) {
      legend(legloc,legend=scenes,col=1:nscen,lwd=3,bty="n",cex=1.0)
    }
    for (i in 2:nscen) {
      lines(yrs,scenquant[[i]][[j]][3,],lwd=3,col=i)
      lines(yrs,scenquant[[i]][[j]][2,],lwd=1,col=i)
      lines(yrs,scenquant[[i]][[j]][4,],lwd=1,col=i)
    }
  }
  mtext(var,side=2,line=-0.2,outer=TRUE,cex=1.2)
  if (nchar(filen) > 0) dev.off()
} # end of plotscene

