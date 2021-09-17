


#' @title combinepopN sums the population numbers-at-size into SAU totals
#'
#' @description combinepopN works with the Nt or catchN from zoneDD, in other
#'     words, with the size-composition data during the conditioning (i.e.
#'     no replicates). It uses the sauindex from the gob object to select the
#'     columns from the in the input matrix and then uses rowSums to combine
#'     the numbers-at-size.
#'
#' @param inN either zoneDD$Nt of zoneDD$catchN
#' @param glb the globals object
#'
#' @return a matrix of Nclass x hyrs x nSAU of numbers-at-size data by SAU
#' @export
#'
#' @examples
#' print ("wait on suitable test data")
combinepopN <- function(inN,glb) {
  nsau <- glb$nSAU
  sauindex <- glb$sauindex
  years <- glb$hyrnames
  hyrs <- glb$hyrs
  mids <- glb$midpts
  saunames <- glb$saunames
  sauNt <- array(0,dim=c(glb$Nclass,hyrs,nsau),
                 dimnames=list(mids,years,saunames))
  for (sau in 1:nsau) { # convert pops to SAU
    for (yr in 1:hyrs) {
      pick <- which(sauindex == sau)
      sauNt[,yr,sau] <- rowSums(inN[,yr,pick])
    }
  }
  return(sauNt=sauNt)
} # end of combinepopN


#' @title saucompcatchN compares the predicted vs observed size composition data
#'
#' @description saucompcatchN provides a comparison of the proportional
#'     distribution of the predicted size-composition with the observed size-
#'     composition for a given SAu for a given year or set of years. This
#'     function could be used with the population numbers-at-size relative
#'     to the numbers-at-size from survey data. It ay get modified for that
#'     purpose
#'
#' @param obs the observed numbers-at-size in the catch for all years for a
#'     given SAU
#' @param pred the predicted numbers-at-size in the catch for all years for
#'     a giuven SAU
#' @param glb the globals object
#' @param years which year or years should be plotted?
#' @param ssc What size class should the plots begin. One should select the
#'     minimum size seen in the observed catch sampling.
#' @param sau which SAU is being examined. For labelling. Default="SAU"
#' @param filen a filename, in preparation for adding these plots to the
#'     internal web-page (still under development). Default=""
#' @param wide how wide should the plot be? Default=8
#' @param tall how tall should the plot be? Default=7
#' @param newdev should a new lot be generated each time? Default=TRUE
#'
#' @return nothing yet, the intention is to include maybe a multinomial
#'     likelihood.
#' @export
#'
#' @examples
#' print ("wait on suitable test data")
#' # obs=cs;pred=sauP;glb=glb;years=1990:1999;ssc=68;sau=label;filen="";wide=10;tall=6;newdev=FALSE
saucompcatchN <- function(obs,pred,glb,years,ssc,sau="SAU",filen="",
                          wide=8,tall=7,newdev=TRUE) {
  pobs <- prop.table(obs,2)
  ppred <- prop.table(pred,2)
  mids <- glb$midpts
  pltrge <- ssc:glb$Nclass
  obsrge <- pltrge - ssc + 1
  testobs <- colSums(ppred)
  plotprep(width=wide,height=tall,newdev=newdev,filename=filen)
  parset(plots=getparplots(length(years)),byrow=FALSE,outmargin=c(1,1,0,0),
         margin=c(0.2,0.375,0.05,0.05))
  for (yr in years) { #  yr=1990
    pickY <- which(colnames(ppred) == yr)
    pickYcs <- which(colnames(pobs) == yr)
    ymax <- getmax(c(ppred[pltrge,pickY],pobs[obsrge,pickYcs]))
    nobs <- sum(obs[obsrge,pickYcs],na.rm=TRUE)
    plot(mids[pltrge],pobs[obsrge,pickYcs],type="l",lwd=1,ylim=c(0,ymax),
         yaxs="i",panel.first = grid(),xlab="",ylab=yr)
    if (testobs[pickY] > 0) lines(mids[pltrge],ppred[pltrge,pickY],lwd=2,col=2)
    text(mids[ssc+3],ymax*0.9,nobs,cex=1.2)
    if (yr == years[1])
      legend("topright",c("Obs","Pred"),col=c(1,2),lwd=3,bty="n",cex=1.2)
  }
  mtext("Shell Length mm",side=1,line=-0.2,outer=TRUE,cex=1.1)
  label <- paste0("Proportional Size Composition in ",sau)
  mtext(label,side=2,line=-0.25,outer=TRUE,cex=1.0)
}
