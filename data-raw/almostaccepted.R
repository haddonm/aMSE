


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



#' @title getavrec pulls out just the AvRec values for each sau for a scenario
#'
#' @description getavrec extracts the AvRec values for each SAU for a given
#'     scenario. It does this by reading in the saudata file and using grep
#'     to search for the correct line and returning the line as is. Care is
#'     required to only take the first record of 'AvRec' so as not to
#'     consider the variability in sAvRec.
#'
#' @param datadir the directory in which one finds the saudata file for the
#'     given scenario
#' @param datafile the exact name of the saudata file
#' @param nsau the number of SAU to be found
#'
#' @return a numeric vector of the average recruitment for each SAU
#' @export
#'
#' @examples
#' \dontrun{
#'   rundir <- "c:/Users/User/DropBox/A_codeUse/aMSEUse/scenarios/MhLML/M1h5/"
#'   getavrec(rundir,"saudataM1h5.csv",8)
#' }
getavrec <- function(datadir,datafile,nsau) {
  filen <- filenametopath(datadir,datafile)
  dat <- readLines(filen)
  pickA <- grep("AvRec",dat)[1] # ignore sAvRec
  avrec <- getConst(dat[pickA],nsau)
  return(avrec)
} # end of getavrec




#' @title sauavrecssq applies the AvRec and returns the ssq from the cpue
#'
#' @description sauavrecssq applies the AvRec vector while adjusting just
#'     one for the picksau, and then compares the predicted CPUE with the
#'     observed and returns a sum-of-squared differences. This is used by
#'     optim while using the Brent option to find an optimum for a single
#'     parameter.
#'
#' @param param is the AvRec for the selected SAU
#' @param rundir the rundir for the scenario being considered
#' @param datadir the directory holding the saudata file, usually = rundir
#' @param controlfile the name of the control file, no path
#' @param datafile the name of the data file, no path
#' @param linenum the linenumber containing the AvRec vector
#' @param calcpopC the HS function that calculates the aspirational catches
#'     for each SAU
#' @param extra the vector of SAU AvRecs minus the selected SAU value
#' @param picksau which SAU is to be worked on as in 1:nsau
#' @param nsau the toal number of SAU
#'
#' @return the SSQ for the selected SAU in a comparison of predicted and
#'     observed CPUE
#' @export
#'
#' @examples
#' print("Wait on suitable internal data files")
sauavrecssq <- function(param,rundir,datadir,controlfile,datafile,linenum,
                        calcpopC,extra,picksau,nsau) {
  nsaum1 <- nsau - 1
  if (picksau == 1)
    replacetxt <- paste0("AvRec,",as.character(param),",",
                         paste0(as.character(extra[1:nsaum1]),collapse=","),
                         collapse=",")
  if ((picksau > 1) & (picksau < nsau))
    replacetxt <- paste0("AvRec,",
                         paste0(as.character(extra[1:(picksau-1)]),collapse=","),
                         ",",as.character(param),",",
                         paste0(as.character(extra[picksau:nsaum1]),collapse=","))
  if (picksau == nsau)
    replacetxt <- paste0("AvRec,",
                         paste0(as.character(extra[1:nsaum1]),collapse=","),
                         ",",as.character(param),collapse=",")
  changeline(datadir,datafile,linenum,replacetxt)
  zone <- makeequilzone(rundir,controlfile,datadir,doproduct=FALSE,
                        cleanslate=FALSE,verbose=FALSE)
  # declare main objects
  glb <- zone$glb
  condC <- zone$zone1$condC
  zoneC <- zone$zoneC
  zoneD <- zone$zoneD
  #Condition on Fishery
  zoneDD <- dohistoricC(zoneD,zoneC,glob=glb,condC,calcpopC=calcpopC,
                        sigR=1e-08,sigB=1e-08)
  hyrs <- glb$hyrs
  sauindex <- glb$sauindex
  popB0 <- getlistvar(zoneC,"B0")
  B0 <- tapply(popB0,sauindex,sum)
  popExB0 <- getlistvar(zoneC,"ExB0")
  ExB0 <- tapply(popExB0,sauindex,sum)
  sauZone <- getsauzone(zoneDD,glb,B0=B0,ExB0=ExB0)
  ssq <- compareCPUE(condC$histCE,sauZone$cpue,glb,rundir,filen="SAU_compareCPUE.png")
  return(ssq[picksau])
} # end of sauavrecssq



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






