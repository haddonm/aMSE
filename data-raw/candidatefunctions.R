
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

# cpue <- zoneDRp$cpue[1:10,1,1]
# years <- 1:10
# getgrad1(cpue)
# getgrad4(cpue)
#
# plotprep()
# parset()
# plot1(years,cpue,lwd=2)



#condcpue <- MCDAdata(resdir,projC$HS,zone1$SAUnames)


vectce <- c(0.867862817,0.885256733,0.877541942,1.008936892,1.169236848,
            1.231699353,1.226890881,1.237515477,1.238744573,1.156193697,
            1.069769077,1.086021344,0.954323058,0.956561516,0.982457905,
            1.001020146,0.951562922,1.040046115,0.970610956,1.016168543,
            0.978065368,0.862886265,0.876245926,0.832028769,0.901777782,
            0.914997995,0.925753845,0.779823254)
names(vectce) <- 1992:2019

quantile(ce,probs=c(0.55))

one <- getgrad4(vectce)
two <- getgrad4(scaleto1(vectce))
cbind(one,two,two-one)

one <- targscore(vectce)
two <- targscore(scaleto1(vectce))
cbind(one$ans[,"targsc"],two$ans[,"targsc"],two$ans[,"targsc"]-one$ans[,"targsc"])


# hcr <- c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2)
# names(hcr) <- c(1:10)
# wid=4
# targqnt=0.55
# pmwts=c(0.65,0.25,0.1)
#
#
# mcda <- mcdahcr(vectce,only=FALSE)
#
# plot(vectce,type="l",ylim=c(0,1.4))
# abline(h=c(limrp,targ,uprp),col=c(2,3,2))
#
  # requirements for doprojection
  zoneCP=zoneCP;zoneDP=zoneDR;glob=glb;ctrl=ctrl;projyrs=projC$projyrs;inityrs=projC$inityrs;

  yr=28
  hcr <- c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2)
  names(hcr) <- c(1:10)
  wid=4
  targqnt=0.55
  pmwts=c(0.65,0.25,0.1)
  only=TRUE


library(microbenchmark)

microbenchmark(
mcdahcr(vectce,28),
oldmcdahcr(vectce,28)
)


zoneDP=zoneDRp
year=11; iter=1

ce <- zoneDP$cpue[1:(year - 1),,iter]
ce
out <- apply(ce,2,mcdahcr)

ces <- apply(ce,2,scaleto1)
outs <- apply(ces,2,mcdahcr)


library(aMSE)
library(rutilsMH)
library(makehtml)

resdir <- "C:/Users/User/Dropbox/A_Code/aMSEUse/conddata/generic"
ctrlzonetemplate(resdir,"control2.csv")


zone <- readctrlfile(resdir,infile="control2.csv")






x <- c(1,2,3,4,5)
wts <- c(1,1,1,1,1)

wtedmean(x,wts)


library(microbenchmark)

microbenchmark(
wtedmean(saucpue,saucatch),
wtmean(saucpue,saucatch)
)



reps=10
hr <- zoneDD$harvestR[1,]
ans <- matrix(0,nrow=reps,ncol=16)
for (pop in 1:16) ans[,pop] <- rnorm(reps,mean=hr[pop],sd=0.003)

ans



# newreadctrlfile ---------------------------------------------------------

readctrlfile2 <- function(datadir,infile="control.csv") {
  # datadir=resdir; infile="control2.csv"
  filenames <- dir(datadir)
  if (length(grep(infile,filenames)) != 1)
    stop(cat(infile," not found in datadir \n"))
  filename <- filenametopath(datadir,infile)
  indat <- readLines(filename)   # reads the whole file as character strings
  begin <- grep("START",indat) + 1
  runlabel <- getStr(indat[begin],1)
  datafile <- getStr(indat[begin+1],1)
  batch <- getsingleNum("batch",indat)
  reps <- getsingleNum("replicates",indat)
  withsigR <- getsingleNum("withsigR",indat)
  withsigB <- getsingleNum("withsigB",indat)
  withsigCE <- getsingleNum("withsigCE",indat)
  Nyrs=40 # to set up equilibrium unfished population; could be read in
  if (length(grep(datafile,filenames)) != 1)
    stop("population data file not found \n")
  cat("All required files appear to be present \n")
  # Now read zone data
  context <- "control_file"
  nSAU <-  getsingleNum("nSAU",indat) # number of spatial management units
  begin <- grep("SAUpop",indat)
  SAUpop <-  getConst(indat[begin],nSAU) # how many populations per SAU
  numpop <- sum(SAUpop)
  SAUnames <- getStr(indat[begin+1],nSAU)
  initdepl <- getConst(indat[begin+2],nSAU)
  minc <-  getsingleNum("minc",indat) # minimum size class
  cw    <- getsingleNum("cw",indat) # class width
  Nclass <- getsingleNum("Nclass",indat) # number of classes
  midpts <- seq(minc,minc+((Nclass-1)*cw),2)
  larvdisp <- getsingleNum("larvdisp",indat)
  randomseed <- getsingleNum("randomseed",indat)
  randomseedP <- getsingleNum("randomseedP",indat)
  initLML <- getsingleNum("initLML",indat)
  projyrs <- getsingleNum("PROJECT",indat)
  firstyear <- getsingleNum("firstyear",indat)
  inityrs <- getsingleNum("inityrs",indat)
  outyear <- c(projyrs,firstyear)
  projLML <- NULL
  HS <- NULL
  histCatch <- NULL
  histyr <- NULL
  histCE <- NULL
  yearCE <- NULL
  compdat=NULL
  if (projyrs > 0) {
    projLML <- numeric(projyrs)
    from <- grep("PROJLML",indat)
    for (i in 1:projyrs) {
      from <- from + 1
      projLML[i] <- getConst(indat[from],1)
    }
    begin <- grep("HARVESTS",indat)
    HS <- getStr(indat[begin],1)
    if (HS == "MCDA") # if HSdetail=1 then get the CPUE calibration data
      HSdetail <- getsingleNum("HSdetail",indat)
    if ((HS == "MCDA") & (HSdetail == 1)) {
      yrce <- getsingleNum("CEYRS",indat)
      if (yrce == 0) {
        warning("CPUE calibration has no data")
      } else {
        begin <- grep("CPUE",indat)
        histCE <- matrix(NA,nrow=yrce,ncol=nSAU)
        yearCE <- numeric(yrce) # of same length as nSAU
        colnames(histCE) <- SAUnames
        for (i in 1:yrce) {
          begin <- begin + 1
          cenum <- as.numeric(unlist(strsplit(indat[begin],",")))
          yearCE[i] <- cenum[1]
          histCE[i,] <- cenum[2:(nSAU+1)]
        }
        rownames(histCE) <- yearCE
      }
    } # end of if HSdetail test
    if (HS == "constantCatch")  # get constant inTAC
      HSdetail <- as.numeric(removeEmpty(unlist(strsplit(indat[begin+1],","))))
  } # end of projyrs if test
  catches <- getsingleNum("CATCHES",indat)
  if (catches > 0) {
     Nyrs <- catches  # don't forget to add an extra year for initiation
     begin <- grep("CondYears",indat)
     histCatch <- matrix(0,nrow=catches,ncol=nSAU)
     colnames(histCatch) <- SAUnames
     histyr <- matrix(0,nrow=Nyrs,ncol=2)
     colnames(histyr) <- c("year","histLML")
     for (i in 1:catches) {
       begin <- begin + 1
       asnum <- as.numeric(unlist(strsplit(indat[begin],",")))
       histyr[i,] <- asnum[1:2]
       histCatch[i,] <- asnum[3:(nSAU+2)]
     }
  } # end of catches loop
  rownames(histCatch) <- histyr[,1]
  rownames(histyr) <- histyr[,1]
  yrce <- getsingleNum("CEYRS",indat)
  if (yrce > 0) {
    begin <- grep("CPUE",indat)
    histCE <- matrix(NA,nrow=yrce,ncol=nSAU)
    yearCE <- numeric(yrce) # of same length as nSAU
    colnames(histCE) <- SAUnames
    for (i in 1:yrce) {
      begin <- begin + 1
      cenum <- as.numeric(unlist(strsplit(indat[begin],",")))
      yearCE[i] <- cenum[1]
      histCE[i,] <- cenum[2:(nSAU+1)]
    }
  } # end of ceyrs loop
  rownames(histCE) <- yearCE
  sizecomp <- getsingleNum("SIZECOMP",indat)
  compdat <- NULL
  if (sizecomp > 0) {
    lffiles <- NULL
    locsizecomp <- grep("SIZECOMP",indat)
    lffiles <- c(lffiles,getStr(indat[locsizecomp+i],1))
    compdat <- vector("list",sizecomp)
    for (i in 1:sizecomp) {
      filename <- filenametopath(datadir,lffiles[i])
      compdat[[i]] <- read.csv(file=filename,header=TRUE)
    }
  } # end of sizecomp loop
  condC <- list(histCatch=histCatch,histyr=histyr,
                histCE=histCE,yearCE=yearCE,initdepl=initdepl,
                compdat=compdat,Sel=NULL,SelWt=NULL)
  projC <- list(projLML=projLML,HS=HS,HSdetail=HSdetail,projyrs=projyrs,
                inityrs=inityrs,Sel=NULL,SelWt=NULL,histCE=histCE)
  outctrl <- list(runlabel,datafile,batch,reps,randomseed,randomseedP,
                  withsigR,withsigB,withsigCE,catches,projyrs)
  names(outctrl) <- c("runlabel","datafile","batch","reps","randseed",
                      "randseedP","withsigR","withsigB","withsigCE",
                      "catchyrs","projection")
  globals <- list(numpop=numpop, nSAU=nSAU, midpts=midpts,
                  Nclass=Nclass, Nyrs=Nyrs,larvdisp=larvdisp)
  totans <- list(SAUnames,SAUpop,minc,cw,larvdisp,randomseed,
                 initLML,condC,projC,globals,outctrl,catches,projyrs)
  names(totans) <- c("SAUnames","SAUpop","minc","cw","larvdisp","randomseed",
                     "initLML","condC","projC","globals","ctrl",
                     "catchyrs","projyrs")
  return(totans)
} # end of readctrlzone







zone1 <- readctrlfile2(resdir,infile="control2.csv")
ctrl <- zone1$ctrl
glb <- zone1$globals     # glb without the movement matrix
constants <- readdatafile(glb$numpop,resdir,ctrl$datafile)















