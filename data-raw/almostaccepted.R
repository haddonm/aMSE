



options("show.signif.stars"=FALSE,
        "stringsAsFactors"=FALSE,
        "max.print"=50000,
        "width"=240)
suppressPackageStartupMessages({
  library(aMSE)
  library(TasHS)
  library(hutils)
  library(hplot)
  library(makehtml)
  library(knitr)
})
dropdir <- getDBdir()
prefixdir <- paste0(dropdir,"A_codeUse/aMSEUse/scenarios/")

startime <- Sys.time()
postfixdir <- "HS21"
verbose <- TRUE
rundir <- filenametopath(prefixdir,postfixdir)
controlfile <- paste0("control",postfixdir,".csv")
outdir <- "C:/aMSE_scenarios/M15h7L75/"
confirmdir(rundir)
confirmdir(outdir)


hsargs <- list(mult=0.1,
               wid = 4,
               targqnt = 0.55,
               maxtarg = c(150,150,150,150,150,150,150,150),
               pmwts = c(0.65,0.25,0.1),
               hcr = c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2),
               startCE = 1992)



zone1 <- readctrlfile(rundir,infile=controlfile,verbose=verbose)
ctrl <- zone1$ctrl
glb <- zone1$globals     # glb without the movement matrix
bysau <- zone1$ctrl$bysau
opar <- NULL
parsin <- zone1$condC$parsin
if (parsin) opar <- as.matrix(zone1$condC$optpars)
if (is.null(bysau)) bysau <- 0
if (bysau) {
  saudata <- readsaudatafile(rundir,ctrl$datafile,optpar=opar)
  constants <- saudata$constants
  saudat <- saudata$saudat
  zone1$condC$poprec <- saudata$poprec
} else {
  constants <- readpopdatafile(rundir,ctrl$datafile)
  saudat <- constants
}

if (verbose) cat("Files read, now making zone \n")
out <- setupzone(constants,zone1,doproduct=TRUE,verbose=verbose) # make operating model
zoneC <- out$zoneC
zoneD <- out$zoneD
glb <- out$glb             # glb now has the movement matrix
product <- out$product     # important bits usually saved in rundir
zone1$globals <- glb


zone <- makeequilzone(rundir,controlfile,doproduct=doproduct,verbose=verbose)


out <- do_condition(rundir,controlfile,
                    calcpopC=calcexpectpopC,
                    verbose = TRUE,
                    doproduct = TRUE,
                    dohistoric=TRUE,
                    mincount=120)

makeoutput(out,rundir,postfixdir,controlfile,hsfile="TasHS Package",
                 doproject=FALSE,openfile=TRUE,verbose=FALSE)





plotcpue <- function(zoneC,zoneDD,condC,glb,sau=1){
  dyn <- getdynamics(zoneC,zoneDD,condC,glb,sau=sau)
  parset()
  plot(yrs,dyn[,],type="l",lwd=2)


}


plotprep(width=8, height=5,newdev=FALSE)
plotcpue(glb$hyrnames,)


