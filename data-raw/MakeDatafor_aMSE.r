#rm(list=ls())   # Cleanup the R console if required
# Set up the run ----------

options("show.signif.stars"=FALSE,"stringsAsFactors"=FALSE,
        "max.print"=50000,"width"=240)
#     listFunctions("C:/A_Mal/A_Book/rcode/agestruct/age_utils.r")

library(aMSE)
library(hutils)
library(hplot)
library(makehtml)

if (dir.exists("c:/Users/User/DropBox")) {
  ddir <- "c:/Users/User/DropBox/A_code/"
} else {
  ddir <- "c:/Users/Malcolm/DropBox/A_code/"
}
rundir <- paste0(ddir,"aMSEUse/conddata/generic2")
dirExists(rundir,make=TRUE,verbose=TRUE)

# equilibrium zone -------------------------------------------------------------
# You now need to ensure that there is, at least, a control.csv, zone1.csv
# and region1.csv file in the data directory plus some other data .csv files
# depending on how conditioned you want the model to be. Templates for the
# correct format can be produced using ctrlfiletemplate(), datafiletemplate(),
# and zonefiletemplate.
zone1 <- readctrlfile(rundir,infile="control.csv")
ctrl <- zone1$ctrl
constants <- readdatafile(glb$numpop,rundir,ctrl$datafile)

datadir <- "C:/Users/User/Dropbox/A_Code/aMSE/data/"

save(ctrl,file=paste0(datadir,"ctrl.RData"))
save(zone1,file=paste0(datadir,"zone1.RData"))
save(constants,file=paste0(datadir,"constants.RData"))


# make data reday for twosau single population examples-----------------
options("show.signif.stars"=FALSE,"stringsAsFactors"=FALSE,
        "max.print"=50000,"width"=240)
#     listFunctions("C:/A_Mal/A_Book/rcode/agestruct/age_utils.r")

  library(aMSE)
  library(hutils)
  library(TasHS)
  library(makehtml)
  dropdir <- getDBdir()
  prefixdir <- paste0(dropdir,"A_codeUse/aMSEUse/scenarios/")
  postfixdir <- "S21"
  verbose <- TRUE
  rundir <- filenametopath(prefixdir,postfixdir)
  ctrlfile <- paste0("control",postfixdir,".csv")

  doproduct=TRUE

  zone1 <- readctrlfile(rundir,infile=ctrlfile,verbose=verbose)
  ctrl <- zone1$ctrl
  glb <- zone1$globals     # glb without the movement matrix
  bysau <- zone1$ctrl$bysau
  opar <- NULL
  parsin <- FALSE
  saudata <- readsaudatafile(rundir,ctrl$datafile,optpar=opar)
  constants <- saudata$constants
  saudat <- saudata$saudat
  zone1$condC$poprec <- saudata$poprec

  if (verbose) cat("Files read, now making zone \n")
  out <- setupzone(constants,zone1,doproduct,verbose=verbose) # make operating model
  zoneC <- out$zoneC
  zoneD <- out$zoneD
  glb <- out$glb             # glb now has the movement matrix
  product <- out$product     # important bits usually saved in rundir
  zone1$globals <- glb


  datadir <- paste0(dropdir,"A_Code/aMSE/data/")

  save(zone1,file=paste0(datadir,"zone1.RData"))
  save(saudat,file=paste0(datadir,"saudat.RData"))





# Complete zone objects---------------------------------------------

library(rutilsMH)
library(aMSE)
library(microbenchmark)

datadir <- "./../../rcode2/aMSE/data/"
# read data files ----------------------------------------------------
rundir <- "./../../rcode2/aMSEUse/out/run1"
dirExists(rundir,make=TRUE,verbose=TRUE)
# You now need to ensure that there is a control.csv, zone1sau2pop6.csv
# and zone1.csv file in the data directory
ctrl <- checkctrldat(rundir)
runname <- ctrl$runlabel
zone1 <- readzonefile(rundir,ctrl$zonefile)
glb <- zone1$globals
constants <- readdatafile(glb$numpop,rundir,ctrl$datafile)

out <- setupzone(constants, zone1)
zoneC <- out$zoneC
zoneD <- out$zoneD
product <- out$product
glb <- out$glb
testzoneC <- zoneC
testzoneD <- zoneD

save(testzoneC,file=paste0(datadir,"testzoneC.RData"))
save(testzoneD,file=paste0(datadir,"testzoneD.RData"))
save(product,file=paste0(datadir,"product.RData"))


# some cpue data -----------------------------------------------------

datadir <- "./../../A_code/aMSE/data-raw/"

blockE13 <- read.csv(paste0(datadir,"block13e.csv"),header=TRUE)

save(blockE13,file=paste0(datadir,"blockE13.RData"))


# save zone data

datadir <- "./../../A_code/aMSE/data-raw/"

save(zone,file=paste0(datadir,"zone.RData"))



# save sau8 LF data------------------------------

# first run do_condition on an 8 sau scenario using lf_WZ90-20.csv

lfs <- out$condC$compdat$lfs

datadir <- "./../../A_code/aMSE/data/"

save(lfs,file=paste0(datadir,"lfs.RData"))


# save zoneC for 8 SAu 56pop


zoneC <- out$zoneC

save(zoneC,file=paste0(datadir,"zoneC.RData"))

zone <- makeequilzone(rundir,controlfile,doproduct=doproduct,
                      verbose=verbose)

save(zone,file=paste0(datadir,"zone.RData"))


# check and transfer -------------------------------------------------------






tools::checkRdaFiles(paths=datadir)
tools::resaveRdaFiles(paths=datadir,compress="auto")
tools::checkRdaFiles(paths=datadir)


















