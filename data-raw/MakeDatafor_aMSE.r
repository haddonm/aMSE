#rm(list=ls())   # Cleanup the R console if required
# Set up the run ----------


# library(r4cpue)
library(aMSE)


options("show.signif.stars"=FALSE,"stringsAsFactors"=FALSE,
        "max.print"=50000,"width"=240)
#     listFunctions("C:/A_Mal/A_Book/rcode/agestruct/age_utils.r")

# two block abalone zone with 6 populations -----------


# Setup data for aMSE -----------------------------------------------------

datadir <- "./../../rcode2/aMSE/data-raw/"

ctrlfile <- "control.csv"
#source(filenametopath(datadir,"sourcer.R"))

ctrl <- readctrlfile(datadir,ctrlfile)
zone1 <- readzonefile(datadir,ctrl$zonefile)
glb <- zone1$globals
constants <- readdatafile(glb$numpop,datadir,ctrl$datafile)



save(ctrl,file=paste0(datadir,"ctrl.RData"))
save(zone1,file=paste0(datadir,"zone1.RData"))
save(constants,file=paste0(datadir,"constants.RData"))


# Complete zone objects---------------------------------------------

library(rutilsMH)
library(aMSE)
library(microbenchmark)

datadir <- "./../../rcode2/aMSE/data-raw/"
# read data files ----------------------------------------------------
resdir <- "./../../rcode2/aMSEUse/out/run1"
dirExists(resdir,make=TRUE,verbose=TRUE)
# You now need to ensure that there is a control.csv, zone1sau2pop6.csv
# and zone1.csv file in the data directory
ctrl <- checkctrldat(resdir)
runname <- ctrl$runlabel
zone1 <- readzonefile(resdir,ctrl$zonefile)
glb <- zone1$globals
constants <- readdatafile(glb$numpop,resdir,ctrl$datafile)

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

datadir <- "./../../rcode2/aMSE/data-raw/"

blockE13 <- read.csv(paste0(datadir,"block13e.csv"),header=TRUE)

save(blockE13,file=paste0(datadir,"blockE13.RData"))


# check and transfer -------------------------------------------------------






tools::checkRdaFiles(paths=datadir)
tools::resaveRdaFiles(paths=datadir,compress="auto")
tools::checkRdaFiles(paths=datadir)


















