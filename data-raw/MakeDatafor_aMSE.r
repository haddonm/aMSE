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
region1 <- readregionfile(datadir,ctrl$regionfile)
glb <- region1$globals
constants <- readdatafile(glb$numpop,datadir,ctrl$datafile)



save(ctrl,file=paste0(datadir,"ctrl.RData"))
save(region1,file=paste0(datadir,"region1.RData"))
save(constants,file=paste0(datadir,"constants.RData"))


# Complete region objects---------------------------------------------

library(rutilsMH)
library(aMSE)
library(microbenchmark)

datadir <- "./../../rcode2/aMSE/data-raw/"
# read data files ----------------------------------------------------
resdir <- "./../../rcode2/aMSEUse/out/run1"
dirExists(resdir,make=TRUE,verbose=TRUE)
# You now need to ensure that there is a control.csv, reg1smu2pop6.csv
# and region1.csv file in the data directory
ctrl <- checkresdir(resdir)
runname <- ctrl$runlabel
region1 <- readregionfile(resdir,ctrl$regionfile)
glb <- region1$globals
constants <- readdatafile(glb$numpop,resdir,ctrl$datafile)

out <- setupregion(constants, region1)
regionC <- out$regionC
regionD <- out$regionD
product <- out$product
glb <- out$glb
testregC <- regionC
testregD <- regionD

save(testregC,file=paste0(datadir,"testregC.RData"))
save(testregD,file=paste0(datadir,"testregD.RData"))
save(product,file=paste0(datadir,"product.RData"))


# some cpue data -----------------------------------------------------

datadir <- "./../../rcode2/aMSE/data-raw/"

blockE13 <- read.csv(paste0(datadir,"block13e.csv"),header=TRUE)

save(blockE13,file=paste0(datadir,"blockE13.RData"))


# check and transfer -------------------------------------------------------






tools::checkRdaFiles(paths=datadir)
tools::resaveRdaFiles(paths=datadir,compress="auto")
tools::checkRdaFiles(paths=datadir)


















