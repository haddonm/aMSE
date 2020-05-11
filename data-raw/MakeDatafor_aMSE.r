#rm(list=ls())   # Cleanup the R console if required
# Set up the run ----------


# library(r4cpue)
library(rutilsMH)
library(aMSE)
library(microbenchmark)


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

data(ctrl)
data(region1)
glb <- region1$globals
data(constants)
# Define the region
ans <- makeregionC(region1,constants)
regionC <- ans$regionC
popdefs <- ans$popdefs
ans <- makeregion(glb,regionC)
regionC <- ans$regionC  # region constants
regionD <- ans$regionD  # region dynamics
# adjust for larval dispersal
ans <- findunfished(regionC,regionD,glb)
testregC <- ans$regionC  # region constants
testregD <- ans$regionD  # region dynamics
product <- ans$product

save(testregC,file=paste0(datadir,"testregC.RData"))
save(testregD,file=paste0(datadir,"testregD.RData"))
save(product,file=paste0(datadir,"product.RData"))

# check and transfer -------------------------------------------------------



tools::checkRdaFiles(paths=datadir)
tools::resaveRdaFiles(paths=datadir,compress="auto")
tools::checkRdaFiles(paths=datadir)


















