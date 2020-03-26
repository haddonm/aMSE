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
glb <- reg1$globals
constants <- readdatafile(datadir,ctrl$datafile,glb)



save(ctrl,file=paste0(datadir,"ctrl.RData"))
save(region1,file=paste0(datadir,"region1.RData"))
save(constants,file=paste0(datadir,"constants.RData"))


# check and transfer -------------------------------------------------------



tools::checkRdaFiles(paths=datadir)
tools::resaveRdaFiles(paths=datadir,compress="auto")
tools::checkRdaFiles(paths=datadir)

devtools::use_data(condDat,
                   pkg="C:/A_mal/Rcode/Abalone/AbMSE",
                   internal=FALSE, overwrite=TRUE)

devtools::use_data(condDat3,
                   pkg="C:/A_mal/Rcode/Abalone/AbMSE",
                   internal=FALSE, overwrite=TRUE)

devtools::use_data(abdat,
                   pkg="C:/A_mal/Rcode/Abalone/abspatial",
                   internal=FALSE, overwrite=TRUE)

















