

library(rutilsMH)
library(aMSE)
library(microbenchmark)

# read data files ----------------------------------------------------
rundir <- "./../../rcode2/aMSEuse/run2"
datadir <- file.path(rundir,"data")
ctrlfile <- "control.csv"
ctrl <- readctrlfile(datadir,ctrlfile)
reg1 <- readregionfile(datadir,ctrl$regionfile)
glb <- reg1$globals
constants <- readdatafile(datadir,ctrl$datafile,glb)

# Define the Zone ----------------------------------------------------
ans <- makeregionC(reg1,constants)
regionC <- ans$regionC
popdefs <- ans$popdefs
ans <- makeregion(glb,regionC)

regionC <- ans$regionC  # region constants
regionD <- ans$regionD  # region dynamics

ans2 <- findunfished(regionC,regionD,glb)
str(ans2,max.level=2)

regionC <- ans2$regionC  # region constants
regionD <- ans2$regionD  # region dynamics
product <- ans2$product

sapply(regionC,"[[","MSY")
sapply(regionC,"[[","SaM")
str(regionD)

regionD$exploitB[1,]
round(regionD$Nt[,1,])
