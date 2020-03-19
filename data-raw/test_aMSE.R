

library(rutilsMH)
library(aMSE)
library(microbenchmark)

# read data files ----------------------------------------------------
datadir <- "./../../rcode2/aMSE/data-raw/"
ctrlfile <- "control.csv"
ctrl <- readctrlfile(datadir,ctrlfile)
reg1 <- readregionfile(datadir,ctrl$regionfile)
constants <- readdatafile(datadir,ctrl$datafile,reg1$globals)

# Define the Zone ----------------------------------------------------
ans <- makeregionC(ctrl,reg1,constants)
regionC <- ans$regionC
popdefs <- ans$popdefs
ans <- makeregion(reg1$globals,regionC)

regionC <- ans$regionC  # region constants
regionD <- ans$regionD  # region dynamics




sapply(regionC,"[[","bLML")
sapply(regionC,"[[","SaM")
str(regionD)

regionD$exploitB[1,]
round(regionD$Nt[,1,])
