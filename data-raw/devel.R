# file reading ------------------------------------------------------

library(rutilsMH)
library(aMSE)
library(microbenchmark)
setpalette("R4")

# read data files ----------------------------------------------------
datadir <- "./../../rcode2/aMSE/data-raw/"
# ctrlfile <- "control.csv"
source(filenametopath(datadir,"sourcer.R"))
# ctrl <- readctrlfile(datadir,ctrlfile)
# reg1 <- readregionfile(datadir,ctrl$regionfile)
# glb <- reg1$globals
# constants <- readdatafile(datadir,ctrl$datafile,glb)
data(ctrl)
data(region1)
glb <- region1$globals
data(constants)

# Define the Zone ----------------------------------------------------
ans <- makeregionC(region1,constants)
regionC <- ans$regionC
popdefs <- ans$popdefs
ans <- makeregion(glb,regionC)
regionC <- ans$regionC  # region constants
regionD <- ans$regionD  # region dynamics

#str(regionC[[1]])

ans <- findunfished(regionC,regionD,glb)
regionC <- ans$regionC  # region constants
regionD <- ans$regionD  # region dynamics
product <- ans$product

str(regionC[[1]])
str(regionC[[4]])
str(product)

approxMSY <- findmsy(product)
approxMSY

sum(sapply(regionC,"[[","effB0"))
sum(sapply(regionC,"[[","B0"))
sum(sapply(regionC,"[[","effExB0"))
sum(sapply(regionC,"[[","ExB0"))



# testing the equilibrium --------------------------------------------
 regDe <- testequil(regionC,regionD,glb,inH=rep(0.0,6))
# regDe$matureB
# regDe$exploitB



source(filenametopath(datadir,"sourcer.R"))



findF1(1,results,location=FALSE)

numpop <- glb$numpop

plotprep(width=7,height=6,newdev=FALSE)
parset(plots=c(2,1))
xval <- findmsy(product)
plotprod(product,xname="MatB",xlab="Spawning Biomass")
for (pop in 1:numpop) abline(v=xval[pop,"MatB"],lwd=2,col=pop)
plotprod(product,xname="AnnH",xlab="Annual Harvest Rate")
for (pop in 1:numpop) abline(v=xval[pop,"AnnH"],lwd=2,col=pop)










