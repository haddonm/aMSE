# file reading ------------------------------------------------------

library(rutilsMH)
library(aMSE)
library(microbenchmark)
setpalette("R4")

# read data files ----------------------------------------------------
datadir <- "./../../rcode2/aMSE/data-raw/"
ctrlfile <- "control.csv"
source(filenametopath(datadir,"sourcer.R"))

ctrl <- readctrlfile(datadir,ctrlfile)
reg1 <- readregionfile(datadir,ctrl$regionfile)
glb <- reg1$globals
constants <- readdatafile(datadir,ctrl$datafile,glb)

# Define the Zone ----------------------------------------------------
source(filenametopath(datadir,"sourcer.R"))
ans <- makeregionC(reg1,constants)
regionC <- ans$regionC
popdefs <- ans$popdefs
ans <- makeregion(glb,regionC)

regionC <- ans$regionC  # region constants
regionD <- ans$regionD  # region dynamics

#str(regionC[[1]])

ans <- findunfished(regionC,regionD,glb)
regLC <- ans$regionC  # region constants
regLD <- ans$regionD  # region dynamics
product <- ans$product

str(product)

findmsy(production)

inH <- rep(0.0,glb$numpop)


# testing the equilibrium --------------------------------------------
 regDe <- testequil(regLC,regLD,glb,inH)
# regDe$matureB
# regDe$exploitB







findF1(1,results,location=FALSE)

numpop <- glb$numpop

plotprep(width=7,height=6,newdev=FALSE)
parset(plots=c(2,1))
plotprod(production,xname="MatB",xlab="Spawning Biomass")
xval <- findmsy(production)
for (pop in 1:numpop) abline(v=xval[pop,"MatB"],lwd=2,col=pop)
plotprod(production,xname="AnnH",xlab="Annual Harvest Rate")
xval <- findmsy(production)
for (pop in 1:numpop) abline(v=xval[pop,"AnnH"],lwd=2,col=pop)











