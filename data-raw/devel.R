# file reading ------------------------------------------------------

library(rutilsMH)
library(aMSE)
library(microbenchmark)
setpalette("R4")

# read data files ----------------------------------------------------
datadir <- "./../../rcode2/aMSEUse/run2/data"
# ctrlfile <- "control.csv"
source(file.path(getwd(),"data-raw/sourcer.R"))
 ctrl <- readctrlfile(datadir)
 region1 <- readregionfile(datadir,ctrl$regionfile)
 glb <- region1$globals
 constants <- readdatafile(datadir,ctrl$datafile,glb)
#data(ctrl)
#data(region1)
#glb <- region1$globals
#data(constants)

# Define the Zone without production ---------------------------------
ans <- makeregionC(region1,constants)
regionC <- ans$regionC
popdefs <- ans$popdefs
ans <- makeregion(glb,regionC)
regionC <- ans$regionC  # region constants
regionD <- ans$regionD  # region dynamics
#  str(regionC[[1]])


# estimate production and move regionC to equilibrium-----------------
ans <- findunfished(regionC,regionD,glb)
regionC <- ans$regionC  # region constants
regionD <- ans$regionD  # region dynamics
product <- ans$product

str(regionC[[1]])
str(regionC[[4]])
str(product)

# Summarize MSY and related statistics -------------------------------
approxMSY <- findmsy(product)
approxMSY

# Some summaries -----------------------------------------------------
sapply(regionC,"[[","MSY")           # msy by population
sum(sapply(regionC,"[[","MSY"))      # total msy
sum(sapply(regionC,"[[","effB0"))    # total effective B0
sum(sapply(regionC,"[[","B0"))       # total B0
sum(sapply(regionC,"[[","effExB0"))  # total effective exploitable B0
sum(sapply(regionC,"[[","ExB0"))     # total exploitable B0



# testing the equilibrium --------------------------------------------
# This runs the dynamics for Nyrs with zero harvest rate, a constant
# larval dispersal rate (from globals), and negligible recruitment
# variation and expects the calculated components to remain constant
# within three decimal places.
 regDe <- testequil(regionC,regionD,glb)
 str(regDe)

 data("product")
 xval <- findmsy(product)
 for (pop in 1:glb$numpop) {
         regionC[[pop]]$MSY <- xval[pop,"Catch"]
         regionC[[pop]]$MSYDepl <- xval[pop,"Deplet"]
 }



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










