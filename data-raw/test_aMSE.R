


library(rutilsMH)
library(aMSE)
library(microbenchmark)
datadir <- "./../../rcode2/aMSE/data-raw/"

ctrlfile <- "control.csv"

ctrl <- readctrlfile(datadir,ctrlfile)



data(condDat)
str(condDat)

glb <- condDat$globals

# Define the Zone ----------------------------------------------------

out <- makeregionC(condDat)
str(out,max.level = 1)
regionC <- out$regionC
popdefs <- out$popdefs



out2 <- makeregion(glb=glb,regC=regionC)
regionC <- out2$regionC
regionD <- out2$regionD

sapply(regionC,"[[","bLML")
sapply(regionC,"[[","SaM")
str(regionD)

regionD$exploitB

