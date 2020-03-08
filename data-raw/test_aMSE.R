


library(rutilsMH)
library(aMSE)
library(microbenchmark)
datadir <- "./../../rcode2/aMSE/data-raw/"

data(condDat)
str(condDat)

glb <- condDat$globals

# Define the Zone ----------------------------------------------------
out <- makeZone(condDat,uplim=0.4)
str(out,max.level = 1)
