



starttime <- (Sys.time())
library(aMSE)
library(rutilsMH)
library(makehtml)
library(knitr)
# Obviously you should modify the resdir to suit your own computer
if (dir.exists("c:/Users/User/DropBox")) {
  ddir <- "c:/Users/User/DropBox/A_code/"
} else {
  ddir <- "c:/Users/Malcolm/DropBox/A_code/"
}
#resdir <- paste0(ddir,"aMSEUse/conddata/generic2")
resdir <- paste0(ddir,"aMSEUse/conddata/generic")
dirExists(resdir,make=TRUE,verbose=TRUE)
#data(zone)
zone <- makeequilzone(resdir,"control2.csv") # normally would read in a file
equiltime <- (Sys.time())
nyrs <- zone$glb$Nyrs
zoneC <- zone$zoneC; zoneD <- zone$zoneD
glb <- zone$glb; zone1 <- zone$zone1
zoneDD <- dohistoricC(zoneD,zoneC,glob=glb,zone1)
x <- getzoneprops(zoneC,zoneDD,glb,year=47)
round(x[c(1,2,6,9),],3)

condce <- zone1$condC$histCE; round(condce,2)


cpue <- zoneDD$cpue
round(cpue,3)


catchsau <- zone$zone1$condC$histCatch
zoneD <- zone$zoneD
glb <- zone$glb
for (yr in 2:nyrs)
  zoneD <- oneyearC(zoneC=zone$zoneC,zoneD=zoneD,Ncl=glb$Nclass,
                  catchp=catchsau[yr,],year=yr,sigmar=1e-08,
                  npop=glb$numpop,movem=glb$move)





