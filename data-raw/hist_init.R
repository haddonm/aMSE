



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

getDD <- function(x,sau) {
  rownames(x) <- 1973:2019
  colnames(x) <- sau
  return(x)
}

cpue <- getDD(zoneDD$cpue,zoneD$SAU)
catch <- getDD(zoneDD$catch,zoneD$SAU)
round(cpue,3)
round(catch,3)


plotprep(width=8,height=7)
parset(plots=c(4,4))
for (i in 1:16) {
  ymax <- getmax(cpue[,i])
  plot(catch[,i],cpue[,i],type="p",pch=1,xlab=zoneD$SAU[i],ylab="",
       ylim=c(0,650))
}




