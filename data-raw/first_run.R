
# Used in readme.Rmd to illustrate the running of aMSE

options("show.signif.stars"=FALSE,
        "stringsAsFactors"=FALSE,
        "max.print"=50000,
        "width"=240)

# use dohistoricC  --------------------------------------------------------
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

# this should be in resdir, but during development is in data-raw
source(paste0(ddir,"aMSE/data-raw/","TasmanianHS.R"))

#data(zone)
starttime <- (Sys.time())
zone <- makeequilzone(resdir,"control2.csv") # normally would read in a file
equiltime <- (Sys.time()); print(equiltime - starttime)

glb <- zone$glb
ctrl <- zone$ctrl
zone1 <- zone$zone1
projC <- zone1$projC
condC <- zone1$condC
zoneC <- zone$zoneC
zoneD <- zone$zoneD

# condition on the historic catches
zoneDD <- dohistoricC(zoneD,zoneC,glob=glb,condC,sigR=1e-08,sigB=1e-08)
midtime <- (Sys.time())
print(equiltime - starttime)
# Illustrate productivity
propD <- getzoneprops(zoneC,zoneDD,glb,year=47)
round(propD,3)
zoneDD$harvestR[45:47,]
popdefs <- getlistvar(zone$zoneC,"popdef")
round(popdefs,2)
# plotprep(width=7, height=6,newdev=FALSE)
# parset()
# plot(propD["B0",1:16],popdefs["L50mat",],type="p",cex=1.1)

# Do the replicates ------------------------------------------------------------
midtime <- (Sys.time()); print(midtime - starttime)

cmcda <- mcdahcr(arrce=condC$histCE,hsargs=hsargs,
                 yearnames=rownames(condC$histCE),saunames=glb$saunames)
str(cmcda)
pms <- cmcda$pms
multTAC <- cmcda$multTAC

out <- prepareprojection(projC,zoneC,glb,zoneDD,ctrl,multTAC)

zoneDP <- out$zoneDP
projC <- out$projC
zoneCP <- out$zoneCP

endtime <- Sys.time(); print(endtime - midtime)

# doprojections
# str(cmcda)
begintime <- Sys.time()

zoneDP <- doTASprojections(ctrl,zoneDP,zoneCP,condC$histCE,glb,mcdahcr,hsargs)

projtime <- Sys.time()
print(projtime - begintime); print(projtime - starttime)

invar <- zoneDP$cesau
lab1 <- "Catch"
plotprep(width=8, height=8,newdev=FALSE)
parset(plots=c(4,2),byrow=FALSE)
label <- glb$saunames
for (sau in 1:8) {
  ymax <- getmax(invar[,sau,])
  plot(1:30,invar[,sau,1],type="l",lwd=1,col="grey",panel.first = grid(),
       ylim=c(0,ymax),yaxs="i",ylab=paste0(lab1,"    ",label[sau]),xlab="year")
  for (i in 1:ctrl$reps) lines(1:30,invar[,sau,i],lwd=1,col="grey")
}

str(zoneDP,max.level = 1)


arrce <- rbind(condC$histCE,zoneDP$cesau[,,1])
rownames(arrce) <- 1992:2049

pm <- mcdahcr(arrce,hsargs,yearnames=1992:2049,saunames=glb$saunames)

f <- 0.05
rge <- range(condC$histCE[,2])
rge
adjust <- (rge[2]-rge[1])*f
rge + c(-adjust,adjust)






#' @title
#'
#' @param a
#' @param b
#' @param d
#' @param f
#'
#' @return
#' @export
#'
#' @examples
extrange <- function(a,b,d,f) {


}


