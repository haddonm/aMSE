
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
equiltime <- (Sys.time())
#print(equiltime - starttime)

# condition on the historic catches
zoneDD <- dohistoricC(zone$zoneD,zone$zoneC,glob=zone$glb,zone$zone1$condC)
midtime <- (Sys.time())
print(equiltime - starttime)


propD <- getzoneprops(zone$zoneC,zoneDD,zone$glb,year=47)
round(propD,3)

zoneDD$harvestR[45:47,]
popdefs <- getlistvar(zone$zoneC,"popdef")
round(popdefs,2)

plotprep(width=7, height=6,newdev=FALSE)
parset()
plot(propD["B0",1:16],popdefs["L50mat",],type="p",cex=1.1)



# Do the replicates ------------------------------------------------------------

midtime <- (Sys.time())
print(midtime - starttime)

glb <- zone$glb
ctrl <- zone$ctrl
zone1 <- zone$zone1
projC <- zone1$projC
condC <- zone1$condC
zoneC <- zone$zoneC

#   source(paste0(ddir,"aMSE/data-raw/","TasmanianHS.R"))


cmcda <- calibrateHCR(histCE=condC$histCE, saunames=zone1$SAUnames,
                      hsargs=hsargs,sauindex=glb$sauindex,projyrs=projC$projyrs)

str(cmcda)
pms <- cmcda$pms
multTAC <- cmcda$yrmultTAC

out <- prepareprojection(zone1,zoneC,glb,zoneDD)

zoneDP <- out$zoneDP
projC <- out$projC
zoneCP <- out$zoneCP


endtime <- Sys.time()
print(endtime - midtime)


# endcatch the final year's catches from the fishery conditioned zone.
#     These are required to generate the aspiration catches by population for
#     the first year of the projections
#
#
# develop projections-----------------------------------------------------------

zoneCP=zoneCP
zoneDP=zoneDP
glob=glb
ctrl=ctrl
projyrs=projC$projyrs
applyHS=mcdahcr
hsargs=hsargs
projpms=cmcda
# get important constants

npop <- glob$numpop
nsau <- glob$nSAU
Ncl <- glob$Nclass
yrce <- nrow(condC$histCE)
nyrs <- projyrs
endyr <- nrow(zoneDD$catch)
movem <- glob$move
reps <- ctrl$reps
sauindex <- glob$sauindex
matb <- numeric(npop)
origTAC <- sum(zoneDD$catch[yrce,]) # mean sum of catches in last year
if (ctrl$randseedP > 0) set.seed(ctrl$randseedP) # set random seed if desired
endcatch <- tapply(zoneDD$catch[endyr,],sauindex,sum,na.rm=TRUE)
acatch <- endcatch * multTAC[yrce,]
exb <- zoneDD$exploitB[endyr,]
inN <- zoneDD$Nt[,endyr,]
sigmar <- ctrl$withsigR # needed to add recruitment variation
sigmab <- ctrl$withsigB


zoneDP <- initiateHS(zoneDP,zoneCP,exb,inN,acatch,sigmar,sigmab,glb)


# doprojections
# str(cmcda)
begintime <- Sys.time()

zoneDP <- doTASprojections(ctrl,zoneDP,zoneCP,condC$histCE,glb,mcdahcr,hsargs)


projtime <- Sys.time()
print(projtime - begintime)
print(projtime - starttime)


invar <- zoneDP$matureB
lab1 <- "MatureB"
plotprep(width=8, height=8,newdev=FALSE)
parset(plots=c(4,2),byrow=FALSE)
label <- glb$saunames
for (sau in 1:8) {
  ymax <- getmax(invar[,sau,])
  plot(1:30,invar[,sau,1],type="l",lwd=1,col="grey",panel.first = grid(),
       ylim=c(0,ymax),yaxs="i",ylab=paste0(lab1,"    ",label[sau]),xlab="year")
  for (i in 1:reps) lines(1:30,invar[,sau,i],lwd=1,col="grey")
}



