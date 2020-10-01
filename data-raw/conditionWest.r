# An initial attempt at full conditioning the model on the west coast.
# this includes reading in all available conditioning data

# Generate a website------------------------------------------------------------
starttime <- as.character(Sys.time())
library(aMSE)
library(rutilsMH)
library(makehtml)  # use as needed.

datdir <- "c:/Users/User/DropBox/A_code/aMSEUse/conddata/sspmrun"
resdir <- "c:/Users/User/DropBox/A_code/aMSEUse/conddata/sspm"
dirExists(datdir,make=TRUE,verbose=TRUE)
# You now need to ensure that there is a control.csv, zone1.csv
# and region1.csv file in the data directory plus some data .csv files

ctrl <- checkctrldat(datdir,ctrlfile="control.csv")
zone1 <- readzonefile(datdir,ctrl$zonefile)
glb <- zone1$globals     # glb without the movement matrix
constants <- readdatafile(glb$numpop,datdir,ctrl$datafile)

out <- setupzone(constants,zone1) # make operating model
zoneC <- out$zoneC
zoneD <- out$zoneD
glb <- out$glb        # glb now has the movement matrix
product <- out$product     # important bits usually saved in resdir
          # did the larval dispersal level disturb the equilibrium?
regDe <- testequil(zoneC,zoneD,glb)

popd <- round(sapply(zoneC,"[[","popdef"),5)
props <- round(getzoneprops(zoneC,zoneD,glb,1),5)
rbind(c(popd["AvRec",],NA),props["MSY",])

dirExists(resdir,make=TRUE,verbose=TRUE)
resfile <- setuphtml(resdir,cleanslate = TRUE)# prepare to save and log results

plotproductivity(resdir,product,glb)
biology_plots(resdir, glb, zoneC)
numbersatsize(resdir, glb, zoneD)
plothistcatch(zone1=zone1,pops=c(1,2,3,8),resdir=resdir)
plothistcatch(zone1=zone1,pops=c(4,5,6,7),resdir=resdir)
plothistCE(zone1=zone1,pops=c(1,2,3,8),resdir=resdir)
plothistCE(zone1=zone1,pops=c(4,5,6,7),resdir=resdir)


endtime <- as.character(Sys.time())

reportlist <- list(starttime=starttime,endtime=endtime,
                   zoneC=zoneC, zoneD=zoneD, product=product,
                   glb=glb,constants=constants)

runnotes <- "This is the first conditioning on the west coast."
# If you unhash this component it will generate a local website inside
# resdir and open it so you can see the results so far.
 make_html(replist=reportlist,resdir=resdir,width=500,
          openfile=TRUE,runnotes=runnotes,verbose=FALSE,
          packagename = "aMSE",htmlname="testrun")



# Apply historical Catches -----------------------------------------------------
 library(aMSE)
 library(rutilsMH)
 library(makehtml)  # use as needed.

 datdir <- "c:/Users/User/DropBox/A_code/aMSEUse/conddata/sspmrun"
 resdir <- "c:/Users/User/DropBox/A_code/aMSEUse/conddata/sspm"
 dirExists(datdir,make=TRUE,verbose=TRUE)
 source("C:/Users/User/Dropbox/A_Code/aMSE/data-raw/candidatefunctions.R")


 # You now need to ensure that there is a control.csv, zone1.csv
 # and region1.csv file in the data directory plus some data .csv files
 ctrl <- checkctrldat(datdir,ctrlfile="control.csv")
 zone1 <- readzonefile(datdir,ctrl$zonefile)
 glb <- zone1$globals     # glb without the movement matrix
 constants <- readdatafile(glb$numpop,datdir,ctrl$datafile)

 out <- setupzone(constants,zone1,inc=0.0025) # make operating model
 zoneC <- out$zoneC
 zoneD <- out$zoneD
 glb <- out$glb        # glb now has the movement matrix
 product <- out$product     # important bits usually saved in resdir
 # did the larval dispersal level disturb the equilibrium?
 regDe <- testequil(zoneC,zoneD,glb)
 origzoneD <- zoneD
 nSAU <- glb$nSAU
 yrs <- zone1$histyr[,"year"]
 histCE <- zone1$histCE
 pickyr <- which(yrs %in% 1992:2019)
 # do the depletions

 zoneD <- origzoneD
 constants["initdepl",] <- c(0.85,0.75,0.8,0.9,0.725,0.95,0.75,0.9)
 origdepl <- constants["initdepl",] #c(0.6,0.6,0.7,0.6,0.6,0.6,0.6,0.6)
 zoneDD <- depleteSAU(zoneC,zoneD,glb,depl=origdepl,product,len=15)
 if (ctrl$hcrname == "historicalC") zoneD <- dohistC(zoneDD,zoneC,zone1,glb)
 predCE <- zoneD$cpue; rownames(predCE) <- yrs
 pCE <- predCE
 for (i in 1:nSAU) pCE[,i] <- predCE[,i]/mean(predCE[pickyr,i])
 ssq <- rep(0.0,nSAU)
 for (i in 1:nSAU) ssq[i] <- sum((log(histCE[,i]) - log(pCE[pickyr,i]))^2,na.rm=TRUE)
 rbind(constants["initdepl",],round(ssq,5))
 sum(ssq)
 plotCPUE(pCE,histCE,pickyr,1:8)



 plotSAUdepl(zoneD$deplsB)
 zprop <- getzoneprops(zoneC,zoneD,glb,year=nyrs)
 rbind(constants["initdepl",],round(zprop["SpBDepl",1:8],4))


 nSAU <- glb$nSAU
 yrs <- zone1$histyr[,"year"]
 nyrs <- length(yrs)
 predCE <- zoneD$cpue; rownames(predCE) <- yrs
 histCE <- zone1$histCE
 pickyr <- which(yrs %in% 1992:2019)
 pCE <- predCE
 for (i in 1:nSAU) pCE[,i] <- predCE[,i]/mean(predCE[pickyr,i])
 ssq <- 0.0
 for (i in 1:nSAU) ssq <- ssq + sum((histCE[,i] - pCE[pickyr,i])^2,na.rm=TRUE)
 ssq



 plotCPUE(pCE,histCE,pickyr,1:8)













