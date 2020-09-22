# An initial attempt at conditioning the model on the west coast.

# a minimal example
starttime <- as.character(Sys.time())
library(aMSE)
library(rutilsMH)
library(makehtml)  # use as needed.

datdir <- "c:/Users/User/DropBox/A_code/aMSEUse/conddata/sspmrun"
resdir <- "c:/Users/User/DropBox/A_code/aMSEUse/conddata/sspm"
dirExists(datdir,make=TRUE,verbose=TRUE)
# You now need to ensure that there is a control.csv, zone1sm\au2pop6.csv
# and region1.csv file in the data directory
# ctrlfiletemplate(resdir)
# zonefiletemplate(resdir)
# datafiletemplate(6,resdir,filename="zone1sau2pop6.csv")

ctrl <- checkresdir(datdir)
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

