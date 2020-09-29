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

 # do the depletions
 origzoneD <- zoneD
 zoneD <- origzoneD
 origdepl <- c(0.6,0.6,0.7,0.6,0.6,0.6,0.6,0.6)
 zoneDD <- depleteSAU(zoneC,zoneD,glb,depl=origdepl,product,len=15)
 zprop <- getzoneprops(zoneC,zoneDD,glb,year=1)
 round(zprop,4)
 zoneD <- zoneDD

 nSAU <- glb$nSAU
 str(zone1,max.level = 1)
 if (ctrl$hcrname == "historicalC") {
   histC <- zone1$histCatch
   histCE <- zone1$histCE
   histyr <- zone1$histyr
   yrs <- histyr[,"year"]
   histLML <- histyr[,"histLML"]
   nyrs <- length(yrs)
   for (yr in 2:nyrs) {  # yr=2
     year <- yrs[yr]
     lml <- histLML[yr]
     catchpop <- histC[yr,]
     zoneD <- oneyearC(zoneC=zoneC,zoneD=zoneD,Ncl=glb$Nclass,
                       catchp=catchpop,year=yr,sigmar=1e-08,
                       npop=glb$nSAU,movem=glb$move)
   }

 }
 zpropF <- getzoneprops(zoneC,zoneD,glb,year=nyrs)
 round(zpropF,4)

 outB <- zoneD$deplsB
 plotprep(width=7,height=6,newdev = FALSE)
 parset(plots=c(3,1))
 ymax <- getmax(outB[,c(1,2,3,8)])
 plot(yrs,outB[,1],type="l",lwd=2,ylim=c(0,ymax),panel.first=grid(),ylab="6,7,8,13")
 for (i in c(2,3,8)) lines(yrs,outB[,i],lwd=2,col=i)
 legend("topright",legend=c(6,7,8,13),col=c(1,2,3,8),lwd=3,bty="n",cex=1.0)
 ymax <- getmax(outB[,c(4,5,7)])
 plot(yrs,outB[,4],type="l",lwd=2,ylim=c(0,ymax),col=4,panel.first=grid(),
      ylab=c("9,10,12"))
 for (i in c(5,7)) lines(yrs,outB[,i],lwd=2,col=i)
 legend("topright",legend=c(9,10,12),col=c(4,5,7),lwd=3,bty="n",cex=1.0)
 ymax <- getmax(outB[,6])
 plot(yrs,outB[,6],type="l",lwd=2,ylim=c(0,ymax),col=6,panel.first=grid(),
      ylab="11")
 legend("topright",legend=c(11),col=c(6),lwd=3,bty="n",cex=1.0)





 # deplete SAUs separately--------------------------------------------------

