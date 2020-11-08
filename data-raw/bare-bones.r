# Used in readme.Rmd to illustrate the running of aMSE

# setup run + resdir  --------------------------------------------------------
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
resdir <- paste0(ddir,"aMSEUse/conddata/generic2")
dirExists(resdir,make=TRUE,verbose=TRUE)

# equilibrium zone -------------------------------------------------------------
# You now need to ensure that there is, at least, a control.csv, zone1.csv
# and region1.csv file in the data directory plus some other data .csv files
# depending on how conditioned you want the model to be. Templates for the
# correct format can be produced using ctrlfiletemplate(), datafiletemplate(),
# and zonefiletemplate.
zone1 <- readctrlzone(resdir,infile="control.csv")
ctrl <- zone1$ctrl
glb <- zone1$globals     # glb without the movement matrix
constants <- readdatafile(glb$numpop,resdir,ctrl$datafile)
#zone1$initLML <- 140
out <- setupzone(constants,zone1) # make operating model
zoneC <- out$zoneC
zoneD <- out$zoneD
glb <- out$glb             # glb now has the movement matrix
product <- out$product     # important bits usually saved in resdir
# did the larval dispersal level disturb the equilibrium?
zoneD <- testequil(zoneC,zoneD,glb)
zoneC <- resetexB0(zoneC,zoneD) # rescale exploitB to avexplB after dynamics
initialC <- zoneC
unfishedD <- zoneD         # keep a copy of the unfished zone

equiltime <- (Sys.time())
print(equiltime - starttime)
# deplete generic zone ---------------------------------------------------------

zoneC <- initialC
zoneD <- unfishedD
# set the initial depletion levels if not defined in the data file
# the values below manage to achieve c(0.7,0.65,0.7,0.6,0.5,0.45,0.475,0.45)
zone1$condC$initdepl <- c(0.53,0.53,0.48,0.48,0.49,0.49,0.52,0.52,
                          0.51,0.51,0.45,0.45,0.484,0.484,0.467,0.467)
origdepl <-  zone1$condC$initdepl
zoneDD <- depletepop(zoneC,zoneD,glb,depl=origdepl,product,len=12)
propD <- getzoneprops(zoneC,zoneDD,glb,year=1)
#round(propD,2)

# setup for projection ---------------------------------------------------------
inityrs <- 10
projyrs <- zone1$projC$projyrs + inityrs
reps <- ctrl$reps
projC <- modprojC(zoneC,glb,zone1)
zoneC <- modzoneCSel(zoneC,projC$Sel,projC$SelWt,glb,projyrs)
zoneDR <- makezoneDR(projyrs,reps,glb,zoneDD) # zoneDReplicates
zoneDRp <- addrepvar(zoneC,zoneDR,zoneDR$harvestR,glb,ctrl)
midtime <- (Sys.time())
print(midtime - equiltime)
# prepare the HS --------------------------------------------------------------

if (projC$HS == "constantCatch") {
  hsFunc <- constCatch
  inTAC <- projC$HSdetail
}
if (projC$HS == "MCDA") {
  hsFunc <- doMCDA
  mcdafile <- projC$HSdetail
  optCE <- MCDAdata(resdir,mcdafile,zone1$SAUnames)
}


# Do the replicates ------------------------------------------------------------
inityr <- 10
saunames <- zone1$SAUnames
sauindex <- glb$sauindex
pyrs <- projC$projyrs + inityr
B0 <- tapply(sapply(zoneC,"[[","B0"),sauindex,sum)
exB0 <- tapply(sapply(zoneC,"[[","ExB0"),sauindex,sum)

zoneDP <- constCatch(1100,zoneDRp,glb,ctrl,projC$projyrs,inityrs=10)
sauzoneDP <- asSAU(zoneDP,sauindex,saunames,B0,exB0)

endtime <- (Sys.time())
print(endtime - midtime)

#calculate the relative MSY weighted MSY-depletion level for each SAU
pmsydepl <- sapply(zoneC,"[[","MSYDepl")
pmsy <- sapply(zoneC,"[[","MSY")
smsy <- tapply(pmsy,sauindex,sum)
smsydepl <- tapply((pmsydepl * pmsy) / smsy[sauindex],sauindex,sum)


ans <- sauzoneDP$saudeplsB
label <- "Mature Depletion"
uplim <- 0.7
plotprep(width=7,height=8,newdev=FALSE)
parset(plots=c(4,2),byrow=FALSE,margin=c(0.25,0.4,0.1,0.05),
       outmargin=c(1,1,0,0))
yrs <- 1:pyrs
for (sau in 1:8) {
  maxy <- uplim
  if (is.null(uplim)) maxy <- getmax(ans[,sau,])
  plot(yrs,ans[,sau,1],lwd=1,col="grey",panel.first=grid(),ylim=c(0,maxy),
       xlab="",ylab=saunames[sau])
  for (iter in 1:reps) lines(yrs,ans[,sau,iter],lwd=1,col="grey")
  abline(h=smsydepl[sau],lwd=2,col=4)
  abline(v=inityr,lwd=1,col=2)
}
mtext("Years",side=1,line=-0.1,outer=TRUE,cex=1.0,font=7)
mtext(label,side=2,line=-0.1,outer=TRUE,cex=1.0,font=7)



# plotprep(width=7,height=5)
# parset(plots=c(2,2))
# hist(colSums(zoneDRp$exploitB[1,,]),main="",col=4,xlab="Exploitable B")
# hist(colSums(zoneDRp$matureB[1,,]),main="",col=4,xlab="Mature B")
# hist(colSums(zoneDRp$recruit[1,,]/1e06),main="",col=4,xlab="Recruitment")
# hist(apply(zoneDRp$cpue[1,,],2,mean),main="",col=4,xlab="CPUE")


# equilibrium zone characterization---------------------------------------------
# resfile <- setuphtml(resdir)# prepare to save and log results
#
# plotproductivity(resdir,product,glb)
# biology_plots(resdir, glb, zoneC)
# numbersatsize(resdir, glb, zoneD)
#
# endtime <- as.character(Sys.time())
#
# reportlist <- list(starttime=starttime,endtime=endtime,
#                    zoneC=zoneC, zoneD=zoneD, product=product,
#                    glb=glb,constants=constants)
# runnotes <- "This is a bare-bones example."
# # If you unhash this component it will generate a local website inside
# # resdir and open it so you can see the results so far.
# make_html(replist=reportlist,resdir=resdir,width=500,
#           openfile=TRUE,runnotes=runnotes,verbose=FALSE,
#           packagename = "aMSE",htmlname="testrun")















