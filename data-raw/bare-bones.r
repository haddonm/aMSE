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

zone <- makeequilzone(resdir,"control.csv")
    equiltime <- (Sys.time())
    origdepl <-  c(0.30,0.31,0.29,0.32,0.30,0.31,0.29,0.32)
zoneDD <- depleteSAU(zone$zoneC,zone$zoneD,zone$glb,origdepl,zone$product,len=12)
out <- prepareprojection(zone$zone1,zoneC,zone$glb,zoneDD,zone$ctrl)
zoneDR <- out$zoneDP
projC <- out$projC
zoneCP <- out$zoneC
    midtime <- (Sys.time())

    print(equiltime - starttime)
    print(midtime - equiltime)
    propD <- getzoneprops(zoneC,zoneDD,glb,year=1)
    round(propD,3)
# prepare the HS --------------------------------------------------------------

# if (projC$HS == "constantCatch") {
#   hsFunc <- constCatch
#   inTAC <- projC$HSdetail
# }
# if (projC$HS == "MCDA") {
#   hsFunc <- doMCDA
#   mcdafile <- projC$HSdetail
#   optCE <- MCDAdata(resdir,mcdafile,zone1$SAUnames)
# }
#

# Do the replicates ------------------------------------------------------------
inityr <- zone$zone1$projC$inityrs
saunames <- zone$zone1$SAUnames
sauindex <- zone$glb$sauindex
pyrs <- projC$projyrs + inityr
B0 <- tapply(sapply(zone$zoneC,"[[","B0"),sauindex,sum)
exB0 <- tapply(sapply(zone$zoneC,"[[","ExB0"),sauindex,sum)

midtime <- (Sys.time())
Rprof()
zoneDP <- applymcda(zoneCP,zoneDRp,glb,ctrl,projC$projyrs,inityrs=10)
sauzoneDP <- asSAU(zoneDP,sauindex,saunames,B0,exB0)
Rprof(NULL)
endtime <- (Sys.time())
print(endtime - midtime)

#calculate the relative MSY weighted MSY-depletion level for each SAU
pmsydepl <- sapply(zoneC,"[[","MSYDepl")
pmsy <- sapply(zoneC,"[[","MSY")
smsy <- tapply(pmsy,sauindex,sum)
smsydepl <- tapply((pmsydepl * pmsy) / smsy[sauindex],sauindex,sum)


ans <- sauzoneDP$matB
label <- "Mature Biomass"
uplim <- NULL
plotprep(width=7,height=8,newdev=FALSE)
parset(plots=c(4,2),byrow=FALSE,margin=c(0.25,0.4,0.1,0.05),
       outmargin=c(1,1,0,0))
yrs <- 1:pyrs
for (sau in 1:8) {
 # maxy <- uplim
 # if (is.null(uplim))
    maxy <- getmax(ans[,sau,])
  plot(yrs,ans[,sau,1],lwd=1,type="l",col="grey",panel.first=grid(),ylim=c(0,maxy),
       xlab="",ylab=saunames[sau])
  for (iter in 1:50) lines(yrs,ans[,sau,iter],lwd=1,col="grey")
 # abline(h=smsydepl[sau],lwd=2,col=4)
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















