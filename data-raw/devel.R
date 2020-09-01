# file reading ------------------------------------------------------
starttime <- as.character(Sys.time())
library(rutilsMH)
library(aMSE)
library(makehtml)
library(microbenchmark)

#Rprof()
# read data files ----------------------------------------------------
 resdir <- "./../../rcode2/aMSEUse/out/run1"
 dirExists(resdir,make=TRUE,verbose=TRUE)
 # You now need to ensure that there is a control.csv, zone1sau2pop6.csv
 # and zone1.csv file in the data directory
 ctrl <- checkresdir(resdir)
 runname <- ctrl$runlabel
 zone1 <- readzonefile(resdir,ctrl$zonefile)
 glb <- zone1$globals
 constants <- readdatafile(glb$numpop,resdir,ctrl$datafile)

 out <- setupzone(constants, zone1)
 zoneC <- out$zoneC
 zoneD <- out$zoneD
 product <- out$product
 glb <- out$glb
 str(zoneD)
# Rprof(NULL)
# outprof <- summaryRprof()
# outprof
# save the primary ojects to the resdir
resfile <- setuphtml(resdir,runname)

unfishedD <- zoneD # save a copy of the unfished state
save(ctrl,file=filenametopath(resdir,"ctrl.RData"))
save(glb,file=filenametopath(resdir,"glb.RData"))
save(product,file=filenametopath(resdir,"product.RData"))
save(zoneC,file=filenametopath(resdir,"zoneC.RData"))
save(unfishedD,file=filenametopath(resdir,"unfishedD.RData"))

# testing the equilibrium --------------------------------------------
# This runs the dynamics for Nyrs with zero harvest rate, a constant
# larval dispersal rate (from globals), and negligible recruitment
# variation and expects the calculated components to remain constant
# within three decimal places.
 zoneDe <- testequil(zoneC=zoneC,zoneD,glb)
 #  str(zoneDe)

# characterize productivity and unfished biology ---------------------
plotproductivity(resdir,runname,product,glb)
biology_plots(resdir, runname, glb, zoneC)
numbersatsize(resdir, runname, glb, zoneD)

# store the initial properties
unfishprops <- getzoneprops(zoneC=zoneC,zoneD=unfishedD,glb=glb,year=1)
filename <- filenametopath(resdir,"unfishprops.csv")
write.table(round(unfishprops,4),file = filename,sep=",")
#  or use tmp <- read.csv(file=filename,header=TRUE,row.names=1)
caption <- paste0("The unfished equilibrium properties of the ",
              "populations and zone before any initial depletion.")
logfilename(filename,resfile=resfile,"Tables",caption=caption)

ctrl$initdepl <-  0.40

if (ctrl$initdepl < 1.0) {
  zoneD <- dodepletion(zoneC=zoneC, unfishedD, glb,
                         depl=ctrl$initdepl, product)
  # store the initial properties after depletion
  initprops <- getzoneprops(zoneC=zoneC,zoneD=zoneD,glb=glb,year=1)
  filename <- filenametopath(resdir,"initprops.csv")
  write.table(round(initprops,4),file = filename,sep=",")
   #  or use tmp <- read.csv(file=filename,header=TRUE,row.names=1)
  caption <- paste0("The equilibrium properties of the populations ",
                "and zone after the initial depletion.")
  logfilename(filename,resfile=resfile,"Tables",caption=caption)
}
save(zoneD,file=filenametopath(resdir,"zoneD.RData"))

# compare the numbers-at-size for unfished and depleted
unfN <- getnas(unfishedD,yr=1,glb,zone1)
depN <- getnas(zoneD,yr=1,glb,zone1)

deplev <- round(initprops["SpBDepl","zone"],4)
filen <- compzoneN(unfN,depN,glb,yr=1,depl=deplev,LML=132,resdir=resdir)
caption <- paste0("Comparison of zoneal Numbers-at-Size before and ",
                  "after depletion.")
logfilename(filen,resfile=resfile,"NumSize",caption=caption)

endtime <- as.character(Sys.time())


 reportlist <- list(
   runname=runname,
   starttime=starttime,endtime=endtime,
   zoneC=zoneC, zoneD=zoneD, product=product,
   glb=glb,constants=constants
 )
 str(reportlist,max.level = 1)

 runnotes <- paste0("The results presented here relate to the included data-sets testzoneC, ",
                    "testzoneD, and product. They are for a zone made up of 2 SAU and 6 population. ",
                    "These results are currently under development and there are many more needed yet.")

#  source(filenametopath(sourcedir,"sourcer.R"))
 make_html(replist=reportlist,resdir=resdir,width=500,
           openfile=TRUE,runnotes=runnotes,verbose=FALSE)



# end of run ---------------------------------------------------------

# cont: characterization ---------------------------------------------

round(head(unfN,20))
round(head(depN,20))



# end of characterization --------------------------------------------
# outline a real run--------------------------------------------------


 storezoneC <- zoneC
 storezoneD <- zoneD

 zoneC <- storezoneC
 zoneD <- storezoneD
 npop <- glb$numpop
 Nc <- glb$Nclass
 nyrs <- glb$Nyrs
 larvdisp <- glb$larvdisp
 catch <- 340.0
 B0 <- getvar(zoneC,"B0") #sapply(zoneC,"[[","B0")
 totB0 <- sum(B0)
 prop <- B0/totB0
 catchpop <- catch * prop

 for (yr in 2:nyrs) {
       zoneD <- oneyearC(zoneC=zoneC,zoneD=zoneD,Ncl=Nc,
                        catchp=catchpop,year=yr,sigmar=1e-08,npop=npop,
                        movem=glb$move)
       catchpop <- catch* (zoneD$cpue[yr,]/sum(zoneD$cpue[yr,]))
 }

 zoneD$catch
 zoneD$cpue
 zoneD$harvestR
summzone <- getsauzone(zoneD)
summzone$catch
summzone$harvestR


catch <- summzone$catch
cpue <- summzone$cpue
nyrs=nrow(catch)
plotprep(width=7,height=5,newdev=FALSE)
parset(plots=c(2,1))
ymax <- getmax(cpue[,1:2])
plot(1:nyrs,cpue[,1],type="l",lwd=2,ylim=c(0,ymax),panel.first=grid())
lines(1:nyrs,cpue[,2],lwd=2,col=2)
ymax <- getmax(catch[,1:2])
plot(1:nyrs,catch[,1],type="l",lwd=2,ylim=c(0,ymax),panel.first=grid())
lines(1:nyrs,catch[,2],lwd=2,col=2)






Nt <- zoneD$Nt
mids <- glb$midpts
plotprep(width=8, height=9,newdev=FALSE)
parset(plots=c(3,2),margin=c(0.3,0.3,0.05,0.05),outmargin=c(1,1,0,0))
for (pop in 1:numpop) {
  plot(mids[5:105],Nt[5:105,1,pop],type="l",lwd=2,panel.first=grid())
  lines(mids[5:105],Nt[5:105,20,pop],lwd=2,col=3)
  lines(mids[5:105],Nt[5:105,40,pop],lwd=2,col=2)
  abline(v=132,col="grey")
  mtext(paste0("pop",pop),side=3,outer=FALSE,line=-1.1)
}

 # Some summaries ----------------------------------------------------
 getvar(zoneC,"MSY") #sapply(zoneC,"[[","MSY")           # msy by population
 sum(getvar(zoneC,"MSY"))      # total msy
 sum(getvar(zoneC,"B0"))       # total B0
 sum(getvar(zoneC,"ExB0"))     # total exploitable B0



 filen <- paste0("product_",runname,".RData")
 save(product,file=filenametopath(tabledir,filen))

 filen <- paste0("msytable_",runname,".csv")
 xval <- findmsy(product)
 write.csv(xval,file = filenametopath(tabledir,filen))

 filen <- filenametopath(tabledir,paste0("F1table_",runname,".csv"))
 xval <- findF1(product)
 write.csv(xval,file = filen)


 tmp <- read.csv(filen)

# barebones ----------------------------------------------------------

 data(constants)
 data(zone1)
 ans <- makezoneC(zone1,constants) # classical equilibrium
 zoneC <- ans$zoneC
 glb <- ans$glb
 ans <- makezone(glb,zoneC) # now make zoneD
 zoneC <- ans$zoneC  # zone constants used as zoneC in oneyearD
 zoneD <- ans$zoneD
 numpop <- glb$numpop
 harvest <- rep(0.2,numpop)
 zoneD <- aMSE:::oneyearD(zoneC=zoneC,zoneD=zoneD,Ncl=glb$Nclass,
                     inHt=harvest,year=2,sigmar=1e-06,npop=numpop,
                     movem=glb$move)
 str(zoneD)
 round(zoneD$catchN[60:105,1:5,1],1)

 # end barebones------------------------------------------------------



 datadir <- "./../../rcode2/aMSE/data-raw/"

 ab <- read.csv(file=paste0(datadir,"block13e.csv"),header=TRUE)
 yrs <- 1992:2019


qs <- quantile(ab$cpue,probs = c(0.5,0.55,0.6))
targCE <- qs[2]
lowCE <- targCE - 0.45 * targCE
upCE <- targCE + 0.45 * targCE
rgeCE <- c(lowCE,targCE,upCE)

plotprep(width=7,height=6,newdev=FALSE)
parset(plots=c(2,1))
plot(yrs,ab$cpue,type="p",pch=16,xlab="",ylab="CPUE",ylim=c(0,100),
     cex=1.0,panel.first=grid())
lines(yrs,ab$cpue,lwd=2)
abline(h=rgeCE,lwd=1,lty=c(2,1,2),col=c(2,1,3))
plot(yrs,ab$catch,type="p",pch=16,xlab="",ylab="CPUE",ylim=c(0,500),
     cex=1.0,panel.first=grid())
lines(yrs,ab$catch,lwd=2)

label <- colnames(ab)
label[8] <- "cpue"
colnames(ab) <- label

datalowSA::getlag(ab)

fish <- cbind(1992:2019,ab$catch,ab$cpue)
colnames(fish) <- c("year","catch","cpue")
rownames(fish) <- fish[,1]
fish
class(fish)


data("blockE13")
ab <- blockE13
grad1 <- getgrad1(ab$cpue)
range(grad1,na.rm=TRUE)
bounds <- round((range(grad1,na.rm=TRUE) * 1.1),2)
low <- seq(bounds[1],0.0,length=6)
high <- seq(0.0,bounds[2],length=6)
xax <- c(low[1:5],high)


plotprep(width=7,height=6,newdev=FALSE)
plot(xax,0:10,type="n",xaxt="n",xlab="Gradient 1",yaxt="n",
     ylab="Gradient1 Score")
axis(side=1,at=xax,label=xax)
axis(side=2,at=0:10,label=0:10)
abline(v=xax,lwd=1,col="grey",lty=3)
abline(h=0:10,lwd=1,col="grey",lty=3)
lines(c(-0.23,0.0),c(0,5),lwd=2)
lines(c(0.0,0.32),c(5,10),lwd=2)

model1 <- lm(0:5 ~ xax[1:6])
vars <- coef(model1)
score1 <- grad1
pickl0 <- which(grad1 <= 0)
score1[pickl0] <- grad1[pickl0]*vars[2] + vars[1]
points(grad1[pickl0],score1[pickl0],pch=16,cex=1.25,col=2)
model2 <- lm(5:10 ~ xax[6:11])
vars2 <- coef(model2)
pickg0 <- which(grad1 > 0)
score1[pickg0] <- grad1[pickg0]*vars2[2] + vars2[1]
points(grad1[pickg0],score1[pickg0],pch=16,cex=1.25,col=2)

cbind(1993:2019,grad1,score1)





grad4 <- getgrad4(grad1)
range(grad4)

plotprep(width=7,height=5,newdev=FALSE)
parset(plots=c(2,1))
plot(ab$year,ab$cpue,type="l",lwd=2,ylim=c(0,100))
plot(0:27,c(NA,grad1),type="l",lwd=2)
lines(0:27,c(NA,NA,NA,NA,grad4),col=2,lwd=2)
abline(h=0.0,col=3)

# targCE -------------------------------------------------------------
data(blockE13)
ab <- blockE13
qs <- quantile(ab$cpue,probs = c(0.5,0.55,0.6))
targCE <- qs[2]
lowCE <- targCE - 0.45 * targCE
upCE <- targCE + 0.45 * targCE
rgeCE <- c(lowCE,targCE,upCE)

ymax <- getmax(ab$cpue)
plotprep(width=7,height=6,newdev=FALSE)
parset(plots=c(2,1))
plot(ab$year,ab$cpue,type="p",pch=16,cex=1.25,panel.first=grid(),
     ylab="CPUE Block13 E",xlab="",ylim=c(0,ymax))
lines(ab$year,ab$cpue,lwd=1)
abline(h=targCE,col=2)
ymax <- getmax(ab$catch)
plot(ab$year,ab$catch,type="l",lwd=2,panel.first=grid(),ylim=c(0,ymax))



# summarize zoneprod ------------------------------------------------------











