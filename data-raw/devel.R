# file reading ------------------------------------------------------
starttime <- as.character(Sys.time())
library(rutilsMH)
library(aMSE)
library(microbenchmark)

#Rprof()
# read data files ----------------------------------------------------
 resdir <- "./../../rcode2/aMSEUse/out/run1"
 dirExists(resdir,make=TRUE,verbose=TRUE)
 # You now need to ensure that there is a control.csv, reg1smu2pop6.csv
 # and region1.csv file in the data directory
 ctrl <- checkresdir(resdir)
 runname <- ctrl$runlabel
 region1 <- readregionfile(resdir,ctrl$regionfile)
 glb <- region1$globals
 constants <- readdatafile(glb$numpop,resdir,ctrl$datafile)

 out <- setupregion(constants, region1)
 regionC <- out$regionC
 regionD <- out$regionD
 product <- out$product
 glb <- out$glb
 str(regionD)
# Rprof(NULL)
# outprof <- summaryRprof()
# outprof
# save the primary ojects to the resdir
resfile <- setuphtml(resdir,runname)

unfishedD <- regionD # save a copy of the unfished state
save(ctrl,file=filenametopath(resdir,"ctrl.RData"))
save(glb,file=filenametopath(resdir,"glb.RData"))
save(product,file=filenametopath(resdir,"product.RData"))
save(regionC,file=filenametopath(resdir,"regionC.RData"))
save(unfishedD,file=filenametopath(resdir,"unfishedD.RData"))

# testing the equilibrium --------------------------------------------
# This runs the dynamics for Nyrs with zero harvest rate, a constant
# larval dispersal rate (from globals), and negligible recruitment
# variation and expects the calculated components to remain constant
# within three decimal places.
 regDe <- testequil(regionC=regionC,regionD,glb)
 #  str(regDe)

# characterize productivity and unfished biology ---------------------
plotproductivity(resdir,runname,product,glb)
biology_plots(resdir, runname, glb, regionC)
numbersatsize(resdir, runname, glb, regionD)

# store the initial properties
unfishprops <- getregionprops(regC=regionC,regD=unfishedD,glb=glb,year=1)
filename <- filenametopath(resdir,"unfishprops.csv")
write.table(round(unfishprops,4),file = filename,sep=",")
#  or use tmp <- read.csv(file=filename,header=TRUE,row.names=1)
caption <- paste0("The unfished equilibrium properties of the ",
              "populations and region before any initial depletion.")
logfilename(filename,resfile=resfile,"Tables",caption=caption)

ctrl$initdepl <-  0.40

if (ctrl$initdepl < 1.0) {
  regionD <- dodepletion(regionC=regionC, unfishedD, glb,
                         depl=ctrl$initdepl, product)
  # store the initial properties after depletion
  initprops <- getregionprops(regC=regionC,regD=regionD,glb=glb,year=1)
  filename <- filenametopath(resdir,"initprops.csv")
  write.table(round(initprops,4),file = filename,sep=",")
   #  or use tmp <- read.csv(file=filename,header=TRUE,row.names=1)
  caption <- paste0("The equilibrium properties of the populations ",
                "and region after the initial depletion.")
  logfilename(filename,resfile=resfile,"Tables",caption=caption)
}
save(regionD,file=filenametopath(resdir,"regionD.RData"))

# compare the numbers-at-size for unfished and depleted
unfN <- getnas(unfishedD,yr=1,glb,region1)
depN <- getnas(regionD,yr=1,glb,region1)

deplev <- round(initprops["SpBDepl","region"],4)
filen <- compregionN(unfN,depN,glb,yr=1,depl=deplev,LML=132,resdir=resdir)
caption <- paste0("Comparison of regional Numbers-at-Size before and ",
                  "after depletion.")
logfilename(filen,resfile=resfile,"NumSize",caption=caption)

endtime <- as.character(Sys.time())


 reportlist <- list(
   runname=runname,
   starttime=starttime,endtime=endtime,
   regionC=regionC, regionD=regionD, product=product,
   glb=glb,constants=constants
 )
 str(reportlist,max.level = 1)

 runnotes <- paste0("The results presented here relate to the included data-sets testregC, ",
                    "testregD, and product. They are for a region made up of 2 SMU and 6 population. ",
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


 storeregC <- regionC
 storeregD <- regionD

 regionC <- storeregC
 regionD <- storeregD
 npop <- glb$numpop
 Nc <- glb$Nclass
 nyrs <- glb$Nyrs
 larvdisp <- glb$larvdisp
 catch <- 340.0
 B0 <- getvar(regionC,"B0") #sapply(regionC,"[[","B0")
 totB0 <- sum(B0)
 prop <- B0/totB0
 catchpop <- catch * prop

 for (yr in 2:nyrs) {
       regionD <- oneyearC(regC=regionC,regD=regionD,Ncl=Nc,
                        catchp=catchpop,year=yr,sigmar=1e-08,npop=npop,
                        movem=glb$move)
       catchpop <- catch* (regionD$cpue[yr,]/sum(regionD$cpue[yr,]))
 }

 regionD$catch
 regionD$cpue
 regionD$harvestR
summreg <- getsmureg(regionD)
summreg$catch
summreg$harvestR


catch <- summreg$catch
cpue <- summreg$cpue
nyrs=nrow(catch)
plotprep(width=7,height=5,newdev=FALSE)
parset(plots=c(2,1))
ymax <- getmax(cpue[,1:2])
plot(1:nyrs,cpue[,1],type="l",lwd=2,ylim=c(0,ymax),panel.first=grid())
lines(1:nyrs,cpue[,2],lwd=2,col=2)
ymax <- getmax(catch[,1:2])
plot(1:nyrs,catch[,1],type="l",lwd=2,ylim=c(0,ymax),panel.first=grid())
lines(1:nyrs,catch[,2],lwd=2,col=2)






Nt <- regionD$Nt
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
 getvar(regionC,"MSY") #sapply(regionC,"[[","MSY")           # msy by population
 sum(getvar(regionC,"MSY"))      # total msy
 sum(getvar(regionC,"B0"))       # total B0
 sum(getvar(regionC,"ExB0"))     # total exploitable B0



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
 data(region1)
 ans <- makeregionC(region1,constants) # classical equilibrium
 regionC <- ans$regionC
 glb <- ans$glb
 ans <- makeregion(glb,regionC) # now make regionD
 regionC <- ans$regionC  # region constants used as regC in oneyearD
 regionD <- ans$regionD
 numpop <- glb$numpop
 harvest <- rep(0.2,numpop)
 regionD <- aMSE:::oneyearD(regionC=regionC,regD=regionD,Ncl=glb$Nclass,
                     inHt=harvest,year=2,sigmar=1e-06,npop=numpop,
                     movem=glb$move)
 str(regionD)
 round(regionD$catchN[60:105,1:5,1],1)

 # end barebones------------------------------------------------------



 datadir <- "./../../rcode2/aMSE/data-raw/"

 ab <- read.csv(file=paste0(datadir,"block13e.csv"),header=TRUE)
 yrs <- 1992:2019
 qs <- quantile(ab$scaledgeo,probs = c(0.5,0.55,0.6))
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


grad1 <- grad1PM(ab$cpue)
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

















