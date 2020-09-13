# file reading ------------------------------------------------------
starttime <- as.character(Sys.time())
library(rutilsMH)
library(aMSE)
library(makehtml)
library(microbenchmark)

#Rprof()
# read data files ----------------------------------------------------
 resdir <- "./../../A_code/aMSEUse/out/run1"
 dirExists(resdir,make=TRUE,verbose=TRUE)
 # You now need to ensure that there is a control.csv, zone1sau2pop6.csv
 # and zone1.csv file in the data directory
 ctrlfiletemplate(resdir)
 zonefiletemplate(resdir)
 datafiletemplate(6,resdir,filename="zone1sau2pop6.csv")

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


# HISTORICAL catch data ------------------------------------------------------
library(rutilsMH)
library(aMSE)

filen <- "C:/Users/User/Dropbox/A_Code/aMSE/data-raw/WZ_cbb.csv"

condat <- read_conddata(filen)
head(condat,12)

plotprep(width=7, height=9, newdev=FALSE)
parset(plots=c(4,2),margin = c(0.25,0.25,0.05,0.05),outmargin = c(1,1,0,0))
yrs <- condat[,"year"]
label <- colnames(condat)
for (sau in 2:9) {
  ymax <- getmax(condat[,sau])
  plot(yrs,condat[,sau],type="l",lwd=2,xlab="",ylab="",ylim=c(0,ymax),
       panel.first=grid())
  mtext(label[sau],side=3,line=-1.2,outer=FALSE,cex=1.25)
}
mtext("Year",side=1,outer=TRUE,cex=1.0)
mtext("Catches (t)",side=2,outer=TRUE,cex=1.0)



plotprep(width=7, height=4, newdev=FALSE)
parset()
yrs <- condat[,"year"]
catches <- rowSums(condat[,2:9],na.rm=TRUE)
ymax <- getmax(catches)
plot(yrs,catches,type="l",lwd=2,xlab="Year",ylab="Catches (t)",ylim=c(0,ymax),
     panel.first=grid())
abline(v=1991.5,col=3)

# l ------------------------------------------------------------------------
# SIZE Composition Data------------------------------------------------------

library(rutilsMH)
library(aMSE)
library(makehtml)
library(microbenchmark)

datadir <- "C:/Users/User/Dropbox/A_Code/aMSE/data-raw/"

comp <- readRDS(paste0(datadir,"compiledMM.df.final.rds"))

properties(comp)

blks <- unique(comp$blockno)

# clean up blockno -----------------------------------------------------------
pick <- which(comp$blockno == "3,1")
comp$blockno[pick] <- "1,3"
pick <- which(comp$blockno == "48,5")
comp$blockno[pick] <- "5,48"
pick <- which(comp$blockno == "7,6")
comp$blockno[pick] <- "6,7"
pick <- which(comp$blockno == "12,9")
comp$blockno[pick] <- "9,12"
pick <- which(comp$blockno == "10,7")
comp$blockno[pick] <- "7,10"
pick <- which(comp$blockno == "10,9")
comp$blockno[pick] <- "9,10"
pick <- which(comp$blockno == "11,10")
comp$blockno[pick] <- "10,11"
pick <- which(comp$blockno == "12,10")
comp$blockno[pick] <- "10,12"
pick <- which(comp$blockno == "12,11")
comp$blockno[pick] <- "11,12"
pick <- which(comp$blockno == "13,12")
comp$blockno[pick] <- "12,13"
pick <- which(comp$blockno == "14,13")
comp$blockno[pick] <- "13,14"
pick <- which(comp$blockno == "16,14")
comp$blockno[pick] <- "14,16"
pick <- which(comp$blockno == "19,17")
comp$blockno[pick] <- "17,19"
pick <- which(comp$blockno == "20,17")
comp$blockno[pick] <- "17,20"
pick <- which(comp$blockno == "23,22")
comp$blockno[pick] <- "22,23"
pick <- which(comp$blockno == "24,23")
comp$blockno[pick] <- "23,24"
pick <- which(comp$blockno == "27,25")
comp$blockno[pick] <- "25,27"
pick <- which(comp$blockno == "39,31")
comp$blockno[pick] <- "31,39"
pick <- which(comp$blockno == "38,37")
comp$blockno[pick] <- "37,38"
pick <- which(comp$blockno == "49,48")
comp$blockno[pick] <- "48,49"

pick <- which(comp$blockno == "12,7,6")
comp$blockno[pick] <- "6,7,12"
pick <- which(comp$blockno == "12,10,11")
comp$blockno[pick] <- "10,11,12"
pick <- which(comp$blockno == "11,10,12")
comp$blockno[pick] <- "10,11,12"
pick <- which(comp$blockno == "16,13,14")
comp$blockno[pick] <- "13,14,16"
pick <- which(comp$blockno == "13,16,14")
comp$blockno[pick] <- "13,14,16"

# add month to comp ------------------------------------------------------------
comp$month <- NA
comp$month <- substr(as.character(comp[,"msr.date"]),6,7)
pick <- which(is.na(comp$month))
length(pick)
nodate <- droplevels(comp[pick,])
properties(nodate)
table(nodate$blockno,nodate$fishyear)

blks <- unique(comp$blockno)
sort(blks)

sort(blks[which(nchar(blks) > 3)])

getblkno <- function(txt) {
  return(as.numeric(unlist(strsplit(txt,","))))
}

getblkno(blks)

# get the western zone --------------------------------------------------------
pick <- which(comp$newzone == "W")
wz <- droplevels(comp[pick,])

blkn <- table(wz$blockno)
table(wz$blockno,wz$fishyear)
table(wz$blockno,wz$month)

pick <- which(wz$blockno %in% c("6","7","8","9","10","11","12","13"))
wzl <- droplevels(wz[pick,])
nrow(wz)-nrow(wzl)

table(wzl$blockno,wzl$fishyear)
table(wzl$blockno,wzl$month)

picky <- which(wzl$blockno == "12")
table(wzl$fishyear[picky],wzl$month[picky])



pick <- which((wzl$blockno == "12") & (wzl$fishyear == 2012))
length(pick)
table(wzl[picky,"month"])
byr <- droplevels(wzl[pick,])


plotprep(width=9,height=7,newdev=FALSE)
parset(plots=c(3,4),margin=c(0.3,0.3,0.05,0.05),outmargin = c(1,1,0.05,0.05),
       byrow=FALSE)
for (mth in 1:12) { # mth=1
  pickm <- which(as.numeric(byr$month) == mth)
  hist(byr[pickm,"shell.length"],col=2,breaks=seq(130,225,5),main="",
       xlab="",ylab="")
  label <- paste0(mth,"_",length(pickm))
  mtext(label,side=3,line=-1.1,adj=0.5,cex=1.1)
}
mtext("Shell Length (mm)",side=1,outer=TRUE,cex=1.1)
mtext("Frequency",side=2,outer=TRUE,cex=1.1)


doplot <- function(byr,season,bins,label) { # season=autumn
  pickm <- which(as.numeric(byr$month) %in% season)
  if (length(pickm) > 0) {
    out1 <- hist(byr[pickm,"shell.length"],col=2,breaks=bins,main="",
                 xlab="",ylab="")
    mlabel <- paste0(label,"_",length(pickm))
    mtext(mlabel,side=3,line=-1.1,adj=1,cex=1.1)
    meanl <- mean(byr[pickm,"shell.length"],na.rm=TRUE)
    sdl <- sd(byr[pickm,"shell.length"],na.rm=TRUE)
    mtext(round(meanl,2),side=3,line=-2.5,adj=1,cex=1.1)
    mtext(round(sdl,2),side=3,line=-4,adj=1,cex=1.1)
  } else {
    meanl <- NULL
    sdl <- NULL
    plotnull("No Data")
  }
  return(invisible(list(out1=out1,meanl=meanl,sdl=sdl)))
}

summer <- c(1,2,12)
spring <- c(9,10,11)
winter <- c(6,7,8)
autumn <- c(3,4,5)
bins=seq(100,252,2)

pick <- which((wzl$blockno == "12") & (wzl$fishyear == 2020))
length(pick)
table(wzl[pick,"month"])
byr <- droplevels(wzl[pick,])
range(byr$shell.length)

plotprep(width=5,height=8,newdev=FALSE)
parset(plots=c(4,1),margin=c(0.3,0.3,0.05,0.05),outmargin = c(1,1,0.05,0.05),
       byrow=FALSE)
out1 <- doplot(byr,summer,bins=bins,"summer")
out2 <- doplot(byr,autumn,bins=bins,"autumn")
out3 <- doplot(byr,winter,bins=bins,"winter")
out4 <- doplot(byr,spring,bins=bins,"spring")








b12 <- droplevels(wz[pick,])
model <- lm(log(b12$whole.weight) ~ log(b12$shell.length), data=b12)
abline(model,lwd=2,col=4)
summary(model)
coef(model)
x <- 138:195
y <- exp(coef(model))[1] * (x ^ coef(model)[2])



plotprep(width=7,height=5,newdev=FALSE)
plot1((b12$shell.length),(b12$whole.weight),type="p",pch=1,xlab="Shell Length mm",
      ylab="Total Weight",col=(b12$fishyear - 2018),defpar=FALSE)
lines(x,y,lwd=2,col=4)


# which records have weights?

pick <- which(wz$whole.weight > 0)
length(pick)

table(wz[pick,"blockno"])



pick <- which((wz$blockno == "6"))# & (wz$whole.weight > 0) &
                # (wz$whole.weight < 2000) & (wz$shell.length > 136) &
                # (wz$whole.weight > 380) & (wz$fishyear == 2020))
length(pick)
table(wz[pick,"msr.date"])





plotprep(width=7,height=5,newdev=FALSE)
plot1(wz[pick,"shell.length"],wz[pick,"whole.weight"],type="p",pch=1,
      xlab="Shell Length mm",
      ylab="Total Weight",col=(b12$fishyear - 2018),defpar=FALSE)















