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

 ctrl <- checkctrldat(resdir)
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


# HISTORICAL catch/cpue data ---------------------------------------------------
#library(rutilsMH)
library(MQMF)

filen <- "C:/Users/User/Dropbox/A_Code/aMSE/data-raw/WZ_cbb.csv"

condat <- aMSE::read_conddata(filen)
str(condat)
catches <- condat$catches
cpue <- condat$cpue

plotprep(width=7, height=9, newdev=FALSE)
parset(plots=c(4,2),margin = c(0.25,0.25,0.05,0.05),outmargin = c(1,1,0,0))
yrs <- catches[,"year"]
label <- colnames(catches)
for (sau in 2:9) {
  ymax <- getmax(catches[,sau])
  plot(yrs,catches[,sau],type="l",lwd=2,xlab="",ylab="",ylim=c(0,ymax),
       panel.first=grid())
  mtext(label[sau],side=3,line=-1.2,outer=FALSE,cex=1.25)
}
mtext("Year",side=1,outer=TRUE,cex=1.0)
mtext("Catches (t)",side=2,outer=TRUE,cex=1.0)

plotprep(width=7, height=9, newdev=FALSE)
parset(plots=c(4,2),margin = c(0.25,0.25,0.05,0.05),outmargin = c(1,1,0,0))
yrs <- cpue[,"year"]
label <- colnames(cpue)
for (sau in 2:9) {
  ymax <- getmax(cpue[,sau])
  plot(yrs,cpue[,sau],type="l",lwd=2,xlab="",ylab="",ylim=c(0,ymax),
       panel.first=grid())
  mtext(label[sau],side=3,line=-1.2,outer=FALSE,cex=1.25)
}
mtext("Year",side=1,outer=TRUE,cex=1.0)
mtext("CPUE (t)",side=2,outer=TRUE,cex=1.0)


plotprep(width=7, height=4, newdev=FALSE)
parset()
yrs <- condat[,"year"]
catches <- rowSums(condat[,2:9],na.rm=TRUE)
ymax <- getmax(catches)
plot(yrs,catches,type="l",lwd=2,xlab="Year",ylab="Catches (t)",ylim=c(0,ymax),
     panel.first=grid())
abline(v=1991.5,col=3)

# SPM to WZ data -------------------------------------------------------------
library(MQMF)

filen <- "C:/Users/User/Dropbox/A_Code/aMSE/data-raw/WZ_cbb.csv"

condat <- aMSE::read_conddata(filen)
str(condat)
catches <- condat$catches
cpue <- condat$cpue

yrs <- catches[,"year"]
nyrs <- length(yrs)
sau <- 9
catch <- catches[,sau]
blocks <- colnames(catches)
pickyr <- match(cpue[,"year"],yrs)
addNA <- rep(NA,nyrs-length(pickyr))
fish <- cbind(yrs,catch,c(addNA,cpue[,sau]))
colnames(fish) <- c("year","catch","cpue")
blocks[sau]
getlag(fish)
#fish <- fish[-1,]
fish

param <- log(c(0.35,700,550,0.125))
plotprep(width=7,height=6,newdev=FALSE)

ans <- plotspmmod(inp=param,indat=fish,schaefer=TRUE,
                  addrmse=TRUE,plotprod=FALSE)
negLL(param,simpspm,log(fish[,"cpue"]),indat=fish)

answer <- fitSPM(pars=param,fish=fish,schaefer=TRUE,steptol = 1e-06)
outfit(answer,title = blocks[sau])
ans <- plotspmmod(inp=answer$estimate,indat=fish,schaefer=TRUE,
                  addrmse=TRUE,plotprod=FALSE)



# l ------------------------------------------------------------------------
# SIZE Composition Data------------------------------------------------------

library(rutilsMH)
library(aMSE)
library(makehtml)
library(microbenchmark)

datadir <- "C:/Users/User/Dropbox/A_Code/aMSEUse/conddata/"


pcomp <- readRDS(paste0(datadir,"compiledMM.df.final.rds"))

properties(pcomp)
novel <- c("docket","date","procnum","procname","numprocs","proclist","newzone",
           "numdays","daylist","daylist_max","datediff","numblocks","blocklist",
           "numsubblocks","subblocklist","catch","length","n","meanSL","minSL",
           "species","weight","datasource","block","subblock","year",
           "blocklist_1","blocklist_2","blocklist_3","blocklist_4","blocklist_5",
           "id","region1","region2","region3","region4","region5","same.region")
colnames(pcomp) <- novel
properties(pcomp)
pickcols <- c(1,2,3,7,8,9,12,16,17,18,21,22,24,26)
comp <- droplevels(pcomp[,pickcols])
properties(comp)

blks <- unique(comp$block)
nzone <- unique(comp$newzone)
# clean up blockno -----------------------------------------------------------
changeblk <- function(compdat,inb,outb) {
  pick <- which(compdat$block == inb)
  compdat$block[pick] <- outb
  return(compdat)
}

sort(unique(comp$block))

comp <- changeblk(comp,"3,1","1,3")
comp <- changeblk(comp,"48,5","5,48")
comp <- changeblk(comp,"7,6","6,7")
comp <- changeblk(comp,"12,9","9,12")
comp <- changeblk(comp,"10,7","7,10")
comp <- changeblk(comp,"10,9","9,10")
comp <- changeblk(comp,"11,10","10,11")
comp <- changeblk(comp,"12,10","10,12")
comp <- changeblk(comp,"12,11","11,12")
comp <- changeblk(comp,"13,12","12,13")
comp <- changeblk(comp,"14,13","13,14")
comp <- changeblk(comp,"16,14","14,16")
comp <- changeblk(comp,"19,17","17,19")
comp <- changeblk(comp,"20,17","17,20")
comp <- changeblk(comp,"23,20","20,23")
comp <- changeblk(comp,"24,23","23,24")
comp <- changeblk(comp,"27,25","25,27")
comp <- changeblk(comp,"39,31","31,39")
comp <- changeblk(comp,"38,37","37,38")
comp <- changeblk(comp,"49,48","48,49")
comp <- changeblk(comp,"12,7,6","6,7,12")
comp <- changeblk(comp,"12,10,11","10,11,12")
comp <- changeblk(comp,"11,10,12","10,11,12")
comp <- changeblk(comp,"16,13,14","13,14,16")
comp <- changeblk(comp,"13,16,14","13,14,16")

sort(unique(comp$block))


# add month to comp ------------------------------------------------------------
comp$month <- NA
comp$month <- substr(as.character(comp[,"date"]),6,7)
pick <- which(is.na(comp$month))
length(pick)
nodate <- droplevels(comp[pick,])
properties(nodate)
table(nodate$block,nodate$year)


blks <- unique(comp$block)
sort(blks)

sort(blks[which(nchar(blks) > 3)])

getblkno <- function(txt) {
  return(as.numeric(unlist(strsplit(txt,","))))
}

getblkno(blks)

# get the western zone --------------------------------------------------------
pick <- which(comp$newzone == "W")
wz <- droplevels(comp[pick,])

blkn <- table(wz$block)
table(wz$block,wz$year)
table(wz$block,wz$month)

pick <- which(wz$block %in% c("6","7","8","9","10","11","12","13"))
wzl <- droplevels(wz[pick,])
nrow(wz)-nrow(wzl)  # number of records lost
pickE <- which((wzl[,"length"] < 120) | (wzl[,"length"] > 230))
extreme <- droplevels(wzl[pickE,])
wzl <- wzl[-pickE,]

table(wzl$year,wzl$block)
table(wzl$block,wzl$month)
table(wzl$year,wzl$month)

picky <- which(wzl$block == "12")
table(wzl$year[picky],wzl$month[picky])



pick <- which((wzl$block == "12") & (wzl$year == 2012))
length(pick)
table(wzl[picky,"month"])
byr <- droplevels(wzl[pick,])


plotprep(width=9,height=7,newdev=FALSE)
parset(plots=c(3,4),margin=c(0.3,0.3,0.05,0.05),outmargin = c(1,1,0.05,0.05),
       byrow=FALSE)
for (mth in 1:12) { # mth=1
  pickm <- which(as.numeric(byr$month) == mth)
  hist(byr[pickm,"length"],col=2,breaks=seq(130,225,5),main="",
       xlab="",ylab="")
  label <- paste0(mth,"_",length(pickm))
  mtext(label,side=3,line=-1.1,adj=0.5,cex=1.1)
}
mtext("Shell Length (mm)",side=1,outer=TRUE,cex=1.1)
mtext("Frequency",side=2,outer=TRUE,cex=1.1)


doplot <- function(byr,season,bins,label) { # season=summer
  pickm <- which(as.numeric(byr$month) %in% season)
  if (length(pickm) > 0) {
    out1 <- hist(byr[pickm,"length"],col=2,breaks=bins,main="",
                 xlab="",ylab="")
    mlabel <- paste0(label,"_",length(pickm))
    mtext(mlabel,side=3,line=-1.1,adj=1,cex=1.1)
    meanl <- mean(byr[pickm,"length"],na.rm=TRUE)
    sdl <- sd(byr[pickm,"length"],na.rm=TRUE)
    mtext(round(meanl,2),side=3,line=-2.5,adj=1,cex=1.1)
    mtext(round(sdl,2),side=3,line=-4,adj=1,cex=1.1)
    ans <- list(out1=out1,meanl=meanl,sdl=sdl)
  } else {
    plotnull("No Data")
    ans <- list(out1=NULL,meanl=NULL,sdl=NULL)
  }
  return(invisible(ans))
}

summer <- c(1,2,12)
spring <- c(9,10,11)
winter <- c(6,7,8)
autumn <- c(3,4,5)
all <- c(1:12)
bins=seq(120,230,2)


plotcomp <- function(blk,year,bins,filen="",resdir="") {
  pick <- which((wzl$block == blk) & (wzl$year == year))
  records <- length(pick)
  table(wzl[pick,"month"])
  byr <- droplevels(wzl[pick,])
  if(records > 0) {
    if ((nchar(filen) > 0) & (nchar(resdir) > 0)) {
      filename <- filenametopath(resdir,filen)
    } else { filename <-  "" }
    plotprep(width=7,height=9,newdev=FALSE,filename=filename,verbose=FALSE)
    parset(plots=c(5,1),margin=c(0.3,0.3,0.05,0.05),outmargin = c(1,1,0.05,0.05),
           byrow=FALSE)
    out1 <- doplot(byr,summer,bins=bins,"summer")
    out2 <- doplot(byr,autumn,bins=bins,"autumn")
    out3 <- doplot(byr,winter,bins=bins,"winter")
    out4 <- doplot(byr,spring,bins=bins,"spring")
    out5 <- doplot(byr,all,bins=bins,"year")
    mtext(paste0("Shell Length - Block ",blk," year ",year),side=1,outer=TRUE,
          line=-0.2,cex=1.0)
    if (nchar(filename) > 0) {
      addlab <- paste0(blk,"_",year)
      caption <- paste0("Size composition for block-year ",addlab)
      addplot(filename,resdir=resdir,category=blk,caption)
    }
    ans <- list(summer=out1,autumn=out2,winter=out3,spring=out4,all=out5)
    return(invisible(ans))
  }
} # end of plotcomp

plotyearcomp <- function(dat,blk,year,bins,filen="",resdir="") {
  pick <- which((dat$block == blk) & (dat$year == year))
  records <- length(pick)
  byr <- droplevels(dat[pick,])
  if(records > 0) {
    if ((nchar(filen) > 0) & (nchar(resdir) > 0)) {
      filename <- filenametopath(resdir,filen)
    } else { filename <-  "" }
    plotprep(width=7,height=4,newdev=FALSE,filename=filename,verbose=FALSE)
    parset(plots=c(1,1))
    out5 <- doplot(byr,season=c(1:12),bins=bins,"year")
    mtext(paste0("Shell Length - Block ",blk," year ",year),side=1,outer=TRUE,
          line=-0.2,cex=1.0)
    if (nchar(filename) > 0) {
      addlab <- paste0(blk,"_",year)
      caption <- paste0("Size composition for block-year ",addlab)
      addplot(filename,resdir=resdir,category=blk,caption)
    }
    return(invisible(out5))
  }
} # end of plotyearcomp




resdir <- "C:/Users/User/Dropbox/A_Code/aMSEUse/conddata/yrcomp"
dirExists(resdir)
setuphtml(resdir=resdir,cleanslate = TRUE)
starttime <- as.character(Sys.time())
for (blk in c("6","7","8","9","10","11","12","13")) { # blk="6"
  pick <- which((wzl$block == blk))
  yrs <- sort(unique(wzl$year[pick]))
  nyrs <- length(yrs)
  yrblk <- as.matrix(table(wzl$year[pick],wzl$month[pick]))
  columns <- c(colnames(yrblk),"total")
  yrblkt <- cbind(yrblk,rowSums(yrblk,na.rm=TRUE))
  colnames(yrblkt) <- columns
  tfilen <- paste0("year_blk_",blk,".csv")
  addtable(yrblkt,tfilen,resdir,category=blk,caption="Records by year x month")
  for (yr in 1:nyrs) {
    year <- yrs[yr]
    filen <- paste0("length_comp_",blk,"_",year,".png")
    plotyearcomp(wzl,blk,year,bins=seq(120,230,2),filen,resdir)
  }
}
endtime <- as.character(Sys.time())

reportlist <- list(starttime=starttime,endtime=endtime,runname="Length_Comp")

runnotes <- "This is the first conditioning on the west coast."
# If you unhash this component it will generate a local website inside
# resdir and open it so you can see the results so far.
make_html(replist=reportlist,resdir=resdir,width=500,
          openfile=TRUE,runnotes=runnotes,verbose=FALSE,
          packagename = "aMSE",htmlname="yearlencomp")




pick <- which((wzl$year == 1976)) # & (wzl$block == "8"))

counts <- as.matrix(table(wzl$length[pick]))
x <- cbind(as.numeric(rownames(counts)),counts)
plotprep(width=7,height=4,newdev=FALSE)
inthist(x,col="orange",border=4,width=1,inc=5)


#select data for OM ------------------------------------------------------------

pickY <- which((wzl$year > 1980) & (wzl$length < 227))

records <- as.matrix(table(wzl$year[pickY],wzl$block[pickY]))
blks <- sort(as.numeric(colnames(records)))
nblk <- length(blks)

early <- droplevels(wzl[pickY,c(2,8,9,13,14)])
head(early,20)

dim(early)

early[,"block"] <- as.numeric(early[,"block"])

early <- early[order(early[,"block"]),]
early$len2 <- NA
early$len2 <- round(early$length/2)*2

ans <- table(early$len2,early$year,early$block)
str(ans)
lens <- seq(120,224,2)
lengths <- rep(lens,8)
blocknam <- c(6,7,8,9,10,11,12,13)
blocks <- c(rep(6,53),rep(7,53),rep(8,53),rep(9,53),
            rep(10,53),rep(11,53),rep(12,53),rep(13,53))
years <- as.numeric(dimnames(ans)[[2]])

dimnames(ans)
result <- matrix(0,nrow=424,ncol=34,dimnames=list(lengths,c("block",years)))
result[,"block"] <- blocks
begin <- 1
for (i in 1:8) {
  result[begin:(begin+52),2:34] <- ans[,,i]
  begin <- begin+53
}

filename <- "C:/Users/User/Dropbox/A_Code/aMSEUse/conddata/lateLF-84-20.csv"
write.csv(result,file=filename)

lateLF <- read.csv(file=filename,header=TRUE)




# extent of weight measurements ------------------------------------------------

pickw <- which(wzl$weight > 0)
tmp <- droplevels(wzl[pickw,])
table(tmp$block,tmp$year)







pick <- which((wzl$weight > 400) & (wzl$weight < 3000) & (wzl$length > 135))# &
            #  (wzl$year == 2020) & (wzl$block == "11"))
length(pick)
b12 <- droplevels(wzl[pick,])
# pickd <- which((b12$length > 180) & (b12$weight < 800))
# b12 <- b12[-pickd,]

# pick <- which(wzl$weight > 0)
# b12 <- droplevels(wzl[pick,])


model <- lm(log(b12$weight) ~ log(b12$length), data=b12)
summary(model)
coef(model)
model2 <- lm(weight ~ length, data=b12)
summary(model2)
x <- 100:210
y <- exp(coef(model))[1] * (x ^ coef(model)[2])
table(b12$year)
y2 <- coef(model2)[1] + (x * coef(model2)[2])

plotprep(width=7,height=5,newdev=FALSE)
plot1((b12$length),(b12$weight),type="p",pch=1,xlab="Shell Length mm",
      ylab="Total Weight",col=(b12$year - 2018),defpar=FALSE,xlim=c(100,210))
lines(x,y,lwd=2,col=4)
#lines(x,y2,lwd=2,col=3)


plotprep(width=5,height=8,newdev=FALSE)
parset(plots=c(4,1),margin=c(0.3,0.3,0.05,0.05),outmargin = c(1,1,0.05,0.05),
       byrow=FALSE)
plot(model2)

# which records have weights?

pick <- which(wz$weight > 0)
length(pick)

table(wz[pick,"block"])



pick <- which((wz$block == "13") & (wz$weight > 0))# &
                # (wz$whole.weight < 2000) & (wz$shell.length > 136) &
                # (wz$whole.weight > 380) & (wz$fishyear == 2020))
length(pick)
table(wz[pick,"date"])





plotprep(width=7,height=5,newdev=FALSE)
plot1(wz[pick,"length"],wz[pick,"weight"],type="p",pch=1,
      xlab="Shell Length mm",
      ylab="Total Weight",col=1,defpar=FALSE)


# auto-documentation ---------------------------------------------------------

listFunctions <- function(infile,verbose=TRUE) { # infile=filepath; verbose=TRUE
  content <- readLines(con=infile)
  funLines <- grep("function",content)
  nLine <- length(funLines)
  delF <- NULL
  for (i in 1:nLine) {
    tmpLine <- gsub(" ","",content[funLines[i]])
    if ((length(grep("function\\(",tmpLine)) == 0) |
        (substr(tmpLine,1,2) == "#'") |
        (length(grep("<-function",tmpLine)) == 0) |
        (length(grep("} #",tmpLine)) > 0)) delF <- c(delF,i)
  }
  ndelF <- length(delF)
  if (ndelF > 0) funLines <- funLines[-delF]
  if (ndelF == nLine) {
    txt <- paste0(infile,"  contained no recognizable functions")
    warning(cat(txt,"\n"))
    out=txt
    funnames=NULL
    funLines=NULL
  } else {
    outlines <- sort(c(funLines))
    out <- content[outlines]
    funnames <- out
    n <- length(out)
    for (i in 1:n) {  # i=1
      out[i] <- unlist(strsplit(out[i],"#"))[1]
      funnames[i] <- removeEmpty(unlist(strsplit(out[i],"<-"))[1])
    }
    if (verbose) {
      cat("\n",infile,"contained: \n")
      for (i in 1:n) {
        if (nchar(out[i]) == 2) {
          cat("\n")
        } else {
          cat(outlines[i],out[i],"\n")
        }
      }
      cat("\n")
    }
  }
  return(invisible(list(functions=out,funnames=funnames,funlines=funLines)))
} # end of listFunctions




fundocs <- function(package=".",verbose=TRUE) {
#  package="."; verbose=FALSE
  rfiles <- dir(paste0(package,"/R/"))
  rinfo <- file.info(paste0(package,"/R/",rfiles))
  rinfo <- rinfo[,-c(2,3,7)]
  nrfile <- length(rfiles)
  contents <- vector(mode="list",nrfile)
  names(contents) <- rfiles
  funnames <- NULL
  for (i in 1:nrfile) { # i = 1
    filepath <- paste0(package,"/R/",rfiles[i])
    out1 <- listFunctions(filepath,verbose=verbose)
    if (is.null(out1$funnames)) {
      detail <- out1$functions
    } else {
      nfun <- length(out1$funnames)
      columns <- c("name","linen","syntax")
      detail <- as.data.frame(matrix("",nrow=nfun,ncol=length(columns),
                           dimnames=list(out1$funnames,columns)))
      funnames <- c(funnames,out1$funnames)
      detail$name <- out1$funnames
      detail$linen <- out1$funlines
      ans <- strsplit(out1$functions,"<-")
      syn <- character(nfun)
      for (j in 1:nfun) syn[j] <- ans[[j]][2]
      syn <- gsub(" ","",syn)
      detail$syntax <- substr(syn,9,(nchar(syn)-1))
    }
    contents[[i]] <- detail
  }
  return(list(contents=contents,funnames=sort(funnames)))
} # end of fundocs

out <- fundocs(package="C:/Users/User/Dropbox/A_Code/aMSE",verbose=TRUE)
str(out,max.level = 1)




str(contents)


out <- tryCatch((is.function(scaletoone)),
                error=function(cond) print("not a function"))


