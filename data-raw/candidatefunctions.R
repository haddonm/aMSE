

# prepare size-composition data ------------------------------------------------

B0 <- getvar(zoneC,"B0")
ExB0 <- getvar(zoneC,"ExB0")
zoneDsau <- zonetosau(zoneDDR,NAS,glb,B0,ExB0)
zonePsau <- zonetosau(zoneDP,NAS,glb,B0,ExB0)


plotprep(width=8,height=8,newdev=FALSE)
x <- plotNt(zonePsau$Nt,year=30,glb=glb,start=3,medcol=1)
str(x)


lcomp <- prepareDDNt(zoneDD$Nt,zoneDD$catchN,glb)

plotCNt(lcomp$Nt,glb,vline=c(140),start=3) # plot the conditioning history

plotCNt(zonePsau$Nt[,,,51],glb,vline=140,start=3) # plot a single replicate from projections


zonePsau <- zonetosau(zoneDP,NAS,glb,B0,ExB0)

plotprep(width=6,height=8,newdev=FALSE)

parset(plots=c(4,2))
saunames <- glb$saunames
nsau <- glb$nSAU
cpue <- zonePsau$cpue
exploitb <- zonePsau$exploitB
for (sau in 1:nsau) {
  ymax <- getmax(cpue[,sau,1])
  xmax <- getmax(exploitb[,sau,1])
  plot(exploitb[,sau,1],cpue[,sau,1],type="p",pch=16,cex=1.0,panel.first=grid(),
       )
}


sau <- 1
title <- paste0("SAU ",saunames[sau])
pickyr <- c(1,5,10,15,20,25,30)
start <- 3
Nt <- lcomp$Nt[,,sau]
Nclass <- glb$Nclass
midpts <- glb$midpts
xmax <- getmax(Nt[start:Nclass,])/1000.0
nplot <- length(pickyr)
parset(plots=c(1,nplot),margin=c(0.3,0.1,0.05,0.05),outmargin=c(1,2,0,0),bty="n")
for (i in 1:nplot) {
   plot(Nt[,pickyr[i]]/1000.0,midpts,type="l",lwd=2,col=1,xlim=c(0,xmax),
        ylim=rev(range(midpts)),panel.first=grid(),ylab="",xlab="",xaxs="i")
  text(0.3*xmax,(Nclass * 2) + 5,paste0("yr ",pickyr[i]),cex=1.1,pos=4)
}
mtext("Shell Length mm",side=2,line=0.75,outer=TRUE,cex=1.1)
mtext(paste0(title,"  Numbers-at-Size 000's"),side=1,line=-0.2,outer=TRUE,cex=1.1)



# plot all first years and all last years
start <- 3
Nt <- zonePsau$Nt
Nclass <- glb$Nclass
midpts <- glb$midpts
reps <- dim(Nt)[4]
sc <- start:Nclass
sizes <- midpts[sc]
dat <- Nt[sc,,,]
nsau <- glb$nSAU
saunames <- glb$saunames
endyr <- dim(Nt)[2]
plotprep(width=8,height=8,newdev=FALSE)
parset(plots=c(4,2))
for (sau in 1:nsau) { #  sau=1
  ymax <- getmax(dat[,c(1,endyr),sau,],mult=1.01)
  plot(sizes,dat[,1,sau,1],type="l",lwd=1,col=rgb(.211,.211,.211,1/40),
       panel.first=grid(),
       ylim=c(0,ymax),ylab=saunames[sau],xlab="Shell Length mm")
  for (i in 2:reps) lines(sizes,dat[,1,sau,i],lwd=1,col=rgb(.231,.251,.251,1/40))
  for (i in 1:reps) lines(sizes,dat[,5,sau,i],lwd=1,col=rgb(1,0,0,1/50))
}


# Read the Obs LF-comp data-----------------------------------------------

options("show.signif.stars"=FALSE,
        "stringsAsFactors"=FALSE,
        "max.print"=50000,
        "width"=240)
# declare libraries ------------------------------------------------------------
library(aMSE)
library(rutilsMH)
# Obviously you should modify the rundir and datadir to suit your own setup
if (dir.exists("c:/Users/User/DropBox")) {
  ddir <- "c:/Users/User/DropBox/A_codeUse/aMSEUse/scenarios/"
} else {
  ddir <- "c:/Users/Malcolm/DropBox/A_codeUse/aMSEUse/scenarios/"
}
doproject <- FALSE  # change to FALSE if only conditioning is required
verbose <- TRUE
rundir <- paste0(ddir,"HS652510")
datadir <- paste0(ddir,"tasdata")

zone1 <- readctrlfile(rundir,infile="controlsau.csv",datadir=datadir,verbose=verbose)

compdat <- zone1$condC$compdat


palfs=compdat$palfs
palfs


# copyto -------------------------------------------------------------


rundir <- "c:/Users/User/DropBox/A_codeUse/aMSEUse/scenarios/HS652510"

destdir <- "c:/Users/User/DropBox/A_codeUse/aMSEUse/scenarios/HS81"

copyto(rundir,todir=destdir,filename="controlsau.csv")

  #     zoned=zoneDsau;zonep=zonePsau;glb=glb;startyr=30;
  #     picksau=9; histCE=histCE;CIprobs=c(0.05,0.5,0.95); addCI=TRUE



# selectivity plots -----------------------------------------------------

rundir=rundir; glb=out$glb;select=out$zoneC[[1]]$Select; histyr=out$condC$histyr;projLML=out$projC$projLML
rge=50:90









select <- out$zoneC[[1]]$Select

rge <- 50:105
plotprep(width=7,height=4,newdev=FALSE)


# sort out LF data-------------------------------------------------

datafile <- paste0(ddir,"aMSEUse/condition/productivity/sau10LF.csv")

lf <- read.csv(file=datafile,header=TRUE)
str1(lf)


yrcount <- tapply(lf[,"counts"],lf[,"year"],sum,na.rm=TRUE)
lfprops <- cbind(years,yrcount)
lfprops

range(lf[,"length"])
mids <- seq(136,210,2)
nc <- length(mids)
nc




answer <- makewidedat(lf,mids)
str(answer)

answer <- makewidedat(lf10,mids)



# Auto-regressive recruitment ------------------------------------------------


# R = exp(Xt)   where X is the predicted mean from the Beverton-Holt equation

oneyearrec <- function(steep,R0,B0,Bsp,sigR,devR=-1,rho=0.0) {
  if (devR > 0) {
    epsilon <- devR
  } else {
    epsilon <- exp(rnorm(length(Bsp),mean=0,sd=sigR) - (sigR * sigR)/2)
  }
  rec <- ((4*steep*R0*Bsp)/((1-steep)*B0+(5*steep-1)*Bsp)) * epsilon
  return(rec)
} # end of oneyearrec







# pot Catch vs Productivity function -------------------------------------------------------

popdefs <- getlistvar(out$zoneC,"popdef")
propD <- getzoneprops(out$zoneC,out$zoneDD,out$glb,year=glb$hyrs)
msy <- findmsy(out$zone$product)[,"Catch"]
nsau <- glb$nSAU
saunames <- glb$saunames
plts <- getparplots(length(msy))
histcat <- out$condC$histCatch
yrs <- as.numeric(rownames(histcat))
nyrs <- length(yrs)
plotprep(width=6,height=8,newdev = FALSE)
parset(plots=plts)
for (i in 1:nsau) {
  plot(yrs,histcat[,i],type="l",lwd=2,xlab="",ylab=paste0("SAU ",saunames[i]),
       panel.first=grid())
  abline(h=c(msy[i],geomean(histcat[,i])),lwd=c(2,1),col=c(2,4),
         lty=c(1,2))
 label <- paste0("depletion = ",round(propD["SpBDepl",i],3))
 x <- trunc(nyrs/4)
 text(yrs[x],5,label,cex=1.0,pos=4)
}






x <- matrix(rnorm(25,mean=5,sd=1),nrow=5,ncol=5)
kablerow(x,rowdigits=c(2,3,4,3,2))

rownames(x) <-  c("a","b","c","d","e")
colnames(x) <- c(1:5)
kablerow(x,rowdigits=c(2,3,4,3,2),namerows=TRUE)


x <- matrix(rnorm(25,mean=5,sd=1),nrow=5,ncol=5)
numdig <- c(2,3,4,3,2)
colnames(x) <- 1:5
kable(x,digits=numdig)
rownames(x) <- c("a","b","c","d","e")
kablerow(x,rowdigits=c(2,3,4,3,2),namerows=TRUE)



# compare runs-----------------------------------------------------------------

outzone6 <- out$outzone
TAC6 <- t(outzone6$TAC)
med2016 <- apply(TAC6,2,median)

load(paste0(ddir,"HS652510/ozoneDP.RData"))

paste0(ddir,"aMSEUse/scenarios/HS652510_2019/")
outzone9 <- outzone
TAC9 <- t(outzone9$TAC)
med2019 <- apply(TAC9,2,median)

yrs <- 2017:2049
plotprep(width=7, height=4,newdev = FALSE)
parset()
ymax <- getmax(c(med2016,med2019))
plot(yrs,c(med2016,NA,NA,NA),type="l",lwd=2,xlab="",ylim=c(0,ymax),
     ylab="Median TAC from 100 Replicates",panel.first=grid())
lines(yrs,c(NA,NA,NA,med2019),lwd=2,col=2)


load(paste0(rundir,"/out.RData"))

str1(out)

condC <- out$condC
saunames <- paste0("sau",6:13)
hcatch <- condC$histCatch
hcatch <- hcatch[-1,]
totC <- rowSums(hcatch)


plotprep(width=8,height=6,newdev=FALSE)
plotzonesau(zonetot=totC,saudat=hcatch,saunames=saunames,label="Zone Catch (t)",
            labelsau="sauCatch",side=4,sauscale=FALSE)



hce <- condC$histCE
yrs <- as.numeric(rownames(hcatch))
ceyrs <- as.numeric(rownames(hce))
pick <- match(ceyrs,yrs)
relC <- hcatch[pick,]
zonece <- catchweightCE(cedat=hce,cdat=relC,nsau=out$glb$nSAU)

plotprep(width=8,height=6,newdev=FALSE)
plotzonesau(zonetot=zonece,saudat=relC,saunames=saunames,label="Zone CPUE",
            labelsau="sauCatch",side=3,sauscale=FALSE)


plotprep(width=8,height=6,newdev=FALSE)
plotzonesau(zonetot=zonece,saudat=hce,saunames=saunames,label="Zone CPUE",
            labelsau="sauCPUE",side=1,sauscale=FALSE)



zoneDP <- out$zoneDP
str1(zoneDP)
NAS <- list(Nt=zoneDP$Nt,NumNe=zoneDP$NumNe,catchN=zoneDP$catchN)



str1(zoneDyn)

# change controlsau.csv by linenumber-------------------------------------------


findlinenumber(rundir,"controlsau.csv")

changeline(rundir,"controlsau.csv",210,"1991,-1,-1,-1,-1,-1,-1,-1,-1,")

changeline(rundir,"controlsau.csv",210,"1991,2,2,2,2,2,2,2,2,")


# SAU assessment ---------------------------------------------------------------

options("show.signif.stars"=FALSE,
        "stringsAsFactors"=FALSE,
        "max.print"=50000,
        "width"=240)
# declare libraries
library(aMSE)
library(rutilsMH)
# Obviously you should modify the rundir and datadir to suit your own setup
if (dir.exists("c:/Users/User/DropBox")) {
  ddir <- "c:/Users/User/DropBox/A_codeUse/aMSEUse/scenarios/"
} else {
  ddir <- "c:/Users/Malcolm/DropBox/A_codeUse/aMSEUse/scenarios/"
}
doproject <- FALSE  # change to FALSE if only conditioning is required
verbose <- TRUE
rundir <- paste0(ddir,"HS652510")
datadir <- paste0(ddir,"tasdata")

zone1 <- readctrlfile(rundir,infile="controlsau.csv",datadir=datadir,verbose=verbose)
consts <- getsaudata(datadir,zone1$ctrl$datafile)
condC <- zone1$condC
str(condC)
sau <- 7 #  sau12
LF <- condC$compdat$lfs[,,sau]
catch <- condC$histCatch[,sau]
ce <- condC$histCE[,sau]
LML <- condC$histyr
biol <- consts[,sau]

# identify years with data
pickc <- which(catch>0)
catch <- catch[pickc]
LML <- LML[pickc,]
yrs <- as.numeric(names(catch))
yrce <- as.numeric(names(ce))
pickce <- match(yrce,yrs)
cena <- putNA(ce,(pickce[1]-1),0)
names(cena) <- yrs
cbind(catch,cena,LML[,2])

printV(round(biol,5))





# align projections -------------------------------------------------------------------

histce <- out$sauout$zoneDsau$cpue[58,1,]
histcer <- out$zoneDDR$cesau[58,1,]
projce <- out$sauout$zonePsau$cpue[1,1,]

cat(mean(histce),mean(histcer),mean(projce))

plotprep(width=8,height=4,newdev=FALSE)
parset(plots=c(2,1))
hist(histce,breaks=25)
hist(histcer,breaks=25)



zoneC <- out$zoneC
condC <- out$condC


plotselectivity


saudyn <- poptosauCE(zoneDD$catch[1:finalyr,],zoneDD$cpue[1:finalyr,],sauindex)

for (i in 1:finalyr) { # i = 10
  saudyn <- poptosauCE(zoneDD$catch[i,],zoneDD$cpue[i,],sauindex)
  print(saudyn)
}


















