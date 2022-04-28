




rundir

pop <- read.csv(paste0(rundir,"/zonebiology.csv"),header=TRUE)

plotprep(width=8, height=8,newdev=FALSE)
parset()
pairs(pop[,c("M","MaxDL","L50","L95","AvRec","steep")])



grow <- read.csv("C:/Users/Malcolm/Dropbox/A_CodeUse/aMSEUse/condition/growth/hel_etal_2011.csv")


plotprep(width=8, height=8,newdev=FALSE)
parset()
pairs(grow[,c("MaxDL","L50","L95","isd")])


plot1(grow[,"L50"],grow[,"MaxDL"],type="p",pch=16,cex=1)

model1 <- lm(grow[,"L50"] ~ grow[,"MaxDL"])
summary(model1)
anova(model1)


model <- lm(grow[,"L95"] ~ grow[,"L50"] + grow[,"MaxDL"])
summary(model)
anova(model)

dat2 <- grow[,c("L50","MaxDL")]

y1 <- predict(model,newdata=dat2)
y2 <- predict(model,newdata=dat2,se.fit=TRUE)

y2f <- rnorm(length(y2$fit),mean=y2$fit,sd=y2$se.fit)


plot1(grow[,"L50"],grow[,"L95"],type="p",pch=16,cex=1)
points(grow[,"L50"],y2f,pch=16,col=2,cex=1)
points(grow[,"L50"],y1,pch=16,col=3,cex=1)




library(mvtnorm)

dat <- grow[,c("L50","L95")]
vcov <- cov(dat)
summary(dat)

x <- rmvnorm(100,mean=c(113.1,150.8),sigma=vcov)

# prepare size-composition data ------------------------------------------------

zoneCP <- out$zoneCP
zoneDP <- out$zoneDP
NAS <- out$NAS
glb <- out$glb
B0 <- getvar(zoneCP,"B0")
ExB0 <- getvar(zoneCP,"ExB0")
zoneDsau <- zonetosau(zoneDDR,NAS,glb,B0,ExB0)
zonePsau <- zonetosau(zoneDP,NAS,glb,B0,ExB0)


plotprep(width=8,height=8,newdev=FALSE)
x <- plotNt(zonePsau$Nt, year=58, glb=glb, start=3, medcol=1)
str(x)

numbersatsizeSAU(rundir="",out$glb,zoneC=out$zoneCP,zoneD=zoneDP, sau=4, ssc=5,
                 yr=58, defpar=TRUE, exploit=TRUE, mature=TRUE,filename="")


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







# plot Catch vs Productivity function -------------------------------------------------------

popdefs <- getlistvar(out$zoneC,"popdef")
propD <- getzoneprops(out$zoneC,out$zoneDD,out$glb,year=glb$hyrs)
msy <- findmsy(out$zone$product)[,"Catch"]
nsau <- glb$nSAU
saunames <- glb$saunames
plts <- pickbound(length(msy))
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




# HCR -------------------------------------------------------
str1(out)

glb <- out$glb
nsau <- glb$nSAU
hcr <- out$hcrout
cpue <- out$zoneDP$cesau
saucpue <- matrix(0,nrow=nrow(cpue),ncol=nsau)

for (i in 1:nsau) saucpue[,i] <- apply(cpue[,i,],1,median)

head(saucpue)
details <- hcr$details


rownames(details$scoret) <- 1992:2049
scoret <- details$scoret

invar=scoret
incpue <- saucpue[30:87,]

plotprep(width=8, height=10,newdev=FALSE)
parset(plots=c(8,1),margin=c(0.3,0.3,0.05,0.05))
yrs <- as.numeric(rownames(invar))
for (i in 1:nsau) {
  ymax <- getmax(hcr$multTAC[,i])
  plot(yrs,hcr$multTAC[,i],type="l",lwd=2,xlab="",ylab="",ylim=c(0.3,1.5),yaxs="i",
       panel.first=grid())
  # plot(yrs,scoret[,i],type="l",lwd=2,xlab="",ylab="",ylim=c(0,10),yaxs="i",
  #      panel.first=grid())
  # lines(yrs,details$score4[,i],lwd=2,col=2)
  # lines(yrs,details$score1[,i],lwd=2,col=3)
  # lines(yrs,details$scoretot[,i],lwd=4,col=4)
  # abline(h=5,lwd=1,col=2)
  abline(v=2020.5,lwd=1,col=2)

}



hcr$refpts



# Compare CatchN through time ----------------------------------------------------------------
glb <- out$glb
product <- out$production
filen=""
nsau <- glb$nSAU
npop <- glb$numpop
sauindex=glb$sauindex
saunames <- glb$saunames
catchN <- out$zoneDD$catchN[55:glb$Nclass,,]
mids <- glb$midpts[55:glb$Nclass]
saupops <- as.numeric(table(sauindex))
commonyaxs <- FALSE
count <- 0
plotprep(width=9,height=9,newdev=FALSE,filename=filen,verbose=FALSE)
parset(plots=c(8,7),margin=c(0.01,0.01,0.01,0.01),outmargin=c(3,3,0.1,0.1))
for (i in 1:npop) { # i=1
  if (commonyaxs) {
    ymax <- getmax(catchN[,58,])
  } else {
    ymax <- getmax(catchN[,58,i])
  }
  count <- count + 1
  label <- paste0(saunames[sauindex[i]],"_",count)
  if (sauindex[(i+1)] != sauindex[i]) count <- 0
  plot(mids,catchN[,2,i],type="l",lwd=2,xlab="",ylab="",xaxt="n",yaxt="n",
       yaxs="i",panel.first=grid(),ylim=c(0,ymax))
  lines(mids,catchN[,58,i],lwd=1,col=4)
  text(mids[1],0.9*ymax,label,cex=1.25,pos=4)
}

filen=""
nsau <- glb$nSAU
nyr <- 58
mids <- glb$midpts[55:glb$Nclass]
saunames <- glb$saunames
hyrnames <- glb$hyrnames

result <- array(0,dim=c(51,58,8),dimnames = list(mids,hyrnames,saunames))

for (i in 1:nyr) { # i = 1
  dat <- catchN[,i,]
  for (sau in 1:8) {
    pick <- which(sauindex == sau)
    result[,i,sau] <- rowSums(dat[,pick])
  }
}

years=56:58
plotprep(width=8,height=7,newdev=FALSE,filename=filen,verbose=FALSE)
parset(plots=pickbound(nsau),margin=c(0.3,0.3,0.05,0.05),outmargin=c(1,1,0,0))
for (sau in 1:nsau) { # i=1
  saucatch <- result[,years,sau]
  ymax <- getmax(saucatch)
  plot(mids,saucatch[,1],type="l",lwd=2,xlab="",ylab="",
       panel.first=grid(),ylim=c(0,ymax),yaxs="i")
  lines(mids,saucatch[,2],lwd=2,col=4)
  lines(mids,saucatch[,3],lwd=2,col=3)
  # text(0.7*max(saucpue[,i]),0.92*ymax,msylab,cex=1.2,pos=4)
  # text(0.8*max(saucpue[,i]),0.75*ymax,label[i],cex=1.5,pos=4)
  # text(1.1*msyce,0.15*ymax,round(msyce,2),cex=1.25,pos=4)
}







# SAU and Zone productivity ----------------------------------------------------

glb <- out$glb
nsau <- glb$nSAU
sauindex <- glb$sauindex
saunames <- glb$saunames
maxpop <- max(table(sauindex))
columns <- c(1:12,"SAU")
msybypop <- matrix(NA,nrow=nsau,ncol=(maxpop+1),dimnames=list(saunames,columns))
product <- out$production
popmsy <- findmsy(product)
saumsy <- tapply(popmsy[,"Catch"],out$glb$sauindex,sum)
zonemsy <- sum(saumsy)
for (i in 1:nsau) {
  picks <- which(sauindex == i)
  sauprod <- round(popmsy[picks,"Catch"],3)
  npop <- length(sauprod)
  msybypop[i,1:npop] <- sauprod
  msybypop[i,(maxpop+1)] <- sum(sauprod,na.rm=TRUE)
}
msybypop
popdefs <- t(getlistvar(out$zoneCP,"popdef"))

options(knitr.kable.NA = '')
kable(msybypop)


plotprep(width=8,height=4,newdev=FALSE,filen="")
parset(cex=1.0,margin=c(0.5,0.5,0.1,0.1))
plot(popmsy[,"Catch"],glb$SAUpop,col=sauindex,pch=16,cex=1.25,xlab="MSY by Population (t)",
     ylab="SAU",panel.first=grid())


plotprep(width=8,height=4,newdev=FALSE,filen="")
parset(cex=1.0,margin=c(0.5,0.5,0.1,0.1))
plot(popmsy[,"MatB"],popmsy[,"Catch"],col=sauindex,pch=16,cex=1.25,xlab="MSY by Population (t)",
     ylab="SAU",panel.first=grid())





model1 <- lm(popmsy[,"Catch"] ~ popmsy[,"MatB"])
summary(model1)

model2 <- lm(popmsy[,"Catch"] ~ popmsy[,"MatB"] + popdefs[,"steeph"])
summary(model2)

model3 <- lm(popmsy[,"Catch"] ~ popmsy[,"MatB"] + popdefs[,"steeph"] + popdefs[,"AvRec"])
summary(model3)

model4 <- lm(popmsy[,"Catch"] ~ popmsy[,"MatB"] + popdefs[,"steeph"] + popdefs[,"AvRec"] + popdefs[,"L50mat"])
summary(model4)


#anova(model1,model2)
anova(model2,model3)

AIC(model1,model2,model3,model4)
BIC(model1,model2,model3,model4)




microbenchmark(
  t(vect) %*% vect,
   vect %*% vect,
   crossprod(vect,vect),
   times=1000
)

x <- c(2,2,2,2,2)
y <- c(3,3,3,3,3)

x * y


# HS performance ---------------------------------------------------------

zoneDP <- out$zoneDP
catch <- zoneDP$catsau
glb <- out$glb
nsau <- glb$nSAU
sum10 <- getprojyrC(catsau=catch,glb=glb)





labelnames <- colnames(sum10)
plotprep(width=8, height=7)
parset(plots=c(3,3))
for (i in 1:(nsau+1)) {
  label <- paste0("Tonnes     ",labelnames[i])
  hist(sum5[,i],breaks=15,main="",xlab=label)
}


# Compare HS --------------------------------------------------------------
options("show.signif.stars"=FALSE,
        "stringsAsFactors"=FALSE,
        "max.print"=50000,
        "width"=240)
# declare libraries ------------------------------------------------------------
library(aMSE)
library(rutilsMH)
library(makehtml)
library(knitr)
# Obviously you should modify the rundir to suit your own setup
prefixdir <- "C:/A_Mal/scenarios/"

verbose <- TRUE
postfixdir <- "M15h75"
rundir <- paste0(prefixdir,postfixdir)
alldirExists(rundir,datadir,verbose=verbose)
source(paste0(rundir,"/TasmanianHS.R"))

controlfile <- "controlM15h75.csv"


zoneDP <- out$zoneDP
catch <- zoneDP$catsau
glb <- out$glb
nsau <- glb$nSAU


sum10 <- out$HSstats$sum10
sum5 <- out$HSstats$sum5

pd <- density(sum5[,(nsau+1)])
plotprep(width=8,height=4,newdev=FALSE)
parset(plots=c(1,8),margin=c(0.1,0.2,0.1,0.01),outmargin=c(3.5,1.75,0,0.5))
xmax <- getmax(pd$y)
plot(pd$y,pd$x,type="l",xlim=c(0,xmax),xaxs="i",ylab="")
mtext("M1h75_5",side=1,outer=FALSE,line=1.3,cex=0.9)
mtext("Relative Density",side=1,outer=TRUE,line=2.4)
mtext("10 Year Summmed Catch (t) by Zone",side=2,outer=TRUE,line=0.5)
pd2 <- density(sum10[,(nsau+1)])
xmax <- getmax(pd2$y)
plot(pd2$y,pd2$x,type="l",xlim=c(0,xmax),xaxs="i",ylab="")
mtext("M1h75_10",side=1,outer=FALSE,line=1.3,cex=0.9)



sauout <- out$sauout$zonePsau
object.size(sauout)
sauout <- sauout[-c(11,10)]
object.size(sauout)

str(HSstats)

catchN <- out$sauout$zonePsau$catchN
catchN <- catchN[56:105,,,]
nn <- object.size(catchN)
nn/118298272


load(paste0(rundir,"/sauoutD.RData"))
sauout <- sauoutD

ceCI <- out$sauout$outCI$cpue
glb <- out$glb
targCE <- hsargs$maxtarg

getmaxCE(ceCI=ceCI,glb=glb,targetCE=targCE)


projce <- out$sauout$zonePsau$cpue
glb <- out$glb
targetCE <- hsargs$maxtarg


reachtargCE <- function(projce,glb,targetCE) {
  pyrs <- glb$pyrs
  hyrs <- glb$hyrs
  totyrs <- hyrs + pyrs
  projcek <- projce[(hyrs+1):totyrs,,]
  nsau <- glb$nSAU
  label <- glb$saunames
  yrs <- glb$pyrnames
  reps <- glb$reps
  result <- matrix(0,nrow=reps,ncol=nsau,dimnames=list(1:reps,label))
  for (sau in 1:nsau) {
    for (i in 1:reps) {
      pick <- which(projcek[,sau,i] > targetCE[sau])
      result[i,sau] <- yrs[pick[1]]
    }
  }
}



sizecatchN <- function(catchN,glb,year,cutcatchN) {
  # cutcatchN=70; year=15; catchN <- out$sauout$zonePsau$catchN
  nsau <- glb$nSAU
  label <- glb$saunames
  reps <- glb$reps
  pyrs <- glb$pyrs
  hyrs <- glb$hyrs
  totyrs <- hyrs + pyrs
  yrs <- glb$pyrnames
  scl <- glb$midpts[cutcatchN:glb$Nclass]
  catchN <- catchN[(cutcatchN:glb$Nclass),(hyrs+1):totyrs,,]
  plotprep(width=7,height=6,newdev=FALSE,filen="")
  parset(plots=c(4,2),margin=c(0.3,0.45,0.05,0.05),outmargin=c(1,1,0,0))
  for (sau in 1:nsau) { # iter = 1
    catdat <- catchN[,year,sau,]
    ymax <- getmax(catdat,mult=1.02)
    plot(scl,catdat[,1],type="l",col="grey",ylim=c(0,ymax),panel.first=grid(),
         ylab=label[sau],xlab="")
    for (iter in 2:reps)
      lines(scl,catdat[,iter],col="grey")
    med <- apply(catdat,1,median)
    lines(scl,med,lwd=2,col=2)
  }
}


sausizecatchN <- function(catchN,glb,sau,years,cutcatchN) {
# cutcatchN=70;years=c(1,5,10,15,20,25,30);sau=6; catchN <- out$sauout$zonePsau$catchN
  nsau <- glb$nSAU
  label <- glb$saunames
  reps <- glb$reps
  pyrs <- glb$pyrs
  hyrs <- glb$hyrs
  totyrs <- hyrs + pyrs
  yrs <- glb$pyrnames
  scl <- glb$midpts[cutcatchN:glb$Nclass]
  catchN <- catchN[(cutcatchN:glb$Nclass),(hyrs+1):totyrs,,]
  plotprep(width=7,height=6,newdev=FALSE,filen="")
  parset(plots=c(4,2),margin=c(0.3,0.45,0.1,0.05),outmargin=c(1,1,0,0),
         byrow=FALSE)
  for (yr in years) { # iter = 1
    catdat <- catchN[,yr,sau,]
    ymax <- getmax(catdat,mult=1.02)
    plot(scl,catdat[,1],type="l",col="grey",ylim=c(0,ymax),panel.first=grid(),
         ylab=yrs[yr],xlab="")
    for (iter in 2:reps)
      lines(scl,catdat[,iter],col="grey")
    med <- apply(catdat,1,median)
    lines(scl,med,lwd=2,col=2)
  }

}





# postfixdir <- "SAfirst"
# rundir <- paste0(prefixdir,postfixdir)
# controlfile="controltest_SA.csv"
# hsargs=c(0,0,0,0,0,0)
# hcrfun=consthcr
# sampleCE=constCPUE
# sampleFIS=constFIS
# sampleNaS=constNaS
# getdata=constdata
# calcpopC=calcexpectpopC
# varyrs=7
# startyr=48
# verbose=TRUE
# ndiagprojs=4
# savesauout=TRUE
# makehcrout=makeouthcr
# cutcatchN=56


dat <- out$condout$sauZone$recruit[,1]
glb <- out$glb
loessfit <- loess(dat ~ as.numeric(glb$hyrnames),span=0.625)

plotprep(newdev=FALSE)
parset()
ymax <- getmax(dat)
plot(glb$hyrnames,dat,type="l",ylim=c(0,ymax),panel.first=grid())
lines(x$x,x$fitted,lwd=2,col=2)




# Get SizeComp Data -------------------------------------------------------------

#' @title numbersatsizeSAU plots the numbers-at-size for a given SAU
#'
#' @description numbersatsizeSAU plots the initial unfished numbers-
#'     at-size distribution for a given SAU, omitting the first four size
#'     classes to avoid the recruitment numbers dominating the plot.
#'
#' @param rundir the results directory, if set to "" then plot is sent to
#'     the console instead
#' @param glb the globals list
#' @param zoneC the constant part of the zone structure
#' @param Nt the numbers-at-size for all SAU, all years, and all reps
#' @param sau the SAU name to select for plotting. If multiple populations are
#'     contained in an SAU their numbers-at-size will be combined.
#' @param ssc index for starting size class. thus 1 = 2, 2 = 4, 5 = 10, etc.
#'     default = 5 for it plots size classes from 10mm up
#' @param yr which year of numbers-at-size should be plotted?
#' @param defpar should the plot parameters be defined. Set to FALSE if
#'     numbersatsizeSAU is to be used to add a plot to a multiple plot.
#' @param exploit should exploitable numbers be plotted as well? default=TRUE
#' @param mature should mature numbers be plotted as well? default=TRUE
#' @param filename default='Numbers-at-Size_Year1.png', but can be changed
#'     to suit whichever year is used
#'
#' @return nothing but it adds a plot to the results directory or console
#' @export
#'
#' @examples
#' print("this will be quite long when I get to it")
numbersatsizeSAU <- function(rundir, glb, zoneC, Nt, sau, ssc=5, yr=1,
                             defpar=TRUE, exploit=TRUE, mature=TRUE,
                             filename="Numbers-at-Size_Year1.png") {
  # rundir=""; glb=glb; Nt=Nt; sau="sau8"; defpar=TRUE; zoneC=zoneCP
  # `ssc=5; yr=1; exploit=TRUE; mature=TRUE; filename=""
  mids <- glb$midpts
  nc <- glb$Nclass
  saunames <- glb$saunames
  picksau <- which(saunames == sau)
  Nt <- as.matrix(Nt[,yr,picksau,]/1000.0)
  if (length(picksau) > 1) {
    Ntt <- as.matrix(rowSums(Nt,na.rm=TRUE))  # totals
  } else {
    Ntt <- Nt
  }
  if (nchar(rundir) > 0) {
    filen <- file.path(rundir,filename)
  } else {
    filen <- ""
  }
  if ((exploit) | (mature)) pickzC <- which(sapply(zoneC,"[[","SAU") == sau)
  maxy <- getmax(Ntt[ssc:nc,])
  if (defpar)
    plotprep(width=7,height=4,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  plot(mids[ssc:nc],Ntt[ssc:nc],type="l",lwd=2,xlab="Shell Length mm (5 - 210mm)",
       ylab="Numbers-at_size '000s",panel.first=grid(),ylim=c(0,maxy))
  if (exploit)
    lines(mids[ssc:nc],zoneC[[pickzC]]$Select[ssc:nc,yr] * Ntt[ssc:nc],col=4,lwd=2)
  if (mature)
    lines(mids[ssc:nc],zoneC[[pickzC]]$Maturity[ssc:nc] * Ntt[ssc:nc],col=2,lwd=2)
  legend("topright",legend=as.character(sau),lwd=0,col=0,bty="n",cex=1.2)
  if (nchar(rundir) > 0) {
    addm <- ""; adde <- ""
    if (mature) addm <- " The red line is mature numbers-at-size. "
    if (exploit) adde <- " The blue line is exploitable numbers-at-size."
    caption <- paste0("The numbers-at-size for the SAU ",sau," the recruitment ",
                      "numbers are omitted for clarity.",addm,adde)
    addplot(filen,rundir=rundir,category="NumSize",caption)
  }
} # end of numbersatsizeSAU





options("show.signif.stars"=FALSE,
        "stringsAsFactors"=FALSE,
        "max.print"=50000,
        "width"=240)
# declare libraries ------------------------------------------------------------
library(aMSE)
library(TasHS)
library(rutilsMH)
library(makehtml)
library(knitr)
# Obviously you should modify the rundir to suit your own setup
if (dir.exists("c:/Users/User/DropBox")) {
  prefixdir <- "c:/Users/User/DropBox/A_codeUse/aMSEUse/scenarios/"
} else {
  prefixdir <- "c:/Users/Malcolm/DropBox/A_codeUse/aMSEUse/scenarios/"
}
doproject <- TRUE  # change to FALSE if only conditioning is required
verbose <- TRUE
postfixdir <- "M15h75"
rundir <- paste0(prefixdir,postfixdir)
alldirExists(rundir,datadir,verbose=verbose)
controlfile <- "controlM15h75.csv"
hsfile <- "TasHS1_Tas.R"
source(paste0(rundir,"/",hsfile))

load(file=paste0("C:/aMSE_scenarios/",postfixdir,".RData"))

compdat <- out$condC$compdat

lfs <- compdat$lfs
palfs <- compdat$palfs
glb <- out$glb
zoneCP <- out$zoneCP
Nt <- out$NAS$Nt
sau <- "sau8"

dim(lfs)
label <- dimnames(lfs)
# sum predicted Nt for each sau for the historic period to cpmpare with data





hyrs=glb$hyrs
ssc <- 5
first=1
last=hyrs
Nt <- out$zoneDD$Nt

sauNt <- plotcondsizes(Nt,first,last,glb,ssc=5,filen="")

catchN <- out$zoneDD$catchN
ssc <- 65
first=30
last=hyrs-5

sauCN <- plotcondsizes(Nt=catchN,first,last,glb,ssc=ssc,filen="",
                       legloc="topright",prop=TRUE)


# compare conditoned predict CatchN ----------------------------------------



library(freqdata)


filen <- "c:/Users/Malcolm/DropBox/A_code/freqdata/data-raw/compiledMM.df.final.RDS"
comp <- readRDS(filen)
columns <- colnames(comp)
columns
newnames <- c("docket","msrdate","proc","procname","procnum","proclist","zone",
              "ndays","daylist","day_max","msrdate_diff","nblocks","blocks",
              "nsubblocks","subblocks","catch","length","n","meanSL","minSL",
              "species","weightg","datasource","blockno","subblockno","sampleid",
              "year","block1","block2","block3","block4","block5","id","region1",
              "region2","region3","region4","region5","same.region")
colnames(comp) <- newnames
props <- properties(comp)

blks <- c("6","7","8","9","10","11","12","13")
nblk <- length(blks)
pickb <- which((comp$blocks %in% blks) & (comp$year > 1989) &
                 (comp$year < 2021) & (comp$zone == "W") &
                 (comp$length >= 135) & (comp$length < 211))
cols <- c(7,12,13,16,17,18,21,22,26,27)

wz1 <- comp[pickb,cols]
nbby <- as.matrix(table(wz1$year,wz1$blocks))
nbby1 <- expandmatrix(nbby)
nbby2 <- cbind(nbby1[,5:8],nbby1[,1:4])

nbby2 - palfs[7:37,]

pick <- which((wz1$year == 2014) & (wz1$blocks == 12))
dat <- wz1[pick,]

plotprep(width=7,height=4,newdev=FALSE)
parset()
inthist(dat$length,width=0.8,border=2)

head(dat,20)

dat$len2 <- trunc((dat$length+1)/2)*2
inthist(dat$len2,width=0.8,border=2)


# Characterize APPs or Populations---------------------------------------------
# uses out$zoneCP and out$zoneDD

zoneCP <- out$zoneCP
zoneDD <- out$zoneDD
glb <- out$glb
hyrs <- out$glb$hyrs
propD <- as.data.frame(t(getzoneprops(zoneCP,zoneDD,glb,year=28)))

propD

str1(zoneCP[[1]])





columns <- c("B0","MSY","MSYDepl","bLML","propprot","SpBDepl","catch","harvestR")

plotpopprops(propD,rundir="",glb,columns,startyr=hyrs,bins=21)


plot1(propD[,"MSY"],propD[,"harvestR"],type="p",pch=16,cex=1.25)

parset()
pairs(propD[1:glb$numpop,c(1,5,7,9,11,14,15,16)],pch=16,col=2,cex=1)









makerect <- function (left, xinc, top, yinc, linecol = "grey",col = NULL)
{
  polygon(makevx(left, xinc), makevy(top, yinc), col = col,
          border = linecol)
  centerx <- (left * 2 + xinc)/2
  centery <- (top * 2 - yinc)/2
  return(invisible(c(centerx, centery)))
}
