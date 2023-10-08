
Nt=sauout$Nt; catchN=sauout$catchN;padding=5;prop=TRUE;minSL=90;sau=1;
 console=TRUE;omit=5;legal=TRUE


# Sizecomp in projections --------------------------------------







#min(outfish$yrlml)

predsizecomp(sau=5, NSC=sauout$Nt, glb=glb, minSL=70,interval=5,
             prop=FALSE,console=TRUE)









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






# compare runs-----------------------------------------------------------------


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
plot(popmsy[,"Catch"],glb$SAUnum,col=sauindex,pch=16,cex=1.25,xlab="MSY by Population (t)",
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



getsaucatchN <- function(catchN,glb,sau,year,ss,bysau=TRUE) {
#  catchN=catchN; glb=glb; sau=6; year=48; ss=1000
  reps <- glb$reps
  indexsau <- which(glb$sauindex == sau)
  numC <- catchN[,year,indexsau,]
  pickR <- trunc(runif(1,min=1,max=reps))
  nC <- rowSums(numC[,,pickR])
  nS <- round(ss*nC/sum(nC))
  actualss <- sum(nS)
  bootdat <- numeric(actualss)
  pickP <- which(nS > 0)
  datboot <- nS[pickP]
  sizes <- as.numeric(names(datboot))
  index <- 1
  for (i in 1:length(datboot)) { # i=1
    reps <- datboot[i]
    addat <- rep(sizes[i],reps)
    bootdat[index:(index+reps - 1)] <- addat
    index <- index+reps
  }
  ans <- sample(bootdat,size=actualss,replace=TRUE)
  return(ans)
} # end of getsaucatchN


samp <- getsaucatchN(catchN,glb,sau=6,year=47,ss=500)


sizecomp <- out$condC$compdat$lfs
lfs <- preparesizecomp(sizecomp[,,6],mincount = 120)
plfs <- prop.table(lfs,2)
sizes <- as.numeric(rownames(plfs))




plotprep(width=9, height=5,newdev=FALSE)
parset()
plot1(sizes,plfs[,13],lwd=2)
for (i in 1:4) {
  samp <- getsaucatchN(catchN,glb,sau=6,year=47,ss=500)
  sampt <- table(samp)
  ssizes <- as.numeric(names(sampt))
  counts <- as.numeric(sampt)
  props <- counts/sum(counts)
  lines(ssizes,props,lwd=2,col=(i+1))
}



bootdat <- numeric(tot[yr])
dat <- sizecomp[,yr]
pickC <- which(dat > 0)
index <- 1
for (i in 1:length(pickC)) { # i=1
  reps <- dat[pickC[i]]
  addat <- rep(sizes[pickC[i]],reps)
  bootdat[index:(index+reps - 1)] <- addat
  index <- index+reps
}






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


# plotphase

#' @title plotphase generates a kobe-plot of harvest rate against mature depletion
#'
#' @description plotphase generates a kobe-like plot of harvest rate again the
#'     mature biomass depletion. This can provide a summary of the stock status
#'     although currently, there is no limit refernece point defined,
#'
#' @param dyn a matrix of the predicted dynamics. It must include a column of
#'  h the dynamics are running.
#' @param answer
#' @param maxdepl this can be used to constrain the plot so that details can be
#'     seen which are obscured if depletionsup to 1.0 are used. default=1.0
#' @param startyr which year to start plotting the dynamics. default = 1992,
#'     which suits Tasmania but elsewhere will likely want somwething different.
#'
#' @return nothing, but is does plot a graph
#' @export
#'
#' @examples
#' print("wait on data sets")
plotphase <- function(dyn,answer,maxdepl=1.0,startyr=1992) {
  # dyn <- out$sauout; answer=outprod$answer; maxdepl=0.4; startyr=1995
  yrs <- as.numeric(rownames(dyn))
  picky <- which(yrs >= startyr)
  nyrs <- length(picky)
  matB <- dyn[picky,"matB"]
  depl <- dyn[picky,"depl"]
  H <- dyn[picky,"harvest"]
  rown <- nrow(answer)
  targH <- answer[rown,"HarvestR"]
  targdepl <- answer[rown,"Depletion"]
  plotprep(width=7, height=6, newdev=FALSE)
  parset(cex.lab=1.25)
  ymax <- getmax(H)
  plot(depl,H,type="p",pch=16,xlim=c(0,maxdepl),ylim=c(0,ymax),xaxs="i",yaxs="i",
       xlab="Mature Biomass Depletion",ylab="Harvest Rate")
  polygon(x=c(0,targdepl,targdepl,0,0),y=c(targH,targH,ymax,ymax,targH),
          col=rgb(255,0,0,175,maxColorValue=255))
  polygon(x=c(targdepl,maxdepl,maxdepl,targdepl,targdepl),
          y=c(targH,targH,ymax,ymax,targH),
          col=rgb(255,255,0,170,maxColorValue=255))
  polygon(x=c(0,targdepl,targdepl,0,0),y=c(0,0,targH,targH,0),
          col=rgb(255,255,0,175,maxColorValue=255))
  polygon(x=c(targdepl,maxdepl,maxdepl,targdepl,targdepl),y=c(0,0,targH,targH,0),
          col=rgb(0,255,0,120,maxColorValue=255))
  abline(h=targH,lwd=2,col=1)
  abline(v=targdepl,lwd=2,col=1)
  arrows(x0=depl[1:(nyrs-1)],y0=H[1:(nyrs-1)],x1=depl[2:nyrs],y1=H[2:nyrs],lwd=2,
         length=0.125)
  points(depl,H,pch=16,cex=1.25)
  points(depl[c(1,nyrs)],H[c(1,nyrs)],pch=16,cex=2,col=c(1,4))
} # end of plotphase


library(codeutils)
library(hplot)

#' # steep=steep;R0=R0;B0=B0;Bsp=820;sigR=insigmar;devR=-1;depensate=0.2
oneyearrec <- function(steep,R0,B0,Bsp,sigR,devR=-1,depensate=0) { #
  if (devR[1] > 0) {
    epsilon <- devR
  } else {
    epsilon <- exp(rnorm(length(Bsp),mean=0,sd=sigR) - (sigR * sigR)/2)
  }
  if ((depensate > 0) & (Bsp/B0 < depensate)) {
    bcurr <- B0 * depensate
    thres <- ((4*steep*R0*bcurr)/((1-steep)*B0+(5*steep-1)*bcurr))
    rec <- ((Bsp/B0)/depensate) * thres * epsilon
  } else {
    rec <- ((4*steep*R0*Bsp)/((1-steep)*B0+(5*steep-1)*Bsp)) * epsilon
  }
  return(rec)
} # end of oneyearrec

insigmar <- 1e-07
steep <- 0.7
R0 <- 1720637
B0 <- 4152
depen <- 0.2
Bsp <- seq(20,4000,20)
nB <- length(Bsp)
recs <- numeric(nB); names(recs) <- Bsp
for (i in 1:nB)
  recs[i] <- oneyearrec2(steep,R0,B0,Bsp[i],insigmar,depensate=depen)



plotprep(width=8, height=5,newdev=FALSE)
parset()
plot1(Bsp/B0,recs,defpar = FALSE,xlab="depletion",ylab="recruitment")
rec1=recs
for (i in 1:nB)
  rec1[i] <- oneyearrec2(steep,R0,B0,Bsp[i],insigmar,depensate=0)
lines(Bsp/B0,rec1,col=2)


# compdatatemplate -----------------------------------------------------------

library(aMSE)
library(codeutils)
library(makehtml)




rundir <- "C:/Users/Malcolm/Dropbox/A_CodeUse/aMSEUse/scenarios/example/"

compdatatemplate(rundir)



sau=1; len=1

cat(lens[len],",",saunames[sau],",",paste0(lfs[len,,sau],collapse=",")," \n")



plotprep(width=8,height=6,newdev=FALSE,filename=filen,verbose=FALSE)
parset(plots=c(1,1),margin=c(0.3,0.4,0.05,0.05),outmargin=c(0,1,0,0),
       byrow=FALSE)
maxy <- getmax(varq)
nscen <- length(scenes)
plot(yrnames,varq[3,,1],type="l",lwd=2,col=0,xlab="",ylab=varname,
     panel.first=grid(),ylim=c(0,maxy))

polygon(poldat4,col=RGB(8,alpha=255))


colname <- c("black","red","green3","blue","cyan","magenta","yellow","gray")








rundir=rundir
scenes=scenes
zone=zone
glb=glbc[[1]]
console=TRUE
q90=TRUE


plotzonedyn(rundir,scenes,zone,glb,console=TRUE,q90=TRUE,polys=TRUE,intens=150)






x1 <- finalcondyeardepletion(rundir,ans[[1]]$sauout,ans[[1]]$glb,deplvar="eB",console=TRUE)
x1


# Indiv population Nt ----------------------------

out <- result$ans[[1]]

nas <- out$NAS
Nt <- nas$Nt

sauout <- out$sauout
glb <- out$glb
str1(glb)


start <- glb$hyrs+1
finish <- glb$hyrs+glb$pyrs
invar <- sauout$matureB[start:finish,1,]

med <- apply(invar,1,median)




med <- getmedbysau(sauout$matureB,glb)

apply(med,2,which.max)
apply(med,2,which.min)


med1 <- getmedbysau(result$ans[[1]]$sauout$matureB,glb)
apply(med1,2,which.max)
med2 <- getmedbysau(result$ans[[2]]$sauout$matureB,glb)
apply(med2,2,which.max)

sau11 <- Nt[,start:finish,28:36,]

rge <- 8:glb$Nclass

sc <- glb$midpts[rge]

plotprep(width=9, height=5)
parset()
plot1(sc,sau11[rge,1,1,1])
for (i in 2:15) lines(sc,sau11[rge,1,1,i],lwd=1,col=i)



round(sau11[,1,1,1:3])


# Using quantscen----------------------------

quantscen <- result$quantscen

label <- names(quantscen)
nvar <- length(label)
scenes <- names(quantscen[[1]])


# AAV functioning------------------------------------


# Run compare_scenarios
# extract the catches from dyn and select one set of projections

cs <- catch[["Base_Case"]][59:88,,]

outmed <- plotsceneproj(rundir="",inarr=cs,glb=glbc[[1]],scene="Base_Case",
                        filen="",label="",maxy=0,Q=90,hline=NA)


aav <- getprojyraavc(cs,glbc[[1]])


getaavA <- function(invect) { # invect=x
  nyr <- length(invect)
  totC <- sum(abs(invect),na.rm=T)
  aac <- sum(abs(invect[2:nyr] - invect[1:(nyr-1)]))
  aav <- 0.0
  if (totC > 0.0) aav <- aac/totC
  return(aav)
}

singlesauproj <- function(arr3d,sau,rn,glb,filen="",rundir="") {
  # arr3d=cs;sau=4;rn=4;glb=glbc[[1]];filen=""
  outaav <- numeric(rn)
  sauname <- glb$saunames[sau]
  reps <- glb$reps
  dat <- arr3d[,sau,]
  yrs <- as.numeric(rownames(dat))
  outrep <- matrix(0,nrow=length(yrs),ncol=rn,dimnames=list(yrs,1:rn))
  meds <- apply(dat,1,median)
  maxy <- getmax(dat)
  plotprep(width=9, height=5,newdev=FALSE,filename = filen,verbose=FALSE)
  parset()
  plot(yrs,(dat[,1]),type="l",lwd=1,col="grey",ylab="",xlab="",ylim=c(0,maxy))
  for (i in 2:reps) lines(yrs,dat[,i],lwd=1,col="grey")
  pickr <- trunc(runif(rn,min=1,max=(reps+1)))
  for (i in 1:rn) {
    outrep[,i] <- dat[,pickr[i]]
    lines(yrs,outrep[,i],lwd=2,col=i)
    outaav[i] <- getaavA(outrep[,i])
  }
  return(list(outaav=outaav,outrep=outrep,meds=meds))
} # end of singlesauproj

out <- singlesauproj(cs,sau=4,rn=10,glb=glbc[[1]])
out$outaav


for (i in 1:10) {
  print(getaav(out$outrep[,i]-out$meds))
}


sau=6
cs <- catch[[1]][59:88,,]
dat <- cs[,sau,]
yrs <- as.numeric(rownames(dat))
y <- apply(dat,1,getaav)
plot1(yrs,y,lwd=2,col=1)
for (scen in 2:3) {
  cs <- catch[[scen]][59:88,,]
  dat <- cs[,sau,]
  y <- apply(dat,1,getaav)
  lines(yrs,y,lwd=2,col=scen)
}



plot1(yrs,out$outrep[,4]-out$meds)


getmaxmedian <- function(obj,glb) {
  nscen <- length(obj)
  scenes <- names(obj)
  meds <- getmedbysau(obj[[1]],glb)
  yrs <- as.numeric(rownames(meds))
  maxmed <- matrix(0,nrow=nscen,ncol=nsau,
                   dimnames=list(scenes,glb$saunames))
  maxmed[1,] <- yrs[apply(meds,2,which.max)]
  for (scen in 2:nscen) {
    meds <- getmedbysau(obj[[scen]],glb)
    maxmed[scen,] <- yrs[apply(meds,2,which.max)]
  }
  return(maxmed)
} # end of getmaxmedian


getmaxmedian(catch,glbc[[1]])
getmaxmedian(cpue,glbc[[1]])

#aav,sum5,sum10 by scenario -----------
ans <- calccatchHSPM(catch,glbc,scenes=names(catch))

sapply(ans$aavc,function(x) apply(x,2,median))
sapply(ans$sum5,function(x) apply(x,2,median))
sapply(ans$sum10,function(x) apply(x,2,median))

# Expanded do_comparison-----------------------


rundir=rundir;postfixdir=postfixdir;outdir=outdir;files=files;pickfiles=c(5,6)
verbose=TRUE; intensity=100

files2 <- files[pickfiles]
nfile <- length(pickfiles)
label <- vector(mode="character",length=nfile)
for (i in 1:nfile) label[i] <- unlist(strsplit(files2[i],".",fixed=TRUE))[1]
ans <- makelist(label) # vector(mode="list",length=nfile)
dyn <- makelist(label) # vector(mode="list",length=nfile)
glbc <- makelist(label) # vector(mode="list",length=nfile)
prods <- makelist(label) # vector(mode="list",length=nfile)
scenes <- vector(mode="character",length=nfile)
scores <- makelist(label) # vector(mode="list",length=nfile)
zone <- makelist(label)
for (i in 1:nfile) { # i = 1
  filename <- paste0(outdir,files2[i])
  if (verbose)
    cat("Loading ",files2[i]," which may take time, be patient  \n")
  out <- NULL # so the function knows an 'out' exists
  load(filename)
  ans[[i]] <- out
  dyn[[i]] <- out$sauout
  glbc[[i]] <- out$glb
  prods[[i]] <- out$sauprod
  scenes[i] <- out$ctrl$runlabel
  scores[[i]] <- out$scores
  zone[[i]] <- out$outzone
}
nscenes <- length(scenes)
catch <- makelist(scenes)
cpue <- makelist(scenes)
for (i in 1:nscenes) {
  catch[[i]] <- dyn[[i]]$catch
  cpue[[i]] <- dyn[[i]]$cpue
}
if (verbose) cat("Now doing the comparisons  \n")
setuphtml(rundir=rundir)
quantscen <- comparedynamics(rundir=rundir,dyn,glbc,scenes)
tabulateproductivity(rundir,prods,scenes)

filename <- "compare_final_HSscores.png"
meds <- comparefinalscores(rundir,scores,scenes,legloc="bottomright",
                           filen=filename,category="Scores")
addplot(filen=filename,rundir=rundir,category="Scores",
        caption="The HS final scores for each sau.")

tabulatefinalHSscores(rundir,meds,scenes,category="Scores")
outcatchHSPM <- calccatchHSPM(catch,glbc,scenes,aavcyrs=10)
medHSPM <- outcatchHSPM$medians
label <- names(medHSPM)
for (i in 1:length(medHSPM)) {
  filename <- paste0("proj",label[i],".csv")
  addtable(medHSPM[[label[i]]],filen=filename,rundir=rundir,category="HSPM",
           caption=paste0("Median ",label[i]," values by sau and scenario."))
}
outtab <- catchHSPM(rundir,hspm=outcatchHSPM,glbc,scenes,
                    filen="compare_catches_boxplots.png",aavcyrs=10)
boxbysau(rundir,hspm=outcatchHSPM,glbc,scenes,compvar="aavc",
         filen="sau_aavc_boxplots.png",aavcyrs=10,maxval=0.4)
boxbysau(rundir,hspm=outcatchHSPM,glbc,scenes,compvar="sum5",
         filen="sau_sum5_boxplots.png",aavcyrs=10,maxval=0) #
boxbysau(rundir,hspm=outcatchHSPM,glbc,scenes,compvar="sum10",
         filen="sau_sum10_boxplots.png",aavcyrs=10,maxval=0) #
catchbpstats(rundir,outtab)
cpueHSPM(rundir,cpue,glbc,scenes=scenes,filen="sau_maxce_boxplots.png")
outcpue <- cpueboxbysau(rundir,cpue,glbc,scenes,filen="sau_maxceyr_box.png",
                        startyr=0,maxval=0)
cdivmsy <- catchvsMSY(catch,glbc,prods,scenes)
nscen <- length(scenes)
for (i in 1:nscen) {
  plotsceneproj(rundir,cdivmsy[[i]],glbc[[i]],scenes[i],
                filen=paste0("Catch_div_MSY_",scenes[i],".png"),
                label="Catch / MSY",hline=1,Q=90)
}
plotzonedyn(rundir,scenes,zone,glbc[[1]],console=FALSE,q90=TRUE,polys=TRUE,
            intens=intensity)
tabulatezoneprod(rundir,prods,scenes)
# ribbon plots by sau and dynamic variable
# cpue <- scenebyvar(dyn=out$dyn,byvar="cpue",glb=out$glbc[[1]])
glb <- glbc[[1]]
catch <- scenebyvar(dyn,byvar="catch",glb=glb,projonly = TRUE)
catqnts <- sauquantbyscene(catch,glb)
nsau <- glbc[[1]]$nSAU
for (sau in 1:nsau)
  sauribbon(rundir,scenes=scenes,sau=sau,varqnts=catqnts,
            glb=glb,varname="Catch",console=FALSE,
            q90=TRUE,intens=intensity,addleg="bottomright")
cpue <- scenebyvar(dyn,byvar="cpue",glb=glb,projonly = TRUE)
cpueqnts <- sauquantbyscene(cpue,glb)
for (sau in 1:nsau)
  sauribbon(rundir,scenes=scenes,sau=sau,varqnts=cpueqnts,
            glb=glb,varname="cpue",console=FALSE,
            q90=TRUE,intens=intensity,addleg="bottomright")
makecompareoutput(rundir=rundir,glbc,scenes,postfixdir,
                  filesused=files[pickfiles],openfile=TRUE,verbose=FALSE)



yrs <- 1990:2020
nyrs <- length(yrs)
deplsB <- 0.01 + 0.00015*yrs * 0.15 + rnorm(nyrs,mean=0.1,sd=0.01)
reps <- 100
# generate a matrix with years in the rows.
deplsB <- matrix(0,nrow=nyrs,ncol=reps,dimnames=list(yrs,1:reps))
for (i in 1:reps)
  deplsB[,i] <- seq(0.1,0.4,length=nyrs) + rnorm(nyrs,mean=0.1,sd=0.04)
plotprep(width=8, height=3.5)
plot1(yrs,deplsB[,1])


#' @title trajectoryplot generates a plot comparing replicate trajectories
#'
#' @description trajectoryplot fulfils a common task with is to compare groups
#'     of trajectories generated in simulations having multiple scenarios.\
#'     Thus, one might have generated the mature biomass depletion for two
#'     or more scenarios across a selection of years, each with a large number
#'     of replicates. Such data is expected to be contained in a list of lists
#'     with each of the top level lists being one of the scenarios. Each of
#'     those lists is expected to be made up of multiple components. Hence one
#'     might have a list called 'dyn' containing lists for two scenarios 'meta'
#'     and 'refp'. Both 'meta' and 'refp' constitute a list of 'matureB',
#'     'exploitB', 'catch', 'cpue', 'deplsB', etc. Where each of those is a
#'     three dimensional array of 1:nyr years, 1:nsau spatial unuts, and 1:N
#'     replicates.
#'
#' @param x a vector of years to be used as the x-axis
#' @param y A matrix of years by whatever variable is to be plotted.
#'
#' @return
#' @export
#'
#' @examples
#' yrs <- 1990:2020
#' nyrs <- length(yrs)
#' reps=100
#' deplsB <- matrix(0,nrow=nyrs,ncol=reps,dimnames=list(yrs,1:reps))
#' for (i in 1:reps)
#'   deplsB[,i] <- seq(0.1,0.4,length=nyrs) + rnorm(nyrs,mean=0.1,sd=0.04)
#' # x=
trajectoryplot <- function(x,y,varname="",fnt=7,cex=1.0,addmedian=TRUE,bounds=90,
                           distyr=NULL,distvar=NULL) {
  if ((length(distyr) > 0) & (length(distvar) > 0))
      warning(cat("Both distyr and distvar in trajectoryplot have values. \n"))
  if (length(distyr) > 0) layout(rbind(c(1,1,1),c(2,3,4)),heights=c(3,2))
  if (length(distvar) > 0) layout(rbind(c(1,1),c(2,3)),heights=c(3,2))
#  x=yrs;y=deplsB;fnt=7;cex=1.0;addmedian=TRUE;varname="depletion";bounds=90;distvar=c(0.3,0.45)
  maxy <- getmax(y)
  plot(x,y[,1],type="l",lwd=1,col="grey",ylim=c(0,maxy),font=fnt,cex=cex,
       ylab=varname,xlab="")
  for (i in 2:reps) lines(x,y[,i],lwd=1,col="grey")
  medy <- apply(y,1,quants)
  if (addmedian) lines(x,medy[3,],lwd=3,col=2)
  if ((addmedian) & bounds > 0) {
    if (bounds == 90) {
      lines(x,medy[2,],lwd=2,col=2)
      lines(x,medy[4,],lwd=2,col=2)
    } else {
      lines(x,medy[1,],lwd=2,col=2)
      lines(x,medy[5,],lwd=2,col=2)
    }
  }
  if (!is.null(distyr)) {
    abline(v=yrs[distyr],lwd=2,col=1)
    for (i in 1:length(distyr)) {
      label <- paste0(varname," in ",yrs[distyr[i]])
      hist(y[distyr[i],],main="",xlab=label,ylab="Proportion",freq=FALSE,
           col="grey")
      lines(density(y[distyr[i],]),lwd=2,col=4)
    }
  }
  if (!is.null(distvar)) {
    abline(h=distvar,lwd=2,col=1)
    for (i in 1:length(distvar)) {
      label <- paste0("First Year ",varname," at ",distvar[i])
      varyrs <- apply(y,2,function(x) which(x >= distvar[i])[1])
      inthist(yrs[varyrs],xlab=label,ylab="frequency")
    }
  }
} # end of trajectoryplot


plotprep(width=9,height=6)
trajectoryplot(yrs,deplsB,varname="depletion",fnt=2,addmedian=TRUE,bounds=90,
               distvar=c(0.3,0.45),cex=1.0)

dvar <- numeric(reps)
for (i in 1:reps) dvar[i] <- which(deplsB[,i] >= distvar[1])[1]

inthist(yrs[dvar])
i=1

yearCE=condC$yearCE





twophaseplots(dyn=sauout,glb=glb,outhcr=outhcr,sau=6,kobeRP=kobeRP,rundir="",
             startyr=condC$yearCE[1],console=TRUE,fnt=7)



HSphaseplot(dyn=sauout,glb=glb,sau=plotsau,
           rundir=rundir,
           startyr=condC$yearCE[1],
           console=TRUE,setpar=TRUE,targdepl=kobeRP[1],
           limdepl=kobeRP[2],limH=kobeRP[3])

tasphaseplot(proxyB=outhcr$targsc,proxyH=outhcr$g4s,glb=glb,sau=6,rundir="",
             console=TRUE,fnt=7,setpar=TRUE)


# plotcompdat --------------------------------------

options("show.signif.stars"=FALSE,
        "stringsAsFactors"=FALSE,
        "max.print"=50000,
        "width"=240)
suppressPackageStartupMessages({
  library(aMSE)
  library(TasHS)
  library(codeutils)
  library(hplot)
  library(makehtml)
  library(knitr)
})
dropdir <- getDBdir()
prefixdir <- pathtopath(dropdir,"/A_codeUse/aMSEUse/hsargs/")

postfixdir <- "BC"
verbose <- TRUE
rundir <- path.expand(filenametopath(prefixdir,postfixdir)) # #"mult05"
controlfile <- paste0("control",postfixdir,".csv")
outdir <- "C:/aMSE_scenarios/BC/"
confirmdir(rundir,ask=FALSE)
confirmdir(outdir,ask=FALSE)


# current Tas HS
hsargs <- list(mult=0.1, # expansion factor for cpue range when calc the targqnt
               wid = 4, # number of years in the grad4 PM
               targqnt = 0.55, # quantile defining the cpue target
               maxtarg = c(150,150,150,150,150,150,150,150), # max cpue Target
               pmwts = c(0.65,0.25,0.1),  # relative weights of PMs
               hcr = c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2), # hcr multipliers
               startCE = 2000,
               endCE = 2019,
               metRunder = 0,  # should the metarules be used. o =
               metRover = 0,   # use metarules
               decrement=1, # use fishery data up to the end of the time series
               pmwtSwitch = 0, # number of years after reaching the targCE to
               stablewts = c(0.4, 0.5, 0.1), # replace pmwts with stablewts
               hcrname="mcdahcr",     # the name of the HCR used
               printmat=NULL)

# refperiodCE = 2000:2019, # use a vector of years instead

checkhsargs(hsargs)

hcrfun=mcdahcr  #constantrefhcr,    #   #
sampleCE=tasCPUE;sampleFIS=tasFIS;sampleNaS=tasNaS
getdata=tasdata;calcpopC=calcexpectpopC;makeouthcr=makeouthcr
fleetdyn=NULL;scoreplot=plotfinalscores
plotmultflags=plotmultandflags;interimout=""
varyrs=7;startyr=38;verbose=TRUE;ndiagprojs=4;cutcatchN=56
matureL=c(40,170);wtatL=c(50,200);mincount=120
includeNAS = FALSE;depensate=0;kobeRP=c(0.4,0.2,0.15)
nasInterval=5;minsizecomp=c(100,135)  # A NEW LINE


  # generate equilibrium zone -----------------------------------------------
  starttime <- (Sys.time())
  zone <- makeequilzone(rundir,controlfile,verbose=verbose)
  equiltime <- (Sys.time()); if (verbose) print(equiltime - starttime)
  # declare main objects ----------------------------------------------------
  glb <- zone$glb
  ctrl <- zone$ctrl
  zoneC <- zone$zoneC
  zoneD <- zone$zoneD
  zone1 <- zone$zone1
  saudat <- zone$saudat
  constants <- zone$constants
  production <- zone$product
  projC <- zone$zone1$projC
  condC <- zone$zone1$condC



saunames <- glb$saunames

compdat <- condC$compdat$lfs

sau <- 4



saucompdata(allcomp=compdat, glb=glb, horizline = NULL,console = TRUE,
            rundir = "",barcol = "red",bordercol = "black",
            ylabel = "",tabname = "")















