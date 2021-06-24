

# prepare size-composition data ------------------------------------------------

B0 <- getvar(zoneC,"B0")
ExB0 <- getvar(zoneC,"ExB0")
zoneDsau <- zonetosau(zoneDDR,glb,B0,ExB0)
zonePsau <- zonetosau(zoneDP,glb,B0,ExB0)


plotprep(width=8,height=8,newdev=FALSE)
x <- plotNt(zonePsau$Nt,year=30,glb=glb,start=3,medcol=1)
str(x)


lcomp <- prepareDDNt(zoneDD$Nt,zoneDD$catchN,glb)

plotCNt(lcomp$Nt,glb,vline=c(140),start=3) # plot the conditioning history

plotCNt(zonePsau$Nt[,,,51],glb,vline=140,start=3) # plot a single replicate from projections


zonePsau <- zonetosau(zoneDP,glb,B0,ExB0)

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



outLF <- getLFdata(datadir,"lateLF-84-20.csv")

palfs=outLF$palfs
palfs





# copyto -------------------------------------------------------------


rundir <- "c:/Users/User/DropBox/A_codeUse/aMSEUse/scenarios/HS652510"

destdir <- "c:/Users/User/DropBox/A_codeUse/aMSEUse/scenarios/HS81"

copyto(rundir,todir=destdir,filename="controlsau.csv")

  #     zoned=zoneDsau;zonep=zonePsau;glb=glb;startyr=30;
  #     picksau=9; histCE=histCE;CIprobs=c(0.05,0.5,0.95); addCI=TRUE


# Diagnostic plots --------------------------------------------------------



invar <- zonePsau$catch
invar2 <- zonePsau$acatch
filen <- ""
nsau=8
nrep=3
reps=100
saunames=glb$saunames
label <- "_catches"

nyrs <- dim(invar)[1]
reps <- dim(invar)[3]
plotprep(width=7,height=7,filename=filen,cex=1.0,verbose=FALSE,newdev = FALSE)
parset(plots=c(4,2))
for (sau in 1:nsau) {
  pickrep <- sample(1:reps,nrep,replace=FALSE)
  ymax <- getmax(invar[,sau,pickrep])
  ylabel <- paste0(paste0("SAU_",saunames[sau],label))
  plot(1:nyrs,invar[,sau,pickrep[1]],type="l",lwd=2,ylim=c(0,ymax),
       ylab=ylabel,xlab="Years",panel.first=grid())
  for (i in 2:nrep) lines(1:nyrs,invar[,sau,pickrep[i]],lwd=2,col=i)
  for (i in 1:nrep) lines(1:nyrs,invar2[,sau,pickrep[i]],lwd=2,col=i,lty=3)
} # end of actual catches


# selectivity plots -----------------------------------------------------

pop1 <- zoneC[[1]]
sel <- pop1$Select
mids <- glb$midpts
rge <- 50:105

plotprep(width=7,height=3.5,newdev=FALSE)
parset()
plot(mids[rge],sel[rge,14],type="l",lwd=2,col=1,panel.first=grid())
lines(mids[rge],sel[rge,15],lwd=2,col=2)


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






# digits by row ----------------------------------------------------------------
library(knitr)




x <- matrix(c(rnorm(5,mean=5,sd=1),seq(1,10,1)),nrow=3,ncol=5,byrow=TRUE,
            dimnames=list(1:3,1:5))
digitsbyrow(x, c(3,0,0))



digitsbyrow
function(df, digits) {

  df <- propDD
  digits=c(0,2,3,3,3,3,2,3,3,3,3,3,0,3,3,1)

  tmp0 <- data.frame(t(df))






  tmp1 <- mapply(
    function(df0, digits0) {
      formatC(df0, format="f", digits=digits0)
    },
    df0=tmp0, digits0=digits
  )
  tmp1 <- data.frame(t(tmp1))
  rownames(tmp1) <- rownames(df)
  colnames(tmp1) <- colnames(df)
  if (class(df)[1] == "matrix") tmp1 <- as.matrix(tmp1)
  return(tmp1)
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


