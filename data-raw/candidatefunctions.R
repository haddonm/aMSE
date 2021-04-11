



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



plotprep(width=10,height=6,newdev=FALSE)

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




destdir <- "c:/Users/Malcolm/DropBox/A_code/aMSEUse/scenarios/tasHS653510"

copyto(rundir,todir=destdir,filename="controlsau.csv")

  #     zoned=zoneDsau;zonep=zonePsau;glb=glb;startyr=30;
  #     picksau=9; histCE=histCE;CIprobs=c(0.05,0.5,0.95); addCI=TRUE




plotprep(width=12,height=8,newdev=FALSE)
doonesau(prerep=zoneDsau,postrep=zonePsau,glb=glb,startyr=30,picksau=11,
         addCI=TRUE,histCE=histCE)




plotprep(width=8,height=8,newdev=FALSE)
dosau(sauZone,glb,picksau=12,histCE=condC$histCE,yrnames=1973:2019)








