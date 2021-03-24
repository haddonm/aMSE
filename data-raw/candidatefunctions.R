



B0 <- getvar(zoneC,"B0")
ExB0 <- getvar(zoneC,"ExB0")
zoneDsau <- zonetosau(zoneDDR,glb,B0,ExB0)



plotNt <- function(Nt,year,start=3,glb,newdev=FALSE) {
#   Nt <- zoneDsau$Nt; origNt=zoneD$Nt;year=47;start=3; glb=glb; newdev=FALSE
  Nclass <- glb$Nclass
  midpts <- glb$midpts
  sc <- start:Nclass
  sizes <- midpts[sc]
  dat <- Nt[sc,,,]
  reps <- dim(dat)[4]
  nsau <- glb$nSAU
  saunames <- glb$saunames
  plotprep(width=8,height=8,newdev=newdev)
  parset(plots=c(4,2),byrow=FALSE)
  for (sau in 1:nsau) { #  sau=1
    sdat <- dat[,year,sau,]
    ymax <- getmax(sdat)
    plot(sizes,sdat[,1],type="l",lwd=1,col="grey",panel.first=grid(),
         ylim=c(0,ymax),ylab=saunames[sau])
    for (i in 1:reps) lines(sizes,sdat[,i])
    lines(sizes,dat[,sau,2,1],lwd=2,col=2)
  }


}

plotNt(zoneDsau$Nt,year=47,start=3,glb=glb,newdev=FALSE)


