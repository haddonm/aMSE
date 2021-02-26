



addrecvar2 <- function(zoneDD,zoneC,glob,condC,ctrl,varyrs,
                       sigR=1e-08,sigB=1e-08,lastsigR=0.3) {
  zoneC=zone$zoneC; zoneDD=zoneDD;glob=zone$glb;condC=zone$zone1$condC
  sigR=1e-08; sigB=1e-08; lastsigR=0.3; reps=ctrl$reps; varyrs=7

  histC <- condC$histCatch
  yrs <- condC$histyr[,"year"]
  nyrs <- length(yrs)
  zoneDDR <- makezoneDP(nyrs,reps,glob)
  finalyr <- nyrs - varyrs
  zoneDDR$matureB[1:finalyr,,] <- zoneDD$matureB[1:finalyr,]
  zoneDDR$exploitB[1:finalyr,,] <- zoneDD$exploitB[1:finalyr,]
  zoneDDR$catch[1:finalyr,,] <- zoneDD$catch[1:finalyr,]
  zoneDDR$harvestR[1:finalyr,,] <- zoneDD$harvestR[1:finalyr,]
  zoneDDR$cpue[1:finalyr,,] <- zoneDD$cpue[1:finalyr,]
  zoneDDR$recruit[1:finalyr,,] <- zoneDD$recruit[1:finalyr,]
  zoneDDR$catchN[,1:finalyr,,] <- zoneDD$catchN[,1:finalyr,]
  zoneDDR$Nt[,1:finalyr,,] <- zoneDD$Nt[,1:finalyr,]
  for (iter in 1:reps) {
    for (year in (finalyr+1):nyrs) {
      catchsau <- histC[year,]
      exb <- zoneDDR$exploitB[year-1,,iter]
      inN <- zoneDDR$Nt[,year-1,,iter]
      out <- oneyearsauC(zoneCC=zoneC,exb=exb,inN=inN,catchsau=catchsau,
                         year=year,Ncl=glob$Nclass,sauindex=glob$sauindex,
                         movem=glob$move,sigmar=lastsigR,sigmab=sigB)
      dyn <- out$dyn
      zoneDDR$exploitB[year,,iter] <- dyn["exploitb",]
      zoneDDR$matureB[year,,iter] <- dyn["matureb",]
      zoneDDR$catch[year,,iter] <- dyn["catch",]
      zoneDDR$harvestR[year,,iter] <- dyn["catch",]/out$dyn["exploitb",]
      zoneDDR$cpue[year,,iter] <- dyn["cpue",]
      zoneDDR$recruit[year,,iter] <- dyn["recruits",]
      zoneDDR$deplsB[year,,iter] <- dyn["deplsB",]
      zoneDDR$depleB[year,,iter] <- dyn["depleB",]
      zoneDDR$Nt[,year,,iter] <- out$NaL
      zoneDDR$catchN[,year,,iter] <- out$catchN
    }
    endcatch <- tapply(zoneDDR$catch[nyrs,,1],glb$sauindex,sum,na.rm=TRUE)
    yrce <- nrow(multTAC)
    acatch <- endcatch * multTAC[yrce,]  # predicted aspirational catches
    sigmar=ctrl$withsigR
    sigmab=ctrl$withsigB
    sauindex <- glob$sauindex
    for (iter in 1:reps) {
      exb=zoneDDR$exploitB[nyrs,,iter]
      inN=zoneDDR$Nt[,nyrs,,iter]
      outy <- oneyearsauC(zoneCC=zoneCP,exb=exb,inN=inN,catchsau=acatch,year=1,
                          Ncl=glb$Nclass,sauindex=sauindex,movem=glb$move,
                          sigmar=sigmar,sigmab=sigmab)
      dyn <- outy$dyn
      saudyn <- popcetosauce(dyn["catch",],dyn["cpue",],sauindex)
      zoneDP$exploitB[1,,iter] <- dyn["exploitb",]
      zoneDP$matureB[1,,iter] <- dyn["matureb",]
      zoneDP$catch[1,,iter] <- dyn["catch",]
      zoneDP$acatch[1,,iter] <- acatch
      zoneDP$catsau[1,,iter] <- saudyn$saucatch
      zoneDP$harvestR[1,,iter] <- dyn["catch",]/dyn["exploitb",]
      zoneDP$cpue[1,,iter] <- dyn["cpue",]
      zoneDP$cesau[1,,iter] <- saudyn$saucpue
      zoneDP$recruit[1,,iter] <- dyn["recruits",]
      zoneDP$deplsB[1,,iter] <- dyn["deplsB",]
      zoneDP$depleB[1,,iter] <- dyn["depleB",]
      zoneDP$Nt[,1,,iter] <- outy$NaL
      zoneDP$catchN[,1,,iter] <- outy$catchN
    }
    return(zoneDP)
  }

} # end of addrecvar2




round(zoneDDR$exploitB[40:47,,1],2)

round(zoneDDR$Nt[,nyrs,1,1:5])
























