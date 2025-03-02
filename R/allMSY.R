


#' @title allMSY estimates MSY by population, SAU, and Zone for each LML in data
#'
#' @description allMSY discovers all values for the LML (nlml) in the
#'     conditioning data and then estimates the productivity and MSY details
#'     at the population, SAU, and Zone scales. It gnerates a webpage output
#'     using makehtml.
#'
#' @param msydir full path to the directory in which the results are to be
#'     stored. Best to use confirmdir(msydir,ask=FALSE) to be sure it exists.
#' @param rundir the rundir from the original analysis
#' @param ctrlfile the name of the controlfile, which will be added to the
#'     rundir so the file can be read
#' @param uplim upper limit to the harvest rates to search for productivity
#' @param inc the increment from 0 up to uplim. The smaller this is the longer
#'     it will take.
#' @param verbose should comments and warnings go to the console. Default=TRUE
#'
#' @returns a list of the glb, ctrl, warning file path, and total time take
#'     (all required by makehtml), plus summarymsy a list of nlml lists of
#'     sauprod, zoneprod, and sauprodzone for each LML. Finally, fro each LML,
#'     inside outproduct is a copy of zoneC, zoneD, product, and glb.
#' @export
#'
#' @examples
#' print("no chance")
allMSY <- function(msydir,rundir,ctrlfile,uplim=0.4,inc=0.005,verbose=TRUE) {
  ##   msydir=msydir;rundir=rundir;ctrlfile=controlfile;uplim=0.9;inc=0.01;verbose=TRUE
  starttime <- Sys.time()
  setuphtml(msydir)
  zone1 <- readctrlfile(rundir,infile=ctrlfile,verbose=verbose)
  ctrl <- zone1$ctrl
  glb <- zone1$globals     # glb without the movement matrix
  bysau <- zone1$ctrl$bysau
  opar <- NULL
  parsin <- zone1$condC$parsin
  if (parsin) opar <- as.matrix(zone1$condC$optpars)
  if (is.null(bysau)) bysau <- 0
  if (bysau) {
    saudata <- readsaudatafile(rundir,ctrl$datafile,optpar=opar,verbose=verbose)
    constants <- saudata$constants
    saudat <- saudata$saudat
    zone1$condC$poprec <- saudata$poprec
  } else {
    constants <- readpopdatafile(rundir,ctrl$datafile)
    saudat <- constants
  }
  if (verbose) cat("Files read, now making zone \n")
  ans <- makezoneC(zone1,constants) # initiates zoneC
  zoneC <- ans$zoneC
  glb <- ans$glb
  ans <- makezone(glob=glb,zoneC=zoneC) # make zoneD, add cpue, qest to zoneC
  zoneC <- ans$zoneC  # zone constants
  zoneD <- ans$zoneD  # zone dynamics
  product <- NULL
  if (verbose)
    cat("Now estimating population productivity between H ",inc," and ",
        uplim, "\n")
  zoneC1 <- zoneC[[1]]
  lml <- sort(unique(zoneC1$LML))
  nlml <- length(lml)
  outz <- makelist(lml)
  outmsy <- makelist(lml)
  summarymsy <- makelist(lml)
  for (i in 1:nlml) { # i = 1
    cat("Using LML = ",lml[i],"\n")
    pickyr <- which(zoneC1$LML == lml[i])
    ans <- modzoneC(zoneC=zoneC,zoneD=zoneD,glob=glb,selectyr=pickyr[1],
                    uplim=uplim,inc=inc)
    zoneC <- ans$zoneC  # zone constants
    product <- ans$product  # productivity by population
    outz[[i]] <- list(zoneC=zoneC, zoneD=zoneD, product=product,glb=glb)
  }
  for (i in 1:nlml) { # i = 1 # tables of app productivity
    product <- outz[[i]]$product
    msymat <- (findmsy(product))
    total <- numeric(7)
    total[c(1,2,4)] <- colSums(msymat[,c(1,2,4)])
    msymat1 <- rbind(msymat,total)
    filen <- paste0("production_by_app_for_LML_",lml[i],".csv")
    caption <- paste0("Production by app for LML = ",lml[i],". Total ",
                      "estimate of MSY (catch) is maximum possible if all SAU",
                      "fished at approriate harest rates rather than a zone",
                      "-wide harvets rate.")
    addtable(msymat1,filen,rundir=msydir,category="population",caption)
  }
  for (i in 1:nlml) { # i = 1  # sau productivity
    product <- outz[[i]]$product
    if (verbose) cat("LML = ",lml[i],"\n")
    nsau <- glb$nSAU
    npop <- glb$numpop
    nh <- dim(product)[1]
    label <- dimnames(product)
    harv <- as.numeric(label[[1]])
    sauindex <- glb$sauindex
    saunames <- glb$saunames
    wts <- matrix(0,nrow=nh,ncol=npop,dimnames=list(label[[1]],label[[3]]))
    sauyield <- matrix(0,nrow=nh,ncol=nsau,dimnames=list(label[[1]],saunames))
    saumatB <- matrix(0,nrow=nh,ncol=nsau,dimnames=list(label[[1]],saunames))
    sauexB <- saumatB
    saucpue <- saumatB
    # Now do sau production
    for (h in 1:nh) {
      saumatB[h,] <- tapply(product[h,"MatB",],sauindex,sum,na.rm=TRUE)
      sauexB[h,] <- tapply(product[h,"ExB",],sauindex,sum,na.rm=TRUE)
      sauyield[h,] <- tapply(product[h,"Catch",],sauindex,sum,na.rm=TRUE)
      wts[h,] <- product[h,"Catch",]/sauyield[h,sauindex]
      saucpue[h,] <- tapply((product[h,"RelCE",] * wts[h,]),sauindex,sum,na.rm=TRUE)
    }
    label <- glb$saunames
    rows <- c("B0","Bmsy","MSY","Dmsy","CEmsy","Hmsy","Bexmsy")
    sauprod <- matrix(0,nrow=length(rows),ncol=nsau,dimnames=list(rows,saunames))
    filen <- paste0("Production_by_SAU_vs_CPUE_for_LML_",lml[i],".png")
    fullfilen <- pathtopath(msydir,filen)
    plotprep(width=8,height=7,newdev=FALSE,filename=fullfilen,verbose=FALSE)
    parset(plots=pickbound(nsau),margin=c(0.3,0.3,0.05,0.05),outmargin=c(1,1,0,0))
    for (sau in 1:nsau) { # sau=1
      ymax <- getmax(sauyield[,sau])
      xmax <- getmax(saucpue[2:nh,sau])
      xmin <- getmin(saucpue[2:nh,sau])
      pick <- which.max(sauyield[,sau])
      msyce <- saucpue[pick,sau]
      plot(saucpue[2:nh,sau],sauyield[2:nh,sau],type="l",lwd=2,xlab="",ylab="",
           panel.first=grid(),ylim=c(0,ymax),yaxs="i",xlim=c(xmin,xmax))
      abline(v=msyce,col=2,lwd=2)
      msylab <- paste0("MSY = ",round(sauyield[pick,sau],1))
      text(0.7*max(saucpue[,sau]),0.92*ymax,msylab,cex=1.2,pos=4)
      text(0.8*max(saucpue[,sau]),0.75*ymax,label[sau],cex=1.5,pos=4)
      text(1.1*msyce,0.15*ymax,round(msyce,2),cex=1.25,pos=4)
      sauprod[,sau] <- c(saumatB[1,sau],saumatB[pick,sau],sauyield[pick,sau],
                         (saumatB[pick,sau]/saumatB[1,sau]),saucpue[pick,sau],
                         harv[pick],sauexB[pick,sau])
    }
    addtxt <- paste0("CPUE at Bmsy (note different scales), LML = ",lml[i])
    mtext(addtxt,side=1,outer=TRUE,cex=1.0,
          line = -0.1)
    mtext("Equilibrium Yield  (note different scales)",side=2,outer=TRUE,cex=1.0,
          line=-0.2)
    caption <- paste0("The equilibrium production vs CPUE for each SAU. The red ",
                      "vertical lines and related number represent the expected",
                      "CPUE when biomass is at Bmsy. ")
    addplot(filen,rundir=msydir,category="SAU",caption)
    # now whole of zone for lml[i]
    yield <- rowSums(product[,"Catch",])
    spb <- rowSums(product[,"MatB",])
    expB <- rowSums(product[,"ExB",])
    Ht <- harv
    depletMSY <- spb/spb[1]
    maxy <- getmax(yield)
    pickmsy <- which.max(yield)
    filen <- pathtopath(msydir,paste0("production_SpB_Zone","_LML_",lml[i],
                                      ".png"))
    plotprep(width=7,height=6,newdev=FALSE,filename=filen,verbose=FALSE)
    parset(plots=c(3,2),cex=0.9,outmargin=c(0,1.0,0,0))
    plot(spb,yield,type="l",lwd=2,col=1,xlab="Spawning Biomass t",
         ylab="",panel.first = grid(),
         ylim=c(0,maxy),yaxs="i")
    abline(v=spb[pickmsy],col=2,lwd=2)

    plot(Ht,spb,type="l",lwd=2,xlab="Annual Harvest Rate",
         ylab=paste0("SpBiom t","_LML_",lml[i]),panel.first = grid(),
         ylim=c(0,getmax(spb)),yaxs="i")
    abline(h=spb[pickmsy],col=2,lwd=2)
    abline(v=Ht[pickmsy],col=2,lwd=2)

    plot(Ht,yield,type="l",lwd=2,col=1,xlab="Annual Harvest Rate",
         ylab="",panel.first = grid(),
         ylim=c(0,maxy),yaxs="i")
    abline(v=Ht[pickmsy],col=2,lwd=2)

    plot(spb,depletMSY,type="l",lwd=2,xlab="Spawning Biomass t",
         ylab=paste0("Depl Level","_LML_",lml[i]),
         panel.first = grid(),ylim=c(0,1.05),yaxs="i")
    abline(h=depletMSY[pickmsy],col=2,lwd=2)
    abline(v=spb[pickmsy],col=2,lwd=2)

    plot(depletMSY,yield,type="l",lwd=2,col=1,xlab="Total Depletion Level",
         ylab="",panel.first = grid())
    abline(v=depletMSY[pickmsy],col=2,lwd=2)

    plot(Ht,depletMSY,type="l",lwd=2,col=1,xlab="Annual Harvest Rate",
         ylab=paste0("Depl Level","_LML_",lml[i]),
         panel.first = grid(),ylim=c(0,1.05),yaxs="i")
    abline(h=depletMSY[pickmsy],col=2,lwd=2)
    abline(v=Ht[pickmsy],col=2,lwd=2)
    mtext(paste0("Production t","_LML_",lml[i]),side=2,outer=TRUE,cex=1.1,
          line=-0.2)
    caption <- paste0("The production curves for the zone. Also the ",
                      "relationships between spawning biomass depletion and ",
                      "harvest rate.")
    addplot(filen,rundir=msydir,category="Zone",caption)
    # zone outputs
    zoneprod <- cbind(catch=yield,matB=spb,harv=Ht,deplsB=depletMSY,
                      expB=expB)
    # add zone production to sauprod
    zonemsy <- zoneprod[which.max(zoneprod[,"catch"]),]
    Pzone <- numeric(nrow(sauprod)); names(Pzone) <- rownames(sauprod)
    Pzone[c(1:4,6:7)] <- c(B0=zoneprod[1,"matB"],Bmsy=zonemsy["matB"],
                           MSY=zonemsy["catch"],Dmsy=zonemsy["expB"],
                           Hmsy=zonemsy["harv"],Bexmsy=zonemsy["expB"])
    wts <- sauprod["MSY",]/Pzone["MSY"]; wts <- wts/sum(wts)
    Pzone[5] <- sum(sauprod["CEmsy",]*wts)
    # Pzone[c(1:3,6:7)] <- rowSums(sauprod[c(1:3,7),],na.rm=TRUE)
    #
    #  Pzone[4:6] <- rowSums(sauprod[4:6,]*wts)
    sauprodzone <- cbind(sauprod,zone=Pzone)
    # add sauprod table to production tab
    filen <- paste0("Production_by_sau_for_LML_",lml[i],".csv")
    caption <- paste0("Productivity properties: B0, Bmsy, MSY, Dmsy, Cemsy,",
                      " Hmsy, and Bexmsy for each sau and zone. CEmsy is the",
                      " predicted CPUE by SAU at MSY weighted by the sauMSY. ",
                      "They must,therefore be treated as approximate.")
    addtable(sauprodzone,filen,rundir=msydir,category="Zone",caption)
    # add zone production table
    filen <- paste0("production_across_zone_for_LML_",lml[i],".csv")
    caption <- paste0("Zone production: Catch, matureB, HarvestR, and Mature ",
                      "depletion level. MSY at H = ",Ht[pickmsy])
    addtable(zoneprod,filen,rundir=msydir,category="Zone",caption)
    summarymsy[[i]] <- list(sauprod=sauprod,zoneprod=zoneprod,
                            sauprodzone=sauprodzone)
    endtime <- Sys.time()
    tottime <- endtime - starttime
    if (verbose) print(tottime)
  }
  outmsy <- list(glb=glb,ctrl=ctrl,warn=glb$warnfile,tottime=tottime,
                 summarymsy=summarymsy,outproduct=outz)
  return(invisible(outmsy))
} # end of allMSY
