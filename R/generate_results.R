  # These plots use both zoneC and zoneD
  # characterize productivity ----------------------------------------


#' @title biology_plots generates a series of stored plots and tables
#'
#' @description biology_plots generates the yield vs spawning biomass,
#'     weight-at-lenght, and emergence plots for all populations. In
#'     addition, it also tabulates the biological properties of each
#'     population and SAU and the total zone
#'
#' @param resdir the results directory
#' @param runlabel the runlabel, ideally from the ctrl object
#' @param glb the globals list
#' @param zoneC the zonal constants by population, zoneC
#'
#' @return nothing but it does plot 3 plots and put one table into
#'     resdir
#' @export
#'
#' @examples
#' print("this will be quite long when I get to it")
biology_plots <- function(resdir, runlabel, glb, zoneC) {
  #      resdir=resdir; zoneC=zoneC; zoneD=zoneD; product=product; glb=glb; runlabel=ctrl$runlabel
  # some globals
  mids <- glb$midpts
  numpop <- glb$numpop
  resfile <- paste0(resdir,"/resultTable_",runlabel,".csv")
  # Yield vs Spawning biomass-----------------------
  # maturation uses zoneC
  matur <- getlistvar(zoneC,"Maturity")
  rownames(matur) <- mids
  file <- paste0("maturity_v_Length_",runlabel,".png")
  filename <- filenametopath(resdir,file)  #  filename=""
  plotprep(width=7,height=4,newdev=FALSE,filename=filename,cex=0.9,
           verbose=FALSE)
  plot(mids,matur[,1],type="l",lwd=2,xlab="Shell Length mm",
       ylab="Proportion Mature",panel.first=grid(),xlim=c(50,210))
  for (pop in 2:numpop) lines(mids,matur[,pop],lwd=2,col=pop)
  legend("topright",paste0("P",1:numpop),lwd=3,col=c(1:numpop),bty="n",
         cex=1.2)
  if (nchar(filename) > 0) dev.off()
  caption <- "The maturity vs length for each population."
  logfilename(filename,resfile=resfile,"Biology",caption)

  # weight-at-length using zoneC--------------------------------------
  WtL <- getlistvar(zoneC,"WtL")
  rownames(WtL) <- mids
  file <- paste0("Weight_at_Length_",runlabel,".png")
  filename <- filenametopath(resdir,file)  #  filename=""
  plotprep(width=7,height=4,newdev=FALSE,filename=filename,cex=0.9,
           verbose=FALSE)
  plot(mids,WtL[,1],type="l",lwd=2,xlab="Shell Length mm",
       ylab="Weight Kg",panel.first=grid(),xlim=c(110,210))
  for (pop in 2:numpop) lines(mids,WtL[,pop],lwd=2,col=pop)
  legend("topleft",paste0("P",1:numpop),lwd=3,col=c(1:numpop),bty="n",
         cex=1.2)
  if (nchar(filename) > 0) dev.off()
    caption <- paste0("The weight-at-length for each population. ",
                      "The x-axis is constrained to encompass legal sizes.")
  logfilename(filename,resfile=resfile,"Biology",caption)

  # emergence uses zoneC----------------------------------------------
  emerg <- getlistvar(zoneC,"Emergent")
  rownames(emerg) <- mids
  file <- paste0("Emergence_at_Length_",runlabel,".png")
  filename <- filenametopath(resdir,file)  #  filename=""
  plotprep(width=7,height=4,newdev=FALSE,filename=filename,cex=0.9,
           verbose=FALSE)
  plot(mids,emerg[,1],type="l",lwd=2,xlab="Shell Length mm",
       ylab="Weight Kg",panel.first=grid(),xlim=c(110,145))
  for (pop in 2:numpop) lines(mids,emerg[,pop],lwd=2,col=pop)
  legend("topleft",paste0("P",1:numpop),lwd=3,col=c(1:numpop),bty="n",
         cex=1.2)
  if (nchar(filename) > 0) dev.off()

  caption <- "The emergence-at-length for each population. The x-axis is constrained to emphasize differences."
  logfilename(filename,resfile=resfile,"Biology",caption)


  # Tabulate biological properties uses zoneC-------------------------
  getvar <- function(indexvar,zone=zoneC) { # indexvar="B0"; zoneC=zoneC
    sau <- getlistvar(zone,"SAU")
    nSAU <- length(unique(sau))
    saures <- numeric(nSAU)
    invar <- getlistvar(zone,indexvar)
    for (mu in 1:nSAU) {
      pickcol <- which(sau == mu)
      saures[mu] <- sum(invar[pickcol])
    }
    ans <- c(invar,saures,sum(invar))
    return(ans)
  }
  rows <- c("SAU","M","R0","B0","ExB0","MSY","MSYDepl","bLML")
  sau <- getlistvar(zoneC,"SAU")
  nSAU <- length(unique(sau))
  columns <- c(paste0("p",1:numpop),
               paste0("sau",1:nSAU),"zone")
  numrow <- length(rows)
  numcol <- length(columns)
  results <- matrix(0,nrow=numrow,ncol=numcol,
                                  dimnames=list(rows,columns))
  results["SAU",] <- c(getlistvar(zoneC,"SAU"),(1:nSAU),0) # no total
  results["B0",] <- as.numeric(getvar("B0"))
  M <- getlistvar(zoneC,"Me")
  wtr <- (results["B0",1:numpop]/results["B0",(numpop+nSAU+1)])

  results["M",] <- c(M,0,0,sum(M*wtr))
  results["R0",] <- getvar("R0")
  results["ExB0",] <- getvar("ExB0")
  results["MSY",] <- getvar("MSY")
  MSYD <- getlistvar(zoneC,"MSYDepl")
  results["MSYDepl",] <- c(MSYD,0,0,sum(MSYD * wtr))
  bLML <- getlistvar(zoneC,"bLML")
  results["bLML",] <- c(bLML,0,0,sum(bLML * wtr))
  for (mu in 1:nSAU) { # mu = 2
    pickcol <- which(sau == mu)
    wtmu <- length(pickcol)
    wtmu <- results["B0",pickcol]/results["B0",(numpop + mu)]
    results["M",(numpop + mu)] <- sum(M[pickcol]*wtmu)
    results["MSYDepl",(numpop + mu)] <- sum(MSYD[pickcol]*wtmu)
    results["bLML",(numpop + mu)] <- sum(bLML[pickcol]*wtmu)
  }
  res <- round(results,3)
  filen <- paste0("zonebiology_",runlabel,".csv")
  filename <- filenametopath(resdir,filen)
  write.table(res,file = filename,sep=",")
  caption <- paste0("Population SAU and Zonal Biological Properties. ",
              "For 'M' 'MSYDepl' and 'bLML' SAU and total is an average ",
              "weighted relative to the proportion of SAU or total B0.")
  logfilename(filename,resfile=resfile,"Tables",caption)


} # end of biology_plots


#' @title numbersatsize plots details of the numbers-at-size
#'
#' @description numbersatsize plots up the initial unfished numbers-
#'     at-size distribution, omitting the first four size classes to
#'     avoid the recruitment numbers dominating the plot.
#'
#' @param resdir the results directory
#' @param runlabel the runlabel, ideally from the ctrl object
#' @param glb the globals list
#' @param zoneD the dynamic part of the zone, zoneD
#'
#' @return nothing but it does add one plot to the results directory
#' @export
#'
#' @examples
#' print("this will be quite long when I get to it")
numbersatsize <- function(resdir, runlabel, glb, zoneD) {
  #      resdir=resdir; zoneD=zoneD; glb=glb; runlabel=ctrl$runlabel
  # some globals
  mids <- glb$midpts
  numpop <- glb$numpop
  resfile <- paste0(resdir,"/resultTable_",runlabel,".csv")
  # initial numbers-at-size uses zoneD--------------------------------
  Nt <- zoneD$Nt[,1,]/1000.0
  Ntt <- rowSums(zoneD$Nt[,1,])/1000.0  # totals
  file <- paste0("Total_Initial_Numbers-at-Size_",runlabel,".png")
  filename <- filenametopath(resdir,file)  #  filename=""
  plotprep(width=7,height=6,newdev=FALSE,filename=filename,cex=0.9,verbose=FALSE)
  parset(plots=c(2,1))
  plot(mids[5:105],Ntt[5:105],type="l",lwd=2,xlab="Shell Length mm (5 - 210mm)",
       ylab="Numbers-at_size '000s",panel.first=grid())
  maxy <- getmax(Nt[5:105,])
  plot(mids[5:105],Nt[5:105,1],type="l",lwd=2,xlab="Shell Length mm (5 - 210mm)",
       ylab="Numbers-at_size '000s",panel.first=grid(),ylim=c(0,maxy))
  for (pop in 2:numpop) lines(mids[5:105],Nt[5:105,pop],lwd=2,col=pop)
  abline(h=0.0,col="darkgrey")
  legend("topright",paste0("P",1:numpop),lwd=3,col=c(1:numpop),bty="n",
         cex=1.2)
  if (nchar(filename) > 0) dev.off()
    caption <- "The numbers-at-size for the whole zone and for each population separately. The recruitment numbers are omitted for clarity."
  logfilename(filename,resfile=resfile,"NumSize",caption)

} # end of numbersatsize


#' @title compzoneN compares numbers-at-size before/after depletion
#'
#' @description compzoneN generates a plot comparing the unfished
#'     numbers-at-size with those for a given level of depletion.
#'
#' @param unfN the unfished numbers-at-size from getzoneprops
#' @param curN the current numbers-at-size from getzoneprops
#' @param glb the global object
#' @param yr the year of the dynamics
#' @param depl the depletion level of the current n-a-s
#' @param LML the legal minimum length in teh comparison year
#' @param resdir the results directory, default = "" leading to no
#'     .png file, just a plot to the screen.
#'
#' @return invisibly the filename ready for logfilename
#' @export
#'
#' @examples
#' print("still to be developed")
#' # unfN=unfN; curN=depN;glb=glb; yr=1; depl=0.3993; LML=132; resdir=resdir
compzoneN <- function(unfN,curN,glb,yr,depl,LML=0,resdir="") {
  usecl=5:glb$Nclass
  mids <- glb$midpts
  filen <- ""
  if (nchar(resdir) > 0) {
    filen <- paste0("zone_numbers-at-size_yr",yr,".png")
    filen <- filenametopath(resdir,filen)
  }
  plotprep(width=7, height=4.5, filename=filen,verbose=FALSE)
  plot(mids[usecl],unfN[usecl,"zone"]/1000.0,type="l",lwd=2,
       panel.first=grid(),xlab="Shell Length (mm)",
       ylab="Numbers-at-Size '000s")
  lines(mids[usecl],curN[usecl,"zone"]/1000.0,lwd=2,col=2)
  if (LML > 0) abline(v=LML,col=1,lty=2)
  legend("topright",c("Unfished",depl),col=c(1,2),lwd=3,bty="n")
  if (nchar(filen) > 0) dev.off()
  return(invisible(filen))
} # end of compzoneN

#' @title plotproductivity characterizes each population's yield curve
#'
#' @description plotproductivity characterizes each population's yield
#'     curve, it also describes the total productivity of the zone.
#'
#' @param resdir the results directory
#' @param runlabel the runlabel, ideally from the ctrl object
#' @param product the productivity 3-D array
#' @param glb the globals list
#'
#' @return nothing but it does place a number of png files into resdir
#' @export
#'
#' @examples
#' print("this will be quite long when I get to it")
plotproductivity <- function(resdir,runlabel,product,glb) {
  #   product=product; glb=glb; runlabel=ctrl$runlabel; resdir=resdir
  # All these plots only use the product array
  xval <- findmsy(product)
  numpop <- glb$numpop
  resfile <- paste0(resdir,"/resultTable_",runlabel,".csv")
  # Yield vs Spawning biomass --------
  file <- paste0("production_SpB_",runlabel,".png")
  filename <- filenametopath(resdir,file)  #  filename=""
  plotprod(product,xname="MatB",xlab="Spawning Biomass t",
           ylab="Production t",filename = filename)
  caption <- paste0("The production curve relative to each population's ",
                    "spawning biomass. The vertical lines identify the ",
                    "Bmsy values.")
  logfilename(filename,resfile=resfile,"Production",caption)

  # Yield vs Annual Harvest Rate
  file <- paste0("production_AnnH_",runlabel,".png")
  filename <- filenametopath(resdir,file)  #  filename=""
  plotprod(product,xname="AnnH",xlab="Annual Harvest Rate",filename = filename)
  caption <- paste0("The production curve relative to the Annual ",
                    "Harvest Rate applied to each population. The ",
                    "vertical lines identify the Hmsy values.")
  logfilename(filename,resfile=resfile,"Production",caption)

  # plot of Yield vs population depletion
  file <- paste0("production_Deplet_",runlabel,".png")
  filename <- filenametopath(resdir,file)
  plotprod(product,xname="Deplet",xlab="Population Depletion Level",
           filename = filename,devoff=FALSE)
  for (pop in 1:numpop) abline(v=xval[pop,"Deplet"],lwd=2,col=pop)
  if (nchar(filename) > 0) dev.off()
  caption <- paste0("The production curve relative to the depletion ",
              "level of each population. The vertical lines identify ",
              "the Depletion level giving rise to the MSY.")
  logfilename(filename,resfile=resfile,"Production",caption)

  # plot of Yield vs population depletion but constrained to within
  # 0.2 and 0.35 levels, to illustrate nearly flat rpoduction curve
  # and more clearly identify the population depletion at MSY
  file <- paste0("production_Deplet_0.2_0.35_",runlabel,".png")
  filename <- filenametopath(resdir,file)
  plotprod(product,xname="Deplet",xlab="Population Depletion Level",
           xlimit=c(0.2,0.35),filename = filename,devoff=FALSE)
  for (pop in 1:numpop) abline(v=xval[pop,"Deplet"],lwd=2,col=pop)
  if (nchar(filename) > 0) dev.off()
  caption <- paste0("The production curve relative to the depletion ",
                    "level of each population. Here the x-axis is ",
                    "shortened to clarify the flatness of the production ",
                    "curve about the MSY points.")
  logfilename(filename,resfile=resfile,"Production",caption)

  # Now do total production --------------------------------------------
  yield <- rowSums(product[,"Catch",])
  spb <- rowSums(product[,"MatB",])
  Ht <- yield/spb
  depletMSY <- spb/spb[1]
  pickmsy <- which.max(yield)
  maxy <- getmax(yield)

  file <- paste0("production_SpB_Total_",runlabel,".png")
  filename <- filenametopath(resdir,file)  #  filename=""
  plotprep(width=7,height=6,newdev=FALSE,filename=filename,verbose=FALSE)
  parset(plots=c(3,2),cex=0.9)
  plot(spb,yield,type="l",lwd=2,col=1,xlab="Spawning Biomass t",
       ylab="Production t",panel.first = grid(),
       ylim=c(0,maxy),yaxs="i")
  abline(v=spb[pickmsy],col=2,lwd=2)

  plot(Ht,spb,type="l",lwd=2,xlab="Annual Harvest Rate",
       ylab="Spawning Biomass t",panel.first = grid(),
       ylim=c(0,getmax(spb)),yaxs="i")
  abline(h=spb[pickmsy],col=2,lwd=2)
  abline(v=Ht[pickmsy],col=2,lwd=2)

  plot(Ht,yield,type="l",lwd=2,col=1,xlab="Annual Harvest Rate",
       ylab="Production t",panel.first = grid(),
       ylim=c(0,maxy),yaxs="i")
  abline(v=Ht[pickmsy],col=2,lwd=2)

  plot(spb,depletMSY,type="l",lwd=2,ylab="Total Depletion Level",
       xlab="Spawning Biomass t",panel.first = grid(),
       ylim=c(0,1.05),yaxs="i")
  abline(h=depletMSY[pickmsy],col=2,lwd=2)
  abline(v=spb[pickmsy],col=2,lwd=2)

  plot(depletMSY,yield,type="l",lwd=2,col=1,xlab="Total Depletion Level",
       ylab="Production t",panel.first = grid())
  abline(v=depletMSY[pickmsy],col=2,lwd=2)

  plot(Ht,depletMSY,type="l",lwd=2,col=1,xlab="Annual Harvest Rate",
       ylab="Total Depletion Level",panel.first = grid(),
       ylim=c(0,1.05),yaxs="i")
  abline(h=depletMSY[pickmsy],col=2,lwd=2)
  abline(v=Ht[pickmsy],col=2,lwd=2)
  if (nchar(filename) > 0) dev.off()
  caption <- paste0("The production curves for the zone. Also the ",
              "relationships between spawning biomass depletion and ",
              "harvest rate.")
  logfilename(filename,resfile=resfile,"Production",caption)


} # end of plotproductivity
