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
#' @param glb the globals list
#' @param zoneC the zonal constants by population, zoneC
#'
#' @return invisibly returns the biological properties of the populations
#' @export
#'
#' @examples
#' print("this will be quite long when I get to it")
biology_plots <- function(resdir, glb, zoneC) {
  mids <- glb$midpts
  numpop <- glb$numpop
  popdef <- getlistvar(zoneC,"popdef")
  SAU <- unique(popdef["SAU",])
  nSAU <- length(SAU)
  # Yield vs Spawning biomass-----------------------
  # maturation uses zoneC
  matur <- getlistvar(zoneC,"Maturity")
  rownames(matur) <- mids
  filen <- filenametopath(resdir,"maturity_v_Length_pop.png")
  plotprep(width=7,height=4,newdev=FALSE,filename=filen,cex=0.9,
           verbose=FALSE)
  plot(mids,matur[,1],type="l",lwd=2,xlab="Shell Length mm",
       ylab="Proportion Mature",panel.first=grid(),xlim=c(50,210))
  for (pop in 2:numpop) lines(mids,matur[,pop],lwd=2,col=pop)
  legend("topright",paste0("P",1:numpop),lwd=3,col=c(1:numpop),bty="n",
         cex=1.2)
  caption <- "The maturity vs length for each population."
  addplot(filen,resdir=resdir,category="Biology",caption)

  # weight-at-length using zoneC--------------------------------------
  WtL <- getlistvar(zoneC,"WtL")
  rownames(WtL) <- mids
  filen <- filenametopath(resdir,"Weight_at_Length_pop.png")
  plotprep(width=7,height=4,newdev=FALSE,filename=filen,cex=0.9,
           verbose=FALSE)
  plot(mids,WtL[,1],type="l",lwd=2,xlab="Shell Length mm",
       ylab="Weight Kg",panel.first=grid(),xlim=c(110,210))
  for (pop in 2:numpop) lines(mids,WtL[,pop],lwd=2,col=pop)
  legend("topleft",paste0("P",1:numpop),lwd=3,col=c(1:numpop),bty="n",cex=1.2)
  caption <- paste0("The weight-at-length for each population. ",
                      "The x-axis is constrained to encompass legal sizes.")
  addplot(filen,resdir=resdir,category="Biology",caption)

  # emergence uses zoneC----------------------------------------------
  emerg <- getlistvar(zoneC,"Emergent")
  rownames(emerg) <- mids
  filen <- filenametopath(resdir,"Emergence_at_Length_pop.png")
  plotprep(width=7,height=4,newdev=FALSE,filename=filen,cex=0.9,
           verbose=FALSE)
  plot(mids,emerg[,1],type="l",lwd=2,xlab="Shell Length mm",
       ylab="Weight Kg",panel.first=grid(),xlim=c(110,145))
  for (pop in 2:numpop) lines(mids,emerg[,pop],lwd=2,col=pop)
  legend("topleft",paste0("P",1:numpop),lwd=3,col=c(1:numpop),bty="n",cex=1.2)
  caption <- "The emergence-at-length for each population. The x-axis is constrained to emphasize differences."
  addplot(filen,resdir=resdir,category="Biology",caption)
  # Tabulate biological properties uses zoneC-------------------------
  rows <- c("SAU","M","R0","B0","ExB0","MSY","MSYDepl","bLML",
            "MaxDL","L50","L95","AvRec","steep")
  sau <- getlistvar(zoneC,"SAU")
  nSAU <- length(unique(sau))
  columns <- paste0("p",1:numpop)
  numrow <- length(rows)
  numcol <- length(columns)
  resultpop <- matrix(0,nrow=numrow,ncol=numpop,
                                  dimnames=list(rows,columns))
  resultpop["SAU",] <- getlistvar(zoneC,"SAU") # no total
  resultpop["B0",] <- as.numeric(getlistvar(zoneC,"B0"))
 # wtr <- (results["B0",1:numpop]/results["B0",(numpop+nSAU+1)])
  resultpop["M",] <- getlistvar(zoneC,"Me")
  resultpop["R0",] <- getlistvar(zoneC,"R0")
  resultpop["ExB0",] <- getlistvar(zoneC,"ExB0")
  resultpop["MSY",] <- getlistvar(zoneC,"MSY")
  resultpop["MSYDepl",] <- getlistvar(zoneC,"MSYDepl")
  resultpop["bLML",] <- getlistvar(zoneC,"bLML")
  popdefs <- getlistvar(zoneC,"popdef")
  resultpop["MaxDL",] <- popdefs["DLMax",]
  resultpop["L50",] <- popdefs["L50",]
  resultpop["L95",] <- popdefs["L95",]
  resultpop["AvRec",] <- round(popdefs["AvRec",])
  resultpop["steep",] <- popdefs["steeph",]
  # for (mu in 1:nSAU) { # mu = 2
  #   pickcol <- which(sau == mu)
  #   wtmu <- length(pickcol)
  #   wtmu <- results["B0",pickcol]/results["B0",(numpop + mu)]
  #   results["M",(numpop + mu)] <- sum(M[pickcol]*wtmu)
  #   results["MSYDepl",(numpop + mu)] <- sum(MSYD[pickcol]*wtmu)
  #   results["bLML",(numpop + mu)] <- sum(bLML[pickcol]*wtmu)
  # }
  res <- round(resultpop,3)
  filen <- paste0("zonebiology.csv")
  caption <- "Population Biological Properties."
  addtable(res,filen,resdir=resdir,category="Tables",caption)
  return(invisible(resultpop))
} # end of biology_plots

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
  filen <- paste0("zone_n-at-size_yr",yr,".png")
  filen <- filenametopath(resdir,filen)
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

#' @title numbersatsize plots details of the numbers-at-size
#'
#' @description numbersatsize plots up the initial unfished numbers-
#'     at-size distribution, omitting the first four size classes to
#'     avoid the recruitment numbers dominating the plot.
#'
#' @param resdir the results directory, if set to "" then plot is sent to
#'     the console instead
#' @param glb the globals list
#' @param zoneD the dynamic part of the zone, zoneD
#'
#' @return nothing but it does add one plot to the results directory
#' @export
#'
#' @examples
#' print("this will be quite long when I get to it")
numbersatsize <- function(resdir, glb, zoneD) {
  # some globals
  mids <- glb$midpts
  numpop <- glb$numpop
  # initial numbers-at-size uses zoneD--
  Nt <- zoneD$Nt[,1,]/1000.0
  Ntt <- rowSums(zoneD$Nt[,1,])/1000.0  # totals
  if (nchar(resdir) > 0) {
    filen <- file.path(resdir,"Initial_N-at-Size.png")
  } else {
    filen <- ""
  }
  plotprep(width=7,height=6,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
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
  if (nchar(resdir) > 0) {
    caption <- paste0("The numbers-at-size for the whole zone and for each ",
                      "population separately. The recruitment numbers are",
                      " omitted for clarity.")
    addplot(filen,resdir=resdir,category="NumSize",caption)
  }
} # end of numbersatsize

#' @title numbersatsizeSAU plots the numbers-at-size for a given SAU
#'
#' @description numbersatsizeSAU plots the initial unfished numbers-
#'     at-size distribution for a given SAU, omitting the first four size
#'     classes to avoid the recruitment numbers dominating the plot.
#'
#' @param resdir the results directory, if set to "" then plot is sent to
#'     the console instead
#' @param glb the globals list
#' @param zoneC the constant part of the zpne structure
#' @param zoneD the dynamic part of the zone, zoneD
#' @param sau the SAU name to select for plotting. If multiple populations are
#'     contained in an SAU their numbers-at-size will be combined.
#' @param yr which year of numbers-at-size should be plotted?
#' @param defpar should the plot parameters be defined. Set to FALSE if
#'     numbersatsizeSAU is to be used to add a plot to a multiple plot.
#' @param exploit should exploitable numbers be plotted as well? default=TRUE
#' @param mature should mature numbers be plotted as well? default=TRUE
#' @param filename default='Numbers-at-Size_Year1.png', but can be changed
#'     to suit whichever year is used
#'
#' @return nothing but it adda a plot to the results directory or console
#' @export
#'
#' @examples
#' print("this will be quite long when I get to it")
numbersatsizeSAU <- function(resdir, glb, zoneC, zoneD, sau, yr=1,
                             defpar=TRUE, exploit=TRUE, mature=TRUE,
                             filename="Numbers-at-Size_Year1.png") {
  # some globals resdir=""; glb=glb; zoneD=zoneD; sau=8; defpar=TRUE
  mids <- glb$midpts
  picksau <- which(zoneD$SAU == as.character(sau))
  Nt <- as.matrix(zoneD$Nt[,yr,picksau]/1000.0)
  if (length(picksau) > 1) {
    Ntt <- as.matrix(rowSums(Nt,na.rm=TRUE))  # totals
  } else {
    Ntt <- Nt
  }
  if (nchar(resdir) > 0) {
    filen <- file.path(resdir,filename)
  } else {
    filen <- ""
  }
  if ((exploit) | (mature)) pickzC <- which(sapply(zoneC,"[[","SAU") == sau)
  maxy <- getmax(Ntt[5:105,])
  if (defpar)
    plotprep(width=7,height=4,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  plot(mids[5:105],Ntt[5:105],type="l",lwd=2,xlab="Shell Length mm (5 - 210mm)",
       ylab="Numbers-at_size '000s",panel.first=grid(),ylim=c(0,maxy))
  if (exploit)
    lines(mids[5:105],zoneC[[pickzC]]$Select[5:105,yr] * Ntt[5:105],col=4,lwd=2)
  if (mature)
    lines(mids[5:105],zoneC[[pickzC]]$Maturity[5:105] * Ntt[5:105],col=2,lwd=2)
  legend("topright",legend=as.character(sau),lwd=0,col=0,bty="n",cex=1.2)
  if (nchar(resdir) > 0) {
    addm <- ""; adde <- ""
    if (mature) addm <- " The red line is mature numbers-at-size. "
    if (exploit) adde <- " The blue line is exploitable numbers-at-size."
    caption <- paste0("The numbers-at-size for the SAU ",sau," the recruitment ",
                      "numbers are omitted for clarity.",addm,adde)
    addplot(filen,resdir=resdir,category="NumSize",caption)
  }
} # end of numbersatsizeSAU

#' @title plothistcatch plots the historical catches used for conditioning
#'
#' @description plothistcatch provides a visual representation of the catches
#'     taken in each SAU/population used in the conditioning of the
#'     Operating Model.
#'
#' @param zone1 The zonewide properties as in readzonefile(indir,ctrl$zonefile)
#' @param pops default = NULL. The plot can be restricted to a specific set
#'     of populations by providing the indices of the populations to be
#'     included. Thus if there were 8 populations the pops=c(1,2,3,8), would
#'     plot the first three and the last population.
#' @param resdir the results directory used by makehtml to store plots and
#'     tables. If set to "", the default, then it plots to the console
#' @param defpar define the plot parameters. Set to FALSE if using
#'     plothistcatch to add a plot to a multiple plot
#'
#' @return nothing, it adds a plot into resdir and modifies resultTable.csv
#' @export
#'
#' @examples
#' print("wait until I have altered the internals data sets")
plothistcatch <- function(zone1,pops=NULL,resdir="",defpar=TRUE) {
  # zone1=zone1; pops=c(1); defpar=FALSE
  glb <- zone1$globals
  if (length(pops) > 0) {
    histcatch <- as.matrix(zone1$histCatch[,pops])
    numpop <- length(pops)
    addlab <- paste0(pops,"_",collapse="")
    label <- paste0("p",pops)
    numpop <- length(pops)
  } else {
    histcatch <- zone1$histCatch
    addlab <- "all"
    label <- colnames(histcatch)
    numpop <- glb$numpop
  }
  histyr <- zone1$histyr
  yearCE <- zone1$yearCE
  yearC <- histyr[,"year"]
  pngfile <- paste0("Historical_Catches_",addlab,".png")
  if (nchar(resdir) > 0) {
    filen <- filenametopath(resdir,pngfile)
  } else { filen <- "" }
  if (defpar)
      plotprep(width=7,height=4,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  ymax <- getmax(histcatch)
  plot(yearC,histcatch[,1],type="l",lwd=2,xlab="Year",
       ylab="Catch (t)",panel.first=grid(),ylim=c(0,ymax))
  if (numpop > 1) for (pop in 2:numpop)
    lines(yearC,histcatch[,pop],lwd=2,col=pop)
  legend("topright",label,lwd=3,col=c(1:numpop),bty="n",cex=1.2)
  if (nchar(resdir) > 0) {
    caption <- "The historical catch used for conditioning for each population."
    addplot(filen,resdir=resdir,category="history",caption)
  }
} # end of plothistcatch

#' @title plothistCE plots the historical cpue used for conditioning
#'
#' @description plothistCE provides a visual representation of the cpue
#'     taken in each SAU used in the conditioning of the
#'     Operating Model.
#'
#' @param zone1 The zonewide properties as in readzonefile(indir,ctrl$zonefile)
#' @param pops default = NULL. The plot can be restricted to a specific set
#'     of populations by providing the indices of the populations to be
#'     included. Thus if there were 8 populations the pops=c(1,2,3,8), would
#'     plot the first three and the last population.
#' @param resdir the results directory used by makehtml to store plots and
#'     tables.
#'
#' @return nothing, it adds a plot into resdir and modifies resultTable.csv
#' @export
#'
#' @examples
#' print("wait until I have altered the internals data sets")
plothistCE <- function(zone1,pops=NULL,resdir="") {
  # zone1=zone1; pops=c(1,2,3,8);
  glb <- zone1$globals
  numpop <- glb$numpop
  if (length(pops) > 0) {
    histCE <- zone1$histCE[,pops]
    numpop <- length(pops)
    addlab <- paste0(pops,"_",collapse="")
  } else {
    histCE <- zone1$histCE
    addlab="all"
  }
  yearCE <- zone1$yearCE
  pngfile <- paste0("Historical_CPUE_",addlab,".png")
  filen <- filenametopath(resdir,pngfile)
  plotprep(width=7,height=4,newdev=FALSE,filename=filen,cex=0.9,
           verbose=FALSE)
  ymax <- getmax(histCE)
  plot(yearCE,histCE[,1],type="l",lwd=2,xlab="Year",
       ylab="Standardized CPUE",panel.first=grid(),ylim=c(0,ymax))
  if (numpop > 1) for (pop in 2:numpop)
    lines(yearCE,histCE[,pop],lwd=2,col=pop)
  label <- colnames(histCE)
  legend("topright",label,lwd=3,col=c(1:numpop),bty="n",cex=1.2)
  caption <- "The historical CPUE used for conditioning for each population."
  addplot(filen,resdir=resdir,category="history",caption)
} # end of plothistCE

#' @title plotproductivity characterizes each population's yield curve
#'
#' @description plotproductivity characterizes each population's yield curve, it
#'   also describes the total productivity of the zone.
#'
#' @param resdir the results directory
#' @param product the productivity 3-D array
#' @param glb the globals list
#'
#' @return nothing but it does place five png files into resdir
#' @export
#'
#' @examples
#' print("this will be quite long when I get to it")
plotproductivity <- function(resdir,product,glb) {
  # All these plots only use the product array
  xval <- findmsy(product)
  numpop <- glb$numpop
  # Yield vs Spawning biomass --------
  filen <- filenametopath(resdir,"production_SpB.png")
  plotprod(product,xname="MatB",xlab="Spawning Biomass t",
           ylab="Production t",filename = filen,devoff=FALSE)
  caption <- paste0("The production curve relative to each population's ",
                    "spawning biomass. The vertical lines identify the ",
                    "Bmsy values.")
  addplot(filen,resdir=resdir,category="Production",caption)

  # Yield vs Annual Harvest Rate
  filen <- filenametopath(resdir,"production_AnnH.png")
  plotprod(product,xname="AnnH",xlab="Annual Harvest Rate",filename = filen,
           devoff=FALSE)
  caption <- paste0("The production curve relative to the Annual ",
                    "Harvest Rate applied to each population. The ",
                    "vertical lines identify the Hmsy values.")
  addplot(filen,resdir=resdir,category="Production",caption)

  # plot of Yield vs population depletion
  filen <- filenametopath(resdir,"production_Deplet.png")
  plotprod(product,xname="Deplet",xlab="Population Depletion Level",
           filename = filen,devoff=FALSE)
  for (pop in 1:numpop) abline(v=xval[pop,"Deplet"],lwd=2,col=pop)
  caption <- paste0("The production curve relative to the depletion ",
              "level of each population. The vertical lines identify ",
              "the Depletion level giving rise to the MSY.")
  addplot(filen,resdir=resdir,category="Production",caption)

  # plot of Yield vs population depletion but constrained to within
  # 0.2 and 0.35 levels, to illustrate nearly flat rpoduction curve
  # and more clearly identify the population depletion at MSY
  filen <- filenametopath(resdir,"production_Deplet_0.2_0.35.png")
  plotprod(product,xname="Deplet",xlab="Population Depletion Level",
           xlimit=c(0.2,0.35),filename = filen,devoff=FALSE)
  for (pop in 1:numpop) abline(v=xval[pop,"Deplet"],lwd=2,col=pop)
  caption <- paste0("The production curve relative to the depletion ",
                    "level of each population. Here the x-axis is ",
                    "shortened to clarify the flatness of the production ",
                    "curve about the MSY points.")
  addplot(filen,resdir=resdir,category="Production",caption)

  # Now do total production --------------------------------------------
  yield <- rowSums(product[,"Catch",])
  spb <- rowSums(product[,"MatB",])
  Ht <- yield/spb
  depletMSY <- spb/spb[1]
  pickmsy <- which.max(yield)
  maxy <- getmax(yield)

  filen <- filenametopath(resdir,"production_SpB_Total.png")
  plotprep(width=7,height=6,newdev=FALSE,filename=filen,verbose=FALSE)
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
  caption <- paste0("The production curves for the zone. Also the ",
              "relationships between spawning biomass depletion and ",
              "harvest rate.")
  addplot(filen,resdir=resdir,category="Production",caption)
} # end of plotproductivity
