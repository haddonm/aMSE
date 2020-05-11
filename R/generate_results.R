  # These plots use both regionC and regionD
  # characterize productivity ----------------------------------------


#' @title biology_plots generates a series of stored plots and tables
#'
#' @description biology_plots generates the yield vs spawning biomass,
#'     weight-at-lenght, and emergence plots for all populations. In
#'     addition, it also tabulates the biological properties of each
#'     population and SMU
#'
#' @param resdir the results directory
#' @param runlabel the runlabel, ideally from the ctrl object
#' @param glb the globals list
#' @param regC the regional constants by population, regionC
#' @param regD the dynamic part of the region, regionD
#' @param product the productivity 3-D array
#'
#' @return nothing but it does plot 3 plots and put one table into
#'     resdir
#' @export
#'
#' @examples
#' print("this will be quite long when I get to it")
biology_plots <- function(resdir, runlabel, glb, regC, regD, product) {
  #      resdir=resdir; regC=regionC; regD=regionD; product=product; glb=glb; runlabel=ctrl$runlabel
  # some globals
  mids <- glb$midpts
  numpop <- glb$numpop
  resfile <- paste0(resdir,"resultTable_",runlabel,".csv")
  # Yield vs Spawning biomass-----------------------
  # maturation uses regC
  matur <- getlistvar(regC,"Maturity")
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

  # weight-at-length using regC--------------------------------------
  WtL <- getlistvar(regC,"WtL")
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

  # emergence uses regC----------------------------------------------
  emerg <- getlistvar(regC,"Emergent")
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


  # Tabulate biological properties uses regC-------------------------
  getvar <- function(indexvar,reg=regC) { # indexvar="B0"
    invar <- getlistvar(reg,indexvar)
    ans <- c(invar,sum(invar))
    return(ans)
  }
  rows <- c("SMU","M","R0","B0","effB0","ExB0","effExB0","MSY","MSYDepl","bLML")
  columns <- c(paste0("p",1:numpop),"region")
  numrow <- length(rows)
  numcol <- length(columns)
  results <- matrix(0,nrow=numrow,ncol=numcol,
                                  dimnames=list(rows,columns))
  results["SMU",] <- c(getlistvar(regC,"SMUname"),0) # no total
  results["B0",] <- as.numeric(getvar("B0"))
  M <- getlistvar(regC,"Me")
  wtr <- (results["B0",1:numpop]/results["B0",(numpop+1)])
  results["M",] <- c(M,sum(M*wtr))
  results["R0",] <- getvar("R0")
  results["effB0",] <- getvar("effB0")
  results["ExB0",] <- getvar("ExB0")
  results["effExB0",] <- getvar("effExB0")
  results["MSY",] <- getvar("MSY")
  MSYD <- getlistvar(regC,"MSYDepl")
  results["MSYDepl",] <- c(MSYD,sum(MSYD * wtr))
  bLML <- getlistvar(regC,"bLML")
  results["bLML",] <- c(bLML,sum(bLML * wtr))
  res <- round(results,3)
  filen <- paste0("regionbiology_",runlabel,".csv")
  filename <- filenametopath(resdir,filen)
  write.table(res,file = filename,sep=",")
  caption <- paste0("Population amd Regional Biological Properties. ",
                   "For 'M' 'MSYDepl' and 'bLML' total is an average ",
                   "weighted relative to the proportion of total B0.")
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
#' @param regD the dynamic part of the region, regionD
#'
#' @return nothing but it does add one plot to the results directory
#' @export
#'
#' @examples
#' print("this will be quite long when I get to it")
numbersatsize <- function(resdir, runlabel, glb, regD) {
  #      resdir=resdir; regC=regionC; regD=regionD; product=product; glb=glb; runlabel=ctrl$runlabel
  # some globals
  mids <- glb$midpts
  numpop <- glb$numpop
  resfile <- paste0(resdir,"resultTable_",runlabel,".csv")
  # initial numbers-at-size uses regD--------------------------------
  Nt <- regD$Nt[,1,]/1000.0
  Ntt <- rowSums(regD$Nt[,1,])/1000.0  # totals
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
    caption <- "The numbers-at-size for the whole region and for each population separately. The recruitment numbers are omitted for clarity."
  logfilename(filename,resfile=resfile,"NumSize",caption)



} # end of numbersatsize


#' @title plotproductivity characterizes each population's yield curve
#'
#' @description plotproductivity characterizes each population's yield
#'     curve, it also describes the total productivity of the region.
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
  resfile <- paste0(resdir,"resultTable_",runlabel,".csv")
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
  caption <- paste0("The production curves for the region. Also the ",
              "relationships between spawning biomass depletion and ",
              "harvest rate.")
  logfilename(filename,resfile=resfile,"Production",caption)


} # end of plotproductivity
