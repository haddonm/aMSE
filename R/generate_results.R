
#' @title biology_plots generates a series of stored plots and tables
#'
#' @description biology_plots generates the yield vs spawning biomass,
#'     weight-at-length, and emergence plots for all populations. In
#'     addition, it also tabulates the biological properties of each
#'     population and SAU and the total zone. This is only called in
#'     do_MSE.
#'
#' @param rundir the results directory
#' @param glb the globals list
#' @param zoneC the zonal constants by population, zoneC
#' @param matL a vector of two containing the left and right hand size classes
#'     for use in the maturity-at-length plots
#' @param Lwt a vector of two containing the left and right hand size classes
#'     for use in the weight-at-length plots
#'
#' @return invisibly returns the biological properties of the populations
#' @export
#'
#' @examples
#' args(biology_plots)
biology_plots <- function(rundir, glb, zoneC, matL=c(30,210), Lwt=c(80,210)) {
  # rundir=rundir; glb=glb;zoneC=zoneC;matL=c(70,200);Lwt=c(80,200);filen=""
  mids <- glb$midpts
  numpop <- glb$numpop
  popdef <- getlistvar(zoneC,"popdef")
  SAU <- unique(popdef["SAU",])
  nSAU <- glb$nSAU
  saunames <- glb$saunames
  sauindex <- glb$sauindex
  # Yield vs Spawning biomass
  # maturation uses zoneC
  matur <- getlistvar(zoneC,"Maturity")
  rownames(matur) <- mids
  filen <- filenametopath(rundir,"maturity_v_Length_pop.png")
  plotprep(width=7,height=8,newdev=FALSE,filename=filen,cex=0.9,
           verbose=FALSE)
  parset(plots=pickbound(nSAU),margin=c(0.25,0.3,0.05,0.05),
         outmargin=c(1.5,1.5,0,0))
  for (sau in 1:nSAU) {
    pickP <- which(sauindex == sau)
    popsau <- length(pickP)
    plot(mids,matur[,pickP[1]],type="l",lwd=1,xlab="",
         ylab="",panel.first=grid(),xlim=c(matL[1],matL[2]))
    if (popsau > 1)
      for (pop in 2:popsau) lines(mids,matur[,pickP[pop]],lwd=1,col=pop)
    legend("topright",paste0("P",pickP),lwd=3,col=c(1:popsau),bty="n",
           cex=0.85)
    text(68,0.9,saunames[sau],cex=1.5,pos=4)
  }
  mtext("Shell Length (mm)",side=1,outer=TRUE,line = -0.1,cex=1.1)
  mtext("Proportion Mature",side=2,outer=TRUE,line = -0.1,cex=1.1)
  caption <- "The maturity vs length for each population in each SAU."
  addplot(filen,rundir=rundir,category="Biology",caption)
  # weight-at-length using zoneC
  WtL <- getlistvar(zoneC,"WtL")
  rownames(WtL) <- mids
  filen <- filenametopath(rundir,"Weight_at_Length_pop.png")
  plotprep(width=7,height=8,newdev=FALSE,filename=filen,cex=0.9,
           verbose=FALSE)
  parset(plots=pickbound(nSAU),margin=c(0.25,0.3,0.05,0.05),
         outmargin=c(1.5,1.5,0,0))
  ymax <- getmax(WtL,mult=1.01)
  for (sau in 1:nSAU) {
    pickP <- which(sauindex == sau)
    popsau <- length(pickP)
    plot(mids,WtL[,pickP[1]],type="l",lwd=1,xlab="",ylim=c(0,ymax),
         ylab="",panel.first=grid(),xlim=c(Lwt[1],Lwt[2]))
    if (popsau > 1)
      for (pop in 2:popsau) lines(mids,WtL[,pickP[pop]],lwd=1,col=pop)
    legend("topleft",paste0("P",pickP),lwd=3,col=c(1:popsau),bty="n",
           cex=0.85)
    labpos <- round(Lwt[1] + (Lwt[2]-Lwt[1])/2)
    text(labpos,0.85*ymax,saunames[sau],cex=1.5,pos=3)
  }
  mtext("Shell Length (mm)",side=1,outer=TRUE,line = -0.1,cex=1.1)
  mtext("Weight (g)",side=2,outer=TRUE,line = -0.1,cex=1.1)
  caption <- paste0("The weight-at-length for each population in each SAU. ",
                      "The x-axis is constrained to encompass legal sizes.")
  addplot(filen,rundir=rundir,category="Biology",caption)
  # emergence uses zoneC
  emerg <- getlistvar(zoneC,"Emergent")
  rownames(emerg) <- mids
  filen <- filenametopath(rundir,"Emergence_at_Length_pop.png")
  plotprep(width=7,height=8,newdev=FALSE,filename=filen,cex=0.9,
           verbose=FALSE)
  parset(plots=pickbound(nSAU),margin=c(0.25,0.3,0.05,0.05),
         outmargin=c(1.5,1.5,0,0))
  for (sau in 1:nSAU) {
    pickP <- which(sauindex == sau)
    popsau <- length(pickP)
    plot(mids,emerg[,pickP[1]],type="l",lwd=1,xlab="",
         ylab="",panel.first=grid(),xlim=c(matL[1],matL[2]))
    if (popsau > 1)
      for (pop in 2:popsau) lines(mids,emerg[,pickP[pop]],lwd=1,col=pop)
    legend("topright",paste0("P",pickP),lwd=3,col=c(1:popsau),bty="n",
           cex=0.85)
    text(68,0.9,saunames[sau],cex=1.5,pos=4)
  }
  mtext("Shell Length (mm)",side=1,outer=TRUE,line = -0.1,cex=1.1)
  mtext("Proportion Emergent",side=2,outer=TRUE,line = -0.1,cex=1.1)
  caption <- paste0("The emergence-at-length for each population.",
                    "The x-axis is constrained to emphasize differences.")
  addplot(filen,rundir=rundir,category="Biology",caption)
  # Tabulate biological properties uses zoneC
  rows <- c("M","R0","B0","ExB0","MSY","MSYDepl","bLML",
            "MaxDL","L50","L95","AvRec","steep")
  sau <- getlistvar(zoneC,"SAU")
  nSAU <- length(unique(sau))
  columns <- paste0("p",1:numpop)
  numrow <- length(rows)
  numcol <- length(columns)
  resultpop <- matrix(0,nrow=numrow,ncol=numpop,
                      dimnames=list(rows,columns))
  resultpop["B0",] <- as.numeric(getlistvar(zoneC,"B0"))
  resultpop["M",] <- as.numeric(getlistvar(zoneC,"Me"))
  resultpop["R0",] <- as.numeric(getlistvar(zoneC,"R0"))
  resultpop["ExB0",] <- as.numeric(getlistvar(zoneC,"ExB0"))
  resultpop["MSY",] <- as.numeric(getlistvar(zoneC,"MSY"))
  resultpop["MSYDepl",] <- as.numeric(getlistvar(zoneC,"MSYDepl"))
  resultpop["bLML",] <- as.numeric(getlistvar(zoneC,"bLML"))
  popdefs <- getlistvar(zoneC,"popdef")
  resultpop["MaxDL",] <- popdefs["DLMax",]
  resultpop["L50",] <- popdefs["L50",]
  resultpop["L95",] <- popdefs["L95",]
  resultpop["AvRec",] <- round(popdefs["AvRec",])
  resultpop["steep",] <- popdefs["steeph",]
  plotrecruitment(rundir,resultpop,glb)
  res <- as.data.frame(t(round(resultpop,3)))
  res[,"SAU"] <- getlistvar(zoneC,"SAU")
  filen <- "zonebiology.csv"
  caption <- "Population Biological Properties."
  addtable(res,filen,rundir=rundir,category="Tables",caption)
  # histograms of Me,Growth pars, MSY
  filen <- filenametopath(rundir,"Basic_Biology_pop.png")
  plotprep(width=7,height=8,newdev=FALSE,filename=filen,cex=1.0,
           verbose=FALSE)
  parset(plots=c(4,2),margin=c(0.5,0.5,0.05,0.05),cex=1.0)
  hist(resultpop["M",],xlab="M",main="",breaks=16,panel.first=grid())
  hist(resultpop["MSY",],xlab="MSY",main="",breaks=16,panel.first=grid())
  hist(resultpop["MaxDL",],xlab="MaxDL",main="",breaks=16,panel.first=grid())
  hist(resultpop["L50",],xlab="L50",main="",breaks=16,panel.first=grid())
  hist(resultpop["L95",],xlab="L95",main="",breaks=16,panel.first=grid())
  hist(resultpop["steep",],xlab="steep",main="",breaks=16,panel.first=grid())
  hist(resultpop["bLML",],xlab="bLML",main="",breaks=16,panel.first=grid())
  hist(resultpop["AvRec",]/1000,xlab="AvRec '000s",main="",breaks=20,panel.first=grid())
  caption <- paste0("The range of biological properties across the zone.")
  addplot(filen,rundir=rundir,category="Biology",caption)
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
#' @param ssc index for starting size class. thus 1 = 2, 2 = 4, 5 = 10, etc.
#'     default = 5 for it plots size classes from 10mm up
#' @param LML the legal minimum length in the comparison year
#' @param rundir the results directory, default = "" leading to no
#'     .png file, just a plot to the screen.
#'
#' @return invisibly the filename ready for logfilename
#' @export
#'
#' @examples
#' print("still to be developed")
#' # unfN=unfN; curN=depN;glb=glb; yr=1; depl=0.3993; LML=132; rundir=rundir
compzoneN <- function(unfN,curN,glb,yr,depl,ssc=5,LML=0,rundir="") {
  usecl=ssc:glb$Nclass
  mids <- glb$midpts
  filen <- paste0("zone_n-at-size_yr",yr,".png")
  filen <- filenametopath(rundir,filen)
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


#' @title fishery_plots generates a set of plots relating to the fishery properties
#'
#' @description fishery_plots produces a series of plots relating to the fishery
#'     and its properties. The selectivity curves reflect the LML through time
#'     and so this needs to be represented.
#'
#' @param rundir the rundir for the particular scenario
#' @param glb the object contianing the global constants
#' @param select the selectivity matrix from zoneCP, with all years of selectivity
#' @param histyr the object from condC containing the historical LML
#' @param projLML the object containing the LML in the projections
#' @param rge the range of the midpts for the size classes so that the view can
#'     be limited so that the lower arm set to zero and the upper set to 1 need
#'     not all be viewed and to give greater separation to any set of LML, which
#'     will usually be close to each other.
#'
#' @return invisibly, the matrix of unique selectivities
#' @export
#'
#' @examples
#' print("wait on suitable internal data sets")
fishery_plots <- function(rundir,glb,select,histyr,projLML, rge=55:105) {
  # rundir=rundir;glb=glb;select=zoneCP[[1]]$Select;histyr=condC$histyr;
  # projLML=projC$projLML; rge=55:105
  mids <- glb$midpts
  maxlml <-
  lml <- cbind(as.numeric(histyr[,"histLML"]),
               rep(tail(mids,1),length(glb$hyrnames)))
  colnames(lml) <- c("LML","Max")
  lmlall <- rbind(lml,projLML)
  allyrs <- c(glb$hyrnames,glb$pyrnames)
  rownames(lmlall) <- allyrs
  outyrs <- uniquepairs(lmlall,col1=1,col2=2)
  realcomb <- outyrs$realcomb
#  nsel <- length(realcomb)
  combnames <- outyrs$combnames[realcomb]
  ncomb <- length(combnames)
  yrlist <- outyrs$yrlist
  selyrs <- numeric(ncomb)
  yrnames <- numeric(ncomb)
  yrlml <- numeric(ncomb)
  picksel <- 0
  for (i in 1:ncomb) { # i = 1
    yrcomb <- yrlist[[combnames[i]]]
#    if (!is.null(yrcomb)) {
    selyrs[i] <- yrcomb[1]
    yrnames[i] <- allyrs[yrcomb[1]]
    yrlml[i] <- lmlall[yrcomb[1],1]
 #   }
  }
  selplot <- select[,selyrs]
  rownames(selplot) <- mids
  filen <- filenametopath(rundir,"fishery_Selectivity_through_time.png")
  plotprep(width=7,height=4,newdev=FALSE,filename=filen,cex=0.9,
           verbose=FALSE)
  parset(cex=1.0)
  plot(mids[rge],selplot[rge,1],type="l",lwd=2,col=1,panel.first=grid(),
       xlab="Shell Length (mm)",ylab="Selectivity")
  for (i in 2:ncomb) lines(mids[rge],selplot[rge,i],lwd=2,col=i)
  legend("right",legend=paste0("lml ",yrlml," ",yrnames),col=c(1:ncomb),lwd=3,
         bty="n",cex=1.2)
  caption <- paste0("The selectivity vs length through the years of the Dynamics.",
                      "The years after the LML value are the first year in the data",
                      " in which that LML is expressed.")
  addplot(filen,rundir=rundir,category="Fishery",caption)
  return(invisible(list(select=selplot,yrlml=yrlml,yrnames=yrnames)))
} # end of fishery_plots

#' @title makeoutput take the output from do_MSE and generates the HTML files
#'
#' @description makeoutput simplifies taking the output from do_MSE and
#'     producing the HTML files used to display all the results in rundir. Note
#'     that 'runnotes' is now a vector of character strings made up using paste0.
#'
#' @param out the output from do_MSE
#' @param rundir the full path to the directory holding the all scenario files
#' @param postdir the name of the directory holding the results, also used to
#'     name the internal webpage
#' @param controlfile the controlfile used to run the MSE
#' @param hsfile the name of the harvest strategy file, default=NULL
#' @param doproject have the projections been run? default = TRUE.
#' @param openfile should the website be opened automatically? default=TRUE
#' @param verbose should details of producing the HTML files be written to the
#'     console?  default=FALSE
#'
#' @return nothing but it does generate a set of HTML files in rundir. If
#'     verbose=TRUE it also writes text to the console
#' @export
#'
#' @examples
#' print("wait on internal data-sets")
#' # out=out;rundir=rundir;postdir=postfixdir;controlfile=controlfile;
#' # hsfile="TasHS Package"; doproject=FALSE;openfile=TRUE;verbose=FALSE
makeoutput <- function(out,rundir,postdir,controlfile,hsfile=NULL,
                       doproject=TRUE,openfile=TRUE,verbose=FALSE) {
  replist <- list(starttime=as.character(out$starttime),
                  endtime=as.character(out$runtime))
  glb <- out$glb
  ctrl <- out$ctrl
  projy <- ifelse(doproject,glb$pyrs,0)
  warnings <- readLines(glb$warnfile)
  numwarn <- length(warnings)
  addnotes <- NULL
  if (numwarn > 1) {
    for (i in 1:numwarn)
       addnotes <- c(addnotes,paste0("Warning ",i," ",warnings[i+1]))
  } else {
    addnotes <- c(addnotes,"No warnings reported.")
  }
  runnotes <- c(ctrl$runlabel,paste0("RunTime = ",out$tottime),
                paste0("replicates = ",glb$reps),paste0("years projected = ",projy),
                paste0("Populations = ",glb$numpop),paste0("SAU = ",glb$nSAU),
                paste0("Randomseed for conditioning = ",ctrl$randseed),
                paste0("Recruitment variability sigR         = ",ctrl$withsigR),
                paste0("Exploitable Biomass variability sigB = ",ctrl$withsigB),
                paste0("Catch-Rate variability sigCE         = ",ctrl$withsigCE),
                addnotes)

  make_html(replist = replist,  rundir = rundir,
            controlfile=controlfile, datafile=out$ctrl$datafile, hsfile=hsfile,
            width = 500, openfile = openfile,  runnotes = runnotes,
            verbose = verbose, packagename = "aMSE",  htmlname = postdir)
  if (verbose) cat("finished  \n")
} # makeoutput

#' @title numbersatsize plots details of the numbers-at-size
#'
#' @description numbersatsize plots up the initial unfished numbers-
#'     at-size distribution, omitting the first four size classes to
#'     avoid the recruitment numbers dominating the plot.
#'
#' @param rundir the results directory, if set to "" then plot is sent to
#'     the console instead
#' @param glb the globals list
#' @param zoneD the dynamic part of the zone, zoneD
#' @param ssc index for starting size class. thus 1 = 2, 2 = 4, 5 = 10, etc.
#'     default = 5 for it plots size classes from 10mm up
#'
#' @return nothing but it does add one plot to the results directory
#' @export
#'
#' @examples
#' print("this will be quite long when I get to it")
numbersatsize <- function(rundir, glb, zoneD, ssc=5) {
  # some globals
  mids <- glb$midpts
  numpop <- glb$numpop
  nc <- glb$Nclass
  # initial numbers-at-size uses zoneD--
  Nt <- zoneD$Nt[,1,]/1000.0
  Ntt <- rowSums(zoneD$Nt[,1,])/1000.0  # totals
  if (nchar(rundir) > 0) {
    filen <- file.path(rundir,"Initial_N-at-Size.png")
  } else {
    filen <- ""
  }
  plotprep(width=7,height=6,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  parset(plots=c(2,1))
  plot(mids[ssc:nc],Ntt[ssc:nc],type="l",lwd=2,xlab="Shell Length mm (5 - 210mm)",
       ylab="Numbers-at_size '000s",panel.first=grid())
  maxy <- getmax(Nt[ssc:nc,])
  plot(mids[ssc:nc],Nt[ssc:nc,1],type="l",lwd=2,xlab="Shell Length mm (5 - 210mm)",
       ylab="Numbers-at_size '000s",panel.first=grid(),ylim=c(0,maxy))
  for (pop in 2:numpop) lines(mids[ssc:nc],Nt[ssc:nc,pop],lwd=2,col=pop)
  abline(h=0.0,col="darkgrey")
  legend("topright",paste0("P",1:numpop),lwd=3,col=c(1:numpop),bty="n",
         cex=1.2)
  if (nchar(rundir) > 0) {
    caption <- paste0("The numbers-at-size for the whole zone and for each ",
                      "population separately. The recruitment numbers are",
                      " omitted for clarity.")
    addplot(filen,rundir=rundir,category="NumSize",caption)
  }
} # end of numbersatsize


#' @title plotproductivity characterizes each population's yield curve
#'
#' @description plotproductivity characterizes each population's yield curve, it
#'   also describes the total productivity of the zone.
#'
#' @param rundir the results directory
#' @param product the productivity 3-D array
#' @param glb the globals list
#' @param hsargs the maximum CPUE target input is inside hsargs. This is
#'     compared with the range of CPUE found across the produciton curve to
#'     ensure that the target is not attempting to achieve the impossible. If
#'     there is no maxtarg, this implies no comparison and no adding to the plot
#'
#' @return It produces five png files into rundir and returns sauprod, a matrix
#'     of B0, Bmsy, MSY, Depletion-msy, and CEmsy for each sau
#' @export
#'
#' @examples
#' print("this will be quite long when I get to it")
plotproductivity <- function(rundir,product,glb,hsargs) {
  # rundir=rundir; product=production;glb=glb; hsargs=hsargs
  maxtarg <- hsargs$maxtarg
  xval <- findmsy(product)
  numpop <- glb$numpop
  if (numpop <= 16) { # otherwise too much detail to see. Ok for SA not TAS
  # Yield vs Spawning biomass
    filen <- filenametopath(rundir,"production_SpB.png")
    plotprod(product,xname="MatB",xlab="Spawning Biomass t",
             ylab="Production t",filename = filen,devoff=FALSE)
    caption <- paste0("The production curve relative to each population's ",
                      "spawning biomass. The vertical lines identify the ",
                      "Bmsy values.")
    addplot(filen,rundir=rundir,category="Production",caption)

    # Yield vs Annual Harvest Rate
    filen <- filenametopath(rundir,"production_AnnH.png")
    plotprod(product,xname="AnnH",xlab="Annual Harvest Rate",filename = filen,
             devoff=FALSE)
    caption <- paste0("The production curve relative to the Annual ",
                      "Harvest Rate applied to each population. The ",
                      "vertical lines identify the Hmsy values.")
    addplot(filen,rundir=rundir,category="Production",caption)

    # plot of Yield vs population depletion
    filen <- filenametopath(rundir,"production_Deplet.png")
    plotprod(product,xname="Deplet",xlab="Population Depletion Level",
             filename = filen,devoff=FALSE)
    for (pop in 1:numpop) abline(v=xval[pop,"Deplet"],lwd=2,col=pop)
    caption <- paste0("The production curve relative to the depletion ",
                "level of each population. The vertical lines identify ",
                "the Depletion level giving rise to the MSY.")
    addplot(filen,rundir=rundir,category="Production",caption)

    # plot of Yield vs population depletion but constrained to within
    # 0.2 and 0.35 levels, to illustrate nearly flat rpoduction curve
    # and more clearly identify the population depletion at MSY
    filen <- filenametopath(rundir,"production_Deplet_0.2_0.35.png")
    plotprod(product,xname="Deplet",xlab="Population Depletion Level",
             xlimit=c(0.2,0.35),filename = filen,devoff=FALSE)
    for (pop in 1:numpop) abline(v=xval[pop,"Deplet"],lwd=2,col=pop)
    caption <- paste0("The production curve relative to the depletion ",
                      "level of each population. Here the x-axis is ",
                      "shortened to clarify the flatness of the production ",
                      "curve about the MSY points.")
    addplot(filen,rundir=rundir,category="Production",caption)
  } # end of numpop <= 16 if statement
  # what CPUE at Bmsy would be exhibited at MSY
  filen <- filenametopath(rundir,"production_CPUE_at_Bmsy.png")
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
  saucpue <- matrix(0,nrow=nh,ncol=nsau,dimnames=list(label[[1]],saunames))
  rows <- c("B0","Bmsy","MSY","Dmsy","CEmsy","Hmsy")
  sauprod <- matrix(0,nrow=length(rows),ncol=nsau,dimnames=list(rows,saunames))
  # Now do sau production
  for (i in 1:nh) {
    saumatB[i,] <- tapply(product[i,"MatB",],sauindex,sum,na.rm=TRUE)
    sauyield[i,] <- tapply(product[i,"Catch",],sauindex,sum,na.rm=TRUE)
    wts[i,] <- product[i,"Catch",]/sauyield[i,sauindex]
    saucpue[i,] <- tapply((product[i,"RelCE",] * wts[i,]),sauindex,sum,na.rm=TRUE)
  }
  label <- glb$saunames
  plotprep(width=8,height=7,newdev=FALSE,filename=filen,verbose=FALSE)
  parset(plots=pickbound(nsau),margin=c(0.3,0.3,0.05,0.05),outmargin=c(1,1,0,0))
  for (i in 1:nsau) { # i=1
    ymax <- getmax(sauyield[,i])
    if (!is.null(maxtarg)) {
      xmax <- getmax(c(saucpue[2:nh,i],maxtarg[i]))
      xmin <- getmin(c(saucpue[2:nh,i],maxtarg[i]))
    } else {
      xmax <- getmax(saucpue[2:nh,i])
      xmin <- getmin(saucpue[2:nh,i])
    }
    pick <- which.max(sauyield[,i])
    msyce <- saucpue[pick,i]
    plot(saucpue[2:nh,i],sauyield[2:nh,i],type="l",lwd=2,xlab="",ylab="",
         panel.first=grid(),ylim=c(0,ymax),yaxs="i",xlim=c(xmin,xmax))
    abline(v=msyce,col=2,lwd=2)
    if (!is.null(maxtarg)) abline(v=maxtarg[i],lwd=2,col=4)
    msylab <- paste0("MSY = ",round(sauyield[pick,i],1))
    text(0.7*max(saucpue[,i]),0.92*ymax,msylab,cex=1.2,pos=4)
    text(0.8*max(saucpue[,i]),0.75*ymax,label[i],cex=1.5,pos=4)
    text(1.1*msyce,0.15*ymax,round(msyce,2),cex=1.25,pos=4)
    sauprod[,i] <- c(saumatB[1,i],saumatB[pick,i],sauyield[pick,i],
                     (saumatB[pick,i]/saumatB[1,i]),saucpue[pick,i],
                     harv[pick])
  }
  mtext("CPUE at Bmsy (note different scales)",side=1,outer=TRUE,cex=1.0,
        line = -0.1)
  mtext("Equilibrium Yield  (note different scales)",side=2,outer=TRUE,cex=1.0,
        line=-0.2)
  caption <- paste0("The equilibrium production vs CPUE for each SAU. The red ",
                    "vertical lines and related number represent the expected",
                    "CPUE when biomass is at Bmsy. ")
  if (!is.null(maxtarg)) caption <- paste0(caption,"The blue vertical line is ",
                                           "the maxtarg for that SAU.")
  addplot(filen,rundir=rundir,category="Production",caption)
  # Now do total production
  yield <- rowSums(product[,"Catch",])
  spb <- rowSums(product[,"MatB",])
  Ht <- yield/spb
  depletMSY <- spb/spb[1]
  pickmsy <- which.max(yield)
  maxy <- getmax(yield)

  filen <- filenametopath(rundir,"production_SpB_Total.png")
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
  addplot(filen,rundir=rundir,category="Production",caption)
  # add sauprod table to production tab
  filen <- paste0("Production_by_sau.csv")
  caption <- "Productivity properties: B0, Mnsy, MSY, Dmsy, Cemsy by sau."
  addtable(sauprod,filen,rundir=rundir,category="Production",caption)
  zoneprod <- matrix(0,nrow=5,ncol=1)

  return(invisible(sauprod))
} # end of plotproductivity

#' @title plotrecruitment literally plots stock-recruitment relationships
#'
#' @description plotrecruitment plots an approximation to the summarized
#'     stock recruitment relationship for each sau together on the same plot.
#'     It also plots the stock-recruitment relationship for each population
#'     within each sau, which illustrates the distribution of the recruitment
#'     total at the sau level distributed relative to yield from the areas of
#'     persistent production results. Finally, it tabulates the sau results.
#'
#' @param rundir the scenarios results directory
#' @param resultpop a summary of the biological properties of each population
#'     generated in the biology_plots function.
#' @param glb the globals object
#'
#' @return nothing but it does save two figures and a table into rundir
#' @export
#'
#' @examples
#' print("wait on data sets")
#' # rundir=rundir; resultpop=resultpop; glb=glb
plotrecruitment <- function(rundir,resultpop,glb) {
  usevar <- c("R0","B0","ExB0","MSY")
  nvar <- length(usevar)
  allvar=c(usevar,"steep")
  numrow <- length(allvar)
  nsau <- glb$nSAU
  saunames <- glb$saunames
  sauindex <- glb$sauindex
  saures <- matrix(0,nrow=numrow,ncol=nsau,dimnames=list(allvar,saunames))
  for (i in 1:nvar) {
    saures[usevar[i],] <- tapply(resultpop[usevar[i],],sauindex,sum,na.rm=TRUE)
  }
  saures["steep",] <- tapply(resultpop["steep",],sauindex,mean,na.rm=TRUE)
  filen <- "sauproductivity_properties.csv"
  caption <- "sau Biological Properties relating to recruitment and productivity."
  addtable(round(saures,4),filen,rundir=rundir,category="Tables",caption)

  #plot each sau BH relation together
  filen <- filenametopath(rundir,"All_sau_recruitment.png")
  plotprep(width=8, height=4.5,filename=filen,verbose=FALSE)
  parset(cex=0.9)
  ymax <- getmax(saures["R0",])
  pickmax <- which.max(saures["R0",])
  steep <- saures["steep",pickmax]
  R0 <- saures["R0",pickmax]
  B0 <- saures["B0",pickmax]
  Bsp <- seq(0,B0,length=100)
  rec <- ((4*steep*R0*Bsp)/((1-steep)*B0+(5*steep-1)*Bsp))
  plot(Bsp,rec,type="l",lwd=2,col=1,xlab="Spawning Biomass t",
       ylab="Recruitment", panel.first=grid(),ylim=c(0,ymax))
  for (sau in 1:nsau) {
    steep <- saures["steep",sau]
    R0 <- saures["R0",sau]
    B0 <- saures["B0",sau]
    rec <- ((4*steep*R0*Bsp)/((1-steep)*B0+(5*steep-1)*Bsp))
    lines(Bsp,rec,lwd=2,col=sau)
  }
  legend("topleft",legend=saunames,col=1:nsau,lwd=3,bty="n")
  caption <- paste0("The Stock-Recruitment relationships for each sau.")
  addplot(filen,rundir=rundir,category="Recruits",caption)

  filen <- filenametopath(rundir,"Each_Populations_recruitment.png")
  plotprep(width=7, height=9,filename=filen,verbose=FALSE)
  parset(plots=pickbound(nsau),margin=c(0.3,0.4,0.05,0.1),
         outmargin=c(1,1,0,0), cex=0.9,byrow=FALSE)
  for (sau in 1:nsau) { # sau = 1
    pickcol <- which(sauindex == sau)
    nextra <- length(pickcol)
    tmp <- as.matrix(resultpop[,pickcol])
    ymax <- getmax(tmp["R0",])
    pickmax <- which.max(tmp["R0",])
    steep <- tmp["steep",pickmax]
    R0 <- tmp["R0",pickmax]
    B0 <- tmp["B0",pickmax]
    Bsp <- seq(0,B0,length=100)
    rec <- ((4*steep*R0*Bsp)/((1-steep)*B0+(5*steep-1)*Bsp))
    plot(Bsp,rec,type="l",lwd=1,col=0,xlab="",ylab=saunames[sau],
         panel.first=grid(),ylim=c(0,ymax))
    for (i in 1:nextra) {
      steep <- tmp["steep",i]
      R0 <- tmp["R0",i]
      B0 <- tmp["B0",i]
      Bsp <- seq(0,B0,length=100)
      rec <- ((4*steep*R0*Bsp)/((1-steep)*B0+(5*steep-1)*Bsp))
      lines(Bsp,rec,lwd=2,col=i)
    }
  }
  mtext("Spawning Biomass t",side=1,outer=TRUE,cex=1.1,
        line=-0.2)
  mtext("Recruitment by Population",side=2,outer=TRUE,cex=1.1,
        line=-0.1)
  caption <- paste0("The Stock-Recruitment relationships for each ",
                    "population in each sau.")
  addplot(filen,rundir=rundir,category="Recruits",caption)
} # end of plotrecruitment
