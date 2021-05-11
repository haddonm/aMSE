

#' @title initiateHS applies HS to the last years of the historical conditioning
#'
#' @description initiateHS uses the last years of fishery conditioning to run
#'     the HS for the first year of the projections. This fills in the first
#'     row of each of the variables in the dynamic object for each of the
#'     proposed iterations.
#'
#' @param zoneDP the empty dynamics object used in the projections
#' @param zoneCP the constant object after modification for projections
#' @param exb a vector of exploitable biomass for each population from the last
#'     year of the conditioning on the fishery or historical fishery data
#' @param inN an array of the numbers-at-length for each population from the
#'     last year of the conditioning on the fishery or historical fishery data
#' @param acatch the aspirational catches per SAU derived from the conditioning
#'     of the HS based on the historical or fishery conditioning data
#' @param sigmar the recruitment variability (sd) to be used in the projections
#' @param sigmab the variability (sd) applied to the exploitable biomass
#'     estimates to introduce a form of implementation error on the distribution
#'     of catches among SAU
#' @param glb the globals object
#'
#' @return a revised zoneDP object with the first year filled in for each
#'     replicate
#' @export
#'
#' @examples
#' print("wait on suitable internal data sets")
initiateHS <- function(zoneDP,zoneCP,inN,acatch,sigmar,sigmab,glb) {
  reps <- dim(zoneDP$matureB)[3]
  r0 <- getvar(zoneCP,"R0") #sapply(zoneC,"[[","R0")
  b0 <- getvar(zoneCP,"B0") #sapply(zoneC,"[[","B0")
  exb0 <- getvar(zoneCP,"ExB0")
  for (iter in 1:reps) {
    popC <- calcexpectpopC(TAC=0,acatch=acatch,
                           exb=zoneDD$exploitB[year-1,],
                           sauindex,sigmab=sigmab)
    outy <- oneyearsauC(zoneCC=zoneCP,inN=inN,popC=popC,year=1,
                        Ncl=glb$Nclass,sauindex=glb$sauindex,movem=glb$move,
                        sigmar=sigmar)
    dyn <- outy$dyn
    saudyn <- poptosauCE(dyn["catch",],dyn["cpue",],glb$sauindex)
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
} # end on initiateHS


#' @title oneyearC conducts one year's dynamics using catch not harvest
#'
#' @description oneyearC conducts one year's dynamics in the simulation
#'     using catches rather than harvest rates. The harvest rates are
#'     estimated after first estimating the exploitable biomass.
#'     returning the revised zoneD, which will have had a single year
#'     of activity included in each of its components.
#'
#' @param zoneC the constant portion of the zone with a list of
#'     properties for each population
#' @param zoneDP the dynamics portion of the zone, with matrices and
#'     arrays for the dynamic variables of the dynamics of the
#'     operating model
#' @param catchp a vector of catches to be taken in the year from each
#'     population
#' @param year the year of the dynamics, would start in year 2 as year
#'     1 is the year of initiation.
#' @param iter the specific replicate being considered
#' @param sigmar the variation in recruitment dynamics, set to 1e-08
#'     when searching for an equilibria.
#' @param Ncl the number of size classes used to describe size, global Nclass
#' @param npop the number of populations, the global numpop
#' @param movem the larval dispersal movement matrix, global move
#'
#' @return a list containing a revised dynamics list
#' @export
#'
#' @examples
#' print("Wait on new data")
#' # data(zone)
#' # zoneC <- zone$zoneC
#' #  glb <- zone$glb
#'  #zoneD <- zone$zoneD
#'  #Nc <- glb$Nclass
#'  #hyrs <- glb$hyrs
#'  #catch <- 900.0 # larger than total MSY ~ 870t
#'  #B0 <- getvar(zoneC,"B0")
#'  #totB0 <- sum(B0)
#'  #prop <- B0/totB0
#'  #catchpop <- catch * prop
#'  #for (yr in 2:hyrs)
#'  #  zoneD <- oneyearC(zoneC=zoneC,zoneD=zoneD,Ncl=Nc,
#'  #                  catchp=catchpop,year=yr,sigmar=1e-08,
#'  #                  npop=glb$numpop,movem=glb$move)
#'  #str(zoneD)
#'  #round(zoneD$catchN[60:105,1:5,1],1)
oneyearC <- function(zoneC,zoneDP,catchp,year,iter,sigmar,Ncl,npop,movem) {
  matb <- numeric(npop)
  for (popn in 1:npop) {  # year=2
    out <- oneyearcat(inpopC=zoneC[[popn]],inNt=zoneDP$Nt[,year-1,popn,iter],
                      Nclass=Ncl,incat=catchp[popn],yr=year)
    zoneDP$exploitB[year,popn,iter] <- out$ExploitB
    zoneD$matureB[year,popn] <- out$MatureB
    zoneD$catch[year,popn] <- out$Catch
    zoneD$harvestR[year,popn] <- out$Harvest
    zoneD$cpue[year,popn] <- out$ce
    zoneD$Nt[,year,popn] <- out$Nt
    zoneD$catchN[,year,popn] <- out$CatchN
    matb[popn] <- out$MatureB
  }
  steep <- getvect(zoneC,"steeph") #sapply(zoneC,"[[","popdef")["steeph",]
  r0 <- getvar(zoneC,"R0") #sapply(zoneC,"[[","R0")
  b0 <- getvar(zoneC,"B0") #sapply(zoneC,"[[","B0")
  recs <- oneyearrec(steep,r0,b0,matb,sigR=sigmar)
  newrecs <- movem %*% recs
  zoneD$recruit[year,] <- newrecs
  zoneD$Nt[1,year,] <- newrecs
  zoneD$deplsB[year,] <- zoneD$matureB[year,]/b0
  zoneD$depleB[year,] <- zoneD$exploitB[year,]/getvar(zoneC,"ExB0")
  return(zoneD)
} # end of oneyearC



#' @title product is the productivity curve matrix from doproduction
#'
#' @description product is the productivity curve matrix from
#'     doproduction when the example zone is generated using the
#'     inbuilt datasets ctrl, zone1, and constants. The slowest
#'     part of building the whole is to use the modregC function
#'     to adjust the zoneC and generate the production array. To
#'     save that time in the examples (to avoid time limits on
#'     examples should this package go to CRAN), then this dataset can
#'     be used instead. This is a three dimensional array of
#'     productivity variables.
#'
#' @name product
#'
#' @docType data
#'
#' @section contents:
#' \itemize{
#'   \item harvestrate the initial harvest rates applied
#'   \item productivity variables ExB, MatB, AnnH, Catch, Deplet, RelCE
#'   \item population the index of each population
#' }
#'
#' @examples
#'  data(product)
#'  product[1:20,,1]
NULL

#' @title testzoneC is a zone list made up of 6 equilibrium populations
#'
#' @description testzoneC is a zone list made up of 6 equilibrium
#'     populations. These have been run with a laral dispersal rate of
#'     0.03 so the change from B0 to effB0 is not great, but still
#'     required for an initial equilibrium. This is here to simplify
#'     the internal testing of funcitons that require a completed
#'     zone starting at equilibrium. Its name is to avoid conflict
#'     with any actual use of zoneC. use str(testzoneC, max.level=1)
#'     to see its format. It can be expected to be used with testzoneD
#'
#' @name testzoneC
#'
#' @docType data
#'
#' @section Subjects:
#'  \itemize{
#'    \item testing of functions that require a full zone
#'    \item initial equilibrium
#'  }
#'  @export
#'
#' @examples
#'  data(testzoneC)
#'  data(testzoneD)
#'  data(zone1)
#'  glb <- zone1$globals
#'  r0 <- getvar(testzoneC,"R0")
#'  move <- makemove(glb$numpop,r0,glb$larvdisp)
#'  glb$move <- move
#'  ans <- testequil(testzoneC, testzoneD, glb)
#'  str(testzoneC[[1]])
NULL

#' @title testzoneD is a list of 8 matrices and 2 arrays defining the dynamics of a zone
#'
#' @description testzoneD is a list of 8 matrices and 2 arrays defining
#'     the dynamics of a zone. These have been run with a larval
#'     dispersal rate of 0.03 to achieve an initial equilibrium. This
#'     is here to simplify the internal testing of functions that
#'     require a completed zone starting at equilibrium. Its name is
#'     to avoid conflict with any actual use of zoneD. use
#'     str(testzoneD, max.level=1) to see its format. It can be
#'     expected to be used with testzoneC.
#'
#' @name testzoneD
#'
#' @docType data
#'
#' @section Subjects:
#'  \itemize{
#'    \item testing of functions that require a full zone
#'    \item initial equilibrium
#'  }
#'  @export
#'
#' @examples
#'  data(testzoneC)
#'  data(testzoneD)
#'  data(zone1)
#'  glb <- zone1$globals
#'  r0 <- getvar(testzoneC,"R0")
#'  move <- makemove(glb$numpop,r0,glb$larvdisp)
#'  glb$move <- move
#'  ans <- testequil(testzoneC, testzoneD, glb)
#'  str(testzoneD)
NULL


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
#' @param rundir the results directory used by makehtml to store plots and
#'     tables. If set to "", the default, then it plots to the console
#' @param defpar define the plot parameters. Set to FALSE if using
#'     plothistcatch to add a plot to a multiple plot
#'
#' @return nothing, it adds a plot into rundir and modifies resultTable.csv
#' @export
#'
#' @examples
#' print("wait until I have altered the internals data sets")
plothistcatch <- function(zone1,pops=NULL,rundir="",defpar=TRUE) {
  # zone1=zone1; pops=c(1,2,3,8); defpar=FALSE;  rundir=rundir
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
  if (nchar(rundir) > 0) {
    filen <- filenametopath(rundir,pngfile)
  } else { filen <- "" }
  if (defpar)
    plotprep(width=7,height=4,newdev=FALSE,filename=filen,cex=0.9,verbose=FALSE)
  ymax <- getmax(histcatch)
  plot(yearC,histcatch[,1],type="l",lwd=2,xlab="Year",
       ylab="Catch (t)",panel.first=grid(),ylim=c(0,ymax))
  if (numpop > 1) for (pop in 2:numpop)
    lines(yearC,histcatch[,pop],lwd=2,col=pop)
  legend("topright",label,lwd=3,col=c(1:numpop),bty="n",cex=1.2)
  if (nchar(rundir) > 0) {
    caption <- "The historical catch used for conditioning for each population."
    addplot(filen,rundir=rundir,category="history",caption)
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
#' @param rundir the results directory used by makehtml to store plots and
#'     tables.
#'
#' @return nothing, it adds a plot into rundir and modifies resultTable.csv
#' @export
#'
#' @examples
#' print("wait until I have altered the internals data sets")
plothistCE <- function(zone1,pops=NULL,rundir="") {
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
  filen <- filenametopath(rundir,pngfile)
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
  addplot(filen,rundir=rundir,category="history",caption)
} # end of plothistCE



#' @title plotzoneproj plots the replicates of a single zone-wide variable
#'
#' @description plotzoneproj takes the results from the function aszone, which
#'     contains a set of variables zonesB, zoneeB, zoneC, zoneH, zoneR, zonece,
#'     zonedeplsB,and zonedepleB, and plots whichever is selected. It contains
#'     the option of also plotting the median and the inner 90 percent quantiles
#'     of each year's spread of values across replicates. DEPRECATED
#'
#' @param zoneV the variable from within the output of aszone
#' @param reps the number of replicates to plot. Generaly one would plot all of
#'     them. If not then it might be best to turn addqnts to FALSE
#' @param yrs the vector of years as in 1:(inityrs + projyrs)
#' @param label the y-axis label, which should obviously reflect which variable
#'     is chosen from the zone summary list.
#' @param addqnts should the median and inner 90 percent quantiles be added to
#'     the plot; default = TRUE
#' @param miny sets the lower limit of the y-axis, default=0
#'
#' @return if addqnts=TRUE the quantiles are returned invisibly
#' @export
#'
#' @examples
#' print("wait on more time")
#' # zoneV=zoneproj$zoneC;reps=reps;yrs=1:50;miny=0;label="Catches"; addqnts=TRUE
plotzoneproj <- function(zoneV,reps,yrs,label="",addqnts=TRUE,miny=0) {
  maxy <- getmax(zoneV)
  ylabel <- "Variable"
  if (nchar(label) > 0) ylabel <- label
  plot(yrs,zoneV[,1],type="n",panel.first=grid(),ylim=c(miny,maxy),yaxs="i",
       ylab=ylabel,xlab="Years")
  for (iter in 1:reps) lines(yrs,zoneV[,iter],lwd=1,col="grey")
  if (addqnts) {
    CI <- apply(zoneV,1,quantile,probs=c(0.05,0.5,0.95))
    lines(yrs,CI[1,],lwd=2,col=4)
    lines(yrs,CI[2,],lwd=2,col=2)
    lines(yrs,CI[3,],lwd=2,col=4)
  }
  if (addqnts) return(invisible(t(CI)))
} # end of plotzoneproj

#' @title plotproj aids the plotting of a projected variable
#'
#' @description plotproj plots out the projections from the MSE for a selected
#'     variable. iters is included so that details of a few trajectories can
#'     be viewed as well as seeing all replicates. DEPRECATED
#'
#' @param invar the array containing the year x SAU x reps values for a selected
#'     variable out of sauzoneDP
#' @param varlabel the label to place along the outer Y-axis
#' @param plotconst a list containing nsau, saunames, reps, projyrs, and plts
#'     a vector of two describing the layout of plots
#' @param vline default=NULL, which means nothing extra is plotted. If given a
#'     year then a vertical line in red will be added to the plot
#' @param iters default=0, which means all iterations will be plotted. If iters
#'     has a value then only that many trajectories will be plotted
#' @param addqnts will calculate and add the median and 90 percent quantiles
#' @param miny sets the lower limit of the y-axis, default=0
#'
#' @return nothing but it does generate a plot
#' @export
#'
#' @examples
#' print("wait on data files")
plotproj <- function(invar,varlabel,plotconst,miny=0,
                     vline=NULL,iters=0,addqnts=FALSE) {
  nsau <- plotconst$nsau
  yrs <- 1:plotconst$projyrs
  saunames <- plotconst$saunames
  reps <- plotconst$reps
  parset(plots=plotconst$plts,byrow=FALSE,margin=c(0.25,0.4,0.1,0.05),
         outmargin=c(1,1,0,0))
  for (sau in 1:nsau) { # sau=1
    maxy <- getmax(invar[,sau,])
    plot(yrs,invar[,sau,1],lwd=1,type="l",col="grey",panel.first=grid(),
         ylim=c(miny,maxy),xlab="",ylab=saunames[sau])
    trajs <- reps
    if (iters > 0) trajs <- iters
    for (iter in 1:trajs) lines(yrs,invar[,sau,iter],lwd=1,col="grey")
    if (!is.null(vline)) abline(v=vline,lwd=1,col=2)
    if (addqnts) {
      CI <- apply(invar[,sau,],1,quantile,probs=c(0.05,0.5,0.95))
      lines(yrs,CI[1,],lwd=1,col=4)
      lines(yrs,CI[2,],lwd=2,col=4)
      lines(yrs,CI[3,],lwd=1,col=4)
    }
  }
  mtext("Years",side=1,line=-0.1,outer=TRUE,cex=1.0,font=7)
  mtext(varlabel,side=2,line=-0.1,outer=TRUE,cex=1.0,font=7)
} # end of plotproj
