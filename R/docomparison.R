

#' @title do_comparison is a wrapper function that compares scenarios
#'
#' @description do_comparison provides a simplified interface for when making
#'     comparisons between multiple scenarios. It uses the saved .RData files
#'     from each scenario and it is up to the user to put those filenames into
#'     a vector of characters. The juris argument is an attempt at generalizing
#'     this comparison function and will allow for different jurisdictions to
#'     enable future additions to do_comparison that will suit their specific
#'     requirements. Currently, the only option that has been implemented is
#'     for juris='TAS', but others can be added now that the basic framework
#'     is implemented.
#'
#' @param rundir the complete path to the directory into which all results and
#'     plots from the comparisons will be placed.
#' @param postfixdir the name of the final sub-directory into which the results
#'     will be placed. This is used as the name of the website generated to
#'     display the results
#' @param outdir the full path to the directory holding all the required .RData
#'     files. If set to '' then files is assumed to include the full path
#' @param files the vector of RData file names to be used. If outdir set to ''
#'     then files is assumed to include the full path as well as the file name
#' @param pickfiles a vector of indices selecting the files to be used.
#' @param verbose should progress updates be made to the console, default=TRUE
#' @param intensity is the density of the rgb colours used in the ribbon plots.
#'     the default = 100
#' @param zero should the phase plots have an origin at zero, default=FALSE
#' @param Q90 default = TRUE which means the inner 90th quantiles will be
#'     plotted. If FALSE then the 95th quantiles will be used.
#' @param altscenes a vector of alternative scenario names for use in plots and
#'     tables where scenario names are very long and confuse results. default=
#'     NULL, which implies using the base scenario names (runlabels).
#' @param juris which jurisdiction is being examined? default="". This is
#'     an argument so that if there are jurisdiction specific functions for
#'     examining outputs that may be of particular interest to a particular
#'     jurisdiction a suitable function can be written and input, along with
#'     any extra arguments it might have, in the final ellipsis.
#' @param ribbonleg location of the legend in the ribon plots, default=topleft
#' @param scencol what sequence of colours should each scenario have. The default=NULL
#'     which implies the colour will reflect the sequence in which they are
#'     read in.
#' @param ... the ellipsis is here in case a jurisdiction specific function or
#'     functions is/are written to perform extra analyses, plots, and tables.
#'
#' @seealso{
#'    \link{scenebyvar}, \link{scenebyzone}. \link{RGB},
#'    \link{getfilestocompare}
#' }
#'
#' @return nothing but it does conduct a comparison of at least two scenarios
#'     and places tables and plots into a given sub-directory
#' @export
#'
#' @examples
#' \dontrun{
#'   suppressPackageStartupMessages({
#'   library(aMSE)
#'   library(TasHS)
#'   library(codeutils)
#'   library(hplot)
#'   library(makehtml)
#'   library(knitr)
#'   })
#'   dropdir <- getDBdir()
#'   prefixdir <- paste0(dropdir,"A_codeUse/aMSEUse/scenarios/")
#'   postfixdir <- "BC_compare"
#'   rundir <- filenametopath(prefixdir,postfixdir)
#'   # normally one would use code to select the files, outdir is used if all
#'   # RData files are in one directory, otherwise set outdir='' and files
#'   # should contain the full paths as well as the filenames
#'   files=c("BC.RData","BCmeta.RData","BC541.RData")
#'   do_comparison(rundir,postfixdir,outdir,files,pickfiles=c(1,2,3))
#' }
do_comparison <- function(rundir,postfixdir,outdir,files,pickfiles,verbose=TRUE,
                          intensity=100,zero=FALSE,Q90=TRUE,altscenes=NULL,
                          juris="",ribbonleg="topleft",scencol=NULL,...) {
  # rundir=rundir;postfixdir=postfixdir;outdir=outdir;files=files;pickfiles=c(1,2)
  #  verbose=TRUE; intensity=100; zero=FALSE; altscenes=NULL
  #  juris="";ribbonleg="topleft"; Q90=TRUE; scencol=c(1,2)
  # get files -------------------
  files2 <- files[pickfiles]
  nfile <- length(pickfiles)
  label <- vector(mode="character",length=nfile)
  for (i in 1:nfile) label[i] <- unlist(strsplit(files2[i],".",fixed=TRUE))[1]
  ans <- makelist(label) # vector(mode="list",length=nfile)
  dyn <- makelist(label)
  glbc <- makelist(label)
  ctrlc <- makelist(label)
  condCc <- makelist(label)
  prods <- makelist(label)
  scenes <- vector(mode="character",length=nfile)
  scores <- makelist(label)
  zone <- makelist(label)
  for (i in 1:nfile) { # i = 2
    if (nchar(outdir) == 0) {
      filename <- files2[i]
    } else {
      filename <- pathtopath(outdir,files2[i])
    }
    if (verbose)
      cat("Loading ",files2[i]," which may take time, be patient  \n")
    out <- NULL # so the function knows an 'out' exists
    load(filename)
    ans[[i]] <- out
    dyn[[i]] <- out$sauout
    glbc[[i]] <- out$glb
    ctrlc[[i]] <- out$ctrl
    condCc[[i]] <- out$condC
    prods[[i]] <- t(out$sauprod$sauprod)
    scenes[i] <- out$ctrl$runlabel
    scores[[i]] <- out$outhcr
    zone[[i]] <- out$outzone
  }
  cat("\n")
  # end results extraction----------------------
  setuphtml(rundir=rundir)
  scenarionames <- scenes
  nscenes <- length(scenes)
  if (!is.null(altscenes)) {
    if (length(altscenes) != nscenes) {
      label <- paste0("The number of altscenes is different to the number ",
                      "of scenarios \n")
      stop(cat(label))
    }
    scenes <- altscenes
  }
  test <- numeric(nscenes - 1)
  for (i in 2:nscenes) test[i] <- abs(glbc[[i-1]]$pyrs - glbc[[i]]$pyrs)
  if (!all(test < 1)) {
    label <- paste0("Number of projection years differs among scenario. \n\n")
    cat(label)
  }
  # end consistency checks--------------------------
  # extract catches and CPUE-------------------------------
  catch <- makelist(scenes)
  cpue <- makelist(scenes)
  for (i in 1:nscenes) {
    catch[[i]] <- dyn[[i]]$catch
    cpue[[i]] <- dyn[[i]]$cpue
  }
  # scenario properties-----------------------------------
  scenprops <- scenarioproperties(scenes,glbc,ctrlc,condCc)
  if ((verbose) & (any(scenprops[,"same"] == 0))) {
    warnlab <- "At least one scenario property differs betweem scenarios"
    cat(warnlab," \n\n")
    addtext(warnlab,rundir=rundir,filename="warnings.txt",
                      category="scenes")
  }
  if (tolower(juris) == "tas") {
     cat("A specific jurisdiction fnction is required \n")
  }
  if (verbose) print(scenprops)
  if (verbose) cat("Now doing the comparisons  \n")
  filename <- "scenarioproperties.csv"
  addtable(scenprops,filen=filename,rundir=rundir,category="scenes",
           caption=paste0("Important scenario properties for comparability."))
  # dynamics tab-----------------------------------------
  quantscen <- comparedynamics(rundir=rundir,dyn,glbc,scenes)
  # productivity tab---------------------------------
  tabulateproductivity(rundir,prods,scenes)
  # Final Scores tab---------------------------------
  scorenames <- sapply(scores,length)
  allpos <-  (all(scorenames > 0))
  if (is.na(allpos)) allpos <- FALSE
  if ((allpos) & (max(scorenames) == min(scorenames))) {
    scenscore <- lapply(scores,"[[","finalsc")
    if (!is.null(scenscore[[1]])) {
      lengths <- numeric(nscenes)
      for (i in 1:nscenes)
        lengths[i] <- length(apply(scenscore[[i]][,1,],1,median))
      maxlen <- max(lengths)
      doplot <- ifelse(all(lengths == maxlen),TRUE,FALSE)
    } else {
      doplot <-  FALSE
    }
    if (doplot) {
      if (!is.null(scenscore[[1]])) {
        filename <- "compare_final_HSscores.png"
        meds <- comparefinalscores(rundir,scores,scenes,legloc="bottomright",
                                   filen=filename,category="Scores")
        addplot(filen=filename,rundir=rundir,category="Scores",
                caption="The HS final scores for each sau.")
        tabulatefinalHSscores(rundir,meds,scenes,category="Scores")
      } else {
        warning(cat("No Final Score or Different years of projection,",
                    " no Score tab produced \n"))
      }
    }
  } # end of if finalscore if statement
  # HSPM tab--------------------------------------------
  if (verbose) cat("Now doing the HSPM tab \n")
  outcatchHSPM <- calccatchHSPM(catch,glbc,scenes,aavcyrs=10)
  medHSPM <- outcatchHSPM$medians
  label <- names(medHSPM)
  for (i in 1:length(medHSPM)) {
    filename <- paste0("proj",label[i],".csv")
    addtable(medHSPM[[label[i]]],filen=filename,rundir=rundir,category="HSPM",
             caption=paste0("Median ",label[i]," values by sau and scenario."))
  }
  # catches tab-------------------------------------
  outtab <- catchHSPM(rundir,hspm=outcatchHSPM,glbc,scenes,
                      filen="compare_catches_boxplots.png",aavcyrs=10)
  # catchBoxPlots tab------------------------------
  boxbysau(rundir,hspm=outcatchHSPM,glbc,scenes,compvar="aavc",
           filen="sau_aavc_boxplots.png",aavcyrs=10,maxval=0)
  boxbysau(rundir,hspm=outcatchHSPM,glbc,scenes,compvar="sum5",
           filen="sau_sum5_boxplots.png",aavcyrs=10,maxval=0) #
  boxbysau(rundir,hspm=outcatchHSPM,glbc,scenes,compvar="sum10",
           filen="sau_sum10_boxplots.png",aavcyrs=10,maxval=0) #
  catchbpstats(rundir,outtab)
  # cpueBoxPlots tab-------------------------------------------
  cpueHSPM(rundir,cpue,glbc,scenes=scenes,filen="sau_maxce_boxplots.png")
  outcpue <- cpueboxbysau(rundir,cpue,glbc,scenes,filen="sau_maxceyr_box.png",
                          startyr=0,maxval=0)
  # C_vs_MSY tab-----------------------------------------
  if (verbose) cat("Now doing the C_vs_MSY tab  \n")
  outzonecatch <- zonecatchvsmsy(prods,zone,glbc)  # first at zone scale
  zcvsmsy <- outzonecatch$zcvsmsy
  zoneribbon(rundir=rundir,scenes=names(zcvsmsy),invar=zcvsmsy,glbc=glbc,
             varname="Zone-Catch-div-MSY",category="C_vs_MSY",console=FALSE,
             q90=TRUE,intens=127,addleg="topleft",addmedian=0,add1line=TRUE)
  zyrtomsy <- outzonecatch$zyrtomsy
  plotzoneyrtovar(rundir=rundir,invar=zyrtomsy,varname="MSY",glbc=glbc,
                  category="C_vs_MSY",console=FALSE)
  # now at sau scale
  cvMSYout <- catchvsMSY(catch,glbc,prods,scenes)
  cdivmsy <- cvMSYout$cdivmsy
  nscen <- length(scenes)
  for (i in 1:nscen) {
    plotsceneproj(rundir,cdivmsy[[i]],glbc[[i]],scenes[i],
                  filen=paste0("Catch_div_MSY_",scenes[i],".png"),
                  label="Catch / MSY",hline=1,Q=90)
  }
  yrtomsy <- cvMSYout$yrtomsy

  # ScenarioPMs tab-------------------------------------------------
  if (verbose) cat("Now doing the Scenario PMs  \n")
  nscen <- length(scenes)
  glb <- glbc[[1]]
  saunames <- glb$saunames
  nsau <- glb$nSAU
  x <- makelist(scenes)
  whichyrs <- (glb$hyrs + 1):(glb$hyrs + glb$pyrs)
  for (scen in 1:length(scenes)) x[[scen]] <- dyn[[scen]]$catch[whichyrs,,]
  filename <- pathtopath(rundir,"catchprojection_deviates.png")
  devout <- plotdevs(x,scenes,saunames,filen=filename)
  addplot(filen=filename,rundir=rundir,category="ScenarioPMs",
          caption="The deviates from a loess on catches in each replicate.")
  tabfile <- "projmeandeviates_catch.csv"
  addtable(devout$meandevs,filen=tabfile,rundir=rundir,category="ScenarioPMs",
           caption="Mean deviations from a loess on catch values by sau and scenario.")
  tabfile <- "proj_sddeviates_catch.csv"
  addtable(devout$sddevs,filen=tabfile,rundir=rundir,category="ScenarioPMs",
           caption="StDev of deviations from a loess on catch values by sau and scenario.")
  x <- makelist(scenes)
  whichyrs <- (glb$hyrs + 1):(glb$hyrs + glb$pyrs)
  for (scen in 1:length(scenes)) x[[scen]] <- zone[[scen]]$catch[whichyrs,]
  filename <- pathtopath(rundir,"catch-projdevs-by-zone.png")
  zdevout <- plotzonedevs(x,scenes,glb,filen=filename)
  addplot(filen=filename,rundir=rundir,category="ScenarioPMs",
          caption="The deviates from a loess on catches in each replicate by zone by scenario.")
  pickyrs <- c(glb$hyrs+5,glb$hyrs+10) # 5 and 10 years
  nyrs <- length(pickyrs)
  res <- state_by_year(indyn=dyn,whichvar="catch",whichyrs=pickyrs,
                       whichfun=mean,glb=glb)
  for (whichyr in 1:nyrs) {
    filename <- plotdynvarinyear(rundir=rundir,dyn=dyn,whichvar="catch",
                                 whichyr=pickyrs[whichyr],glb,console=FALSE)
    addplot(filen=filename,rundir=rundir,category="ScenarioPMs",
            caption="Histograms of catch for each scen and sau across replicates")
    realyr <- c(glb$hyrnames,glb$pyrnames)[pickyrs[whichyr]]
    tabfile <- paste0("mean_catches_by_scenario_SAU_in_year_",realyr,".csv")
    addtable(res[whichyr,,],filen=tabfile,rundir=rundir,category="ScenarioPMs",
             caption=paste0("Mean catch of distribution of replicates by sau ",
             "and scenario in year ",realyr))
  }
  # Now do rate of change for 4 variables: catch, cpue, matureB, and harvestR
  pickvar <- "catch"
  res <- getrateofchange(dyn=dyn,whichvar=pickvar,glb=glb)
  filename <- plotrateofchange(rundir=rundir,res=res,whichvar=pickvar,glb=glb,
                               console=FALSE)
  addplot(filen=filename,rundir=rundir,category="ScenarioPMs",
          caption="The rate of change of the median catches acroiss scenarios.")
  pickvar <- "cpue"
  res <- getrateofchange(dyn=dyn,whichvar=pickvar,glb=glb)
  filename <- plotrateofchange(rundir=rundir,res=res,whichvar=pickvar,glb=glb,
                               console=FALSE)
  addplot(filen=filename,rundir=rundir,category="ScenarioPMs",
          caption="The rate of change of the median cpue acroiss scenarios.")
  pickvar <- "matureB"
  res <- getrateofchange(dyn=dyn,whichvar=pickvar,glb=glb)
  filename <- plotrateofchange(rundir=rundir,res=res,whichvar=pickvar,glb=glb,
                               console=FALSE)
  addplot(filen=filename,rundir=rundir,category="ScenarioPMs",
          caption="The rate of change of the median matureB across scenarios.")
  pickvar <- "harvestR"
  res <- getrateofchange(dyn=dyn,whichvar=pickvar,glb=glb)
  filename <- plotrateofchange(rundir=rundir,res=res,whichvar=pickvar,glb=glb,
                               console=FALSE)
  addplot(filen=filename,rundir=rundir,category="ScenarioPMs",
          caption="The rate of change of the median harvestR across scenarios.")
  # zone scale rates of change
  invar <- "catch"
  res <- getzonechangerate(zone=zone,whichvar=invar,glb=glb)
  filename <- plotzonechangerate(rundir=rundir,res=res,whichvar=invar,glb=glb,
                                 console=FALSE)
  addplot(filen=filename,rundir=rundir,category="ScenarioPMs",
          caption="Rate of change of the zonal median catch across scenarios.")
  invar <- "cpue"
  res <- getzonechangerate(zone=zone,whichvar=invar,glb=glb)
  filename <- plotzonechangerate(rundir=rundir,res=res,whichvar=invar,glb=glb,
                                 console=FALSE)
  addplot(filen=filename,rundir=rundir,category="ScenarioPMs",
          caption="Rate of change of the zonal median CPUE across scenarios.")
  invar <- "matureB"
  res <- getzonechangerate(zone=zone,whichvar=invar,glb=glb)
  filename <- plotzonechangerate(rundir=rundir,res=res,whichvar=invar,glb=glb,
                                 console=FALSE)
  addplot(filen=filename,rundir=rundir,category="ScenarioPMs",
          caption="Rate of change of the zonal median Mature Biomass across scenarios.")
  invar <- "harvestR"
  res <- getzonechangerate(zone=zone,whichvar=invar,glb=glb)
  filename <- plotzonechangerate(rundir=rundir,res=res,whichvar=invar,glb=glb,
                                 console=FALSE)
  addplot(filen=filename,rundir=rundir,category="ScenarioPMs",
          caption="Rate of change of the zonal median harvestR across scenarios.")
  # zone tab----------------------------------------
  if (verbose) cat("Now doing the Zone tab  \n")
  outzonephaseplot <- CBmsyphaseplot(rundir=rundir,zone=zone,prods=prods,
                                     glbc=glbc,category="zone",console=FALSE)
  outprod <- tabulatezoneprod(rundir,prods,scenes)
  zcatch <- makelist(scenes)
  for (scen in 1:nscen) zcatch[[scen]] <- zone[[scen]]$catch
  zoneribbon(rundir=rundir,scenes,invar=zcatch,glbc=glbc,varname="Catch",
             category="zone",console=FALSE,addmedian=0,addleg="bottomright")
  zcpue <- makelist(scenes)
  for (scen in 1:nscen) zcpue[[scen]] <- zone[[scen]]$cpue
  zoneribbon(rundir=rundir,scenes,invar=zcpue,glbc=glbc,varname="CPUE",
             category="zone",console=FALSE,addmedian=0,addleg="bottomright")
  zdeplsB <- makelist(scenes)
  for (scen in 1:nscen) zdeplsB[[scen]] <- zone[[scen]]$deplsB
  zoneribbon(rundir=rundir,scenes,invar=zdeplsB,glbc=glbc,varname="deplsB",
             category="zone",console=FALSE,addmedian=0,addleg="bottomright")
  plotzonedyn(rundir,scenes,zone,glbc,console=FALSE,q90=Q90,polys=TRUE,
              intens=intensity,hlines=list(catch=outprod[,"MSY"],
              spawnB=outprod[,"Bmsy"],harvestR=0,cpue=outprod[,"CEmsy"],
              expB=outprod[,"Bexmsy"]),scencol=scencol)

  # ribbon plots by sau and dynamic variable
  # Catch ribbon tab -------------------------------------------
  if (verbose) cat("Now doing the Catch tab  \n")
  catch <- scenebyvar(dyn,byvar="catch",glb=glbc,projonly = TRUE)
  catqnts <- sauquantbyscene(catch,glbc)
  sauribbon(rundir,scenes=scenes,varqnts=catqnts,
              glbc=glbc,varname="Catch",console=FALSE,
              q90=Q90,intens=intensity,addleg=ribbonleg)
  # cpue ribbon tab------------------------------------------
  if (verbose) cat("Now doing the cpue tab  \n")
  cpue <- scenebyvar(dyn,byvar="cpue",glb=glbc,projonly = TRUE)
  cpueqnts <- sauquantbyscene(cpue,glbc)
  sauribbon(rundir,scenes=scenes,varqnts=cpueqnts,
              glbc=glbc,varname="cpue",console=FALSE,
              q90=Q90,intens=intensity,addleg=ribbonleg)
  # depletion spawnB ribbon tab--------------------------------------------
  if (verbose) cat("Now doing the depletion spawnB tab  \n")
  deplsB <- scenebyvar(dyn,byvar="deplsB",glb=glbc,projonly=TRUE)
  deplsBqnts <- sauquantbyscene(deplsB,glbc)
  sauribbon(rundir,scenes=scenes,varqnts=deplsBqnts,
              glbc=glbc,varname="deplsB",console=FALSE,
              q90=Q90,intens=intensity,addleg=ribbonleg)
  # depletion exploitB ribbon tab--------------------------------------------
  if (verbose) cat("Now doing the depletion exploitB tab  \n")
  depleB <- scenebyvar(dyn,byvar="depleB",glb=glbc,projonly=TRUE)
  depleBqnts <- sauquantbyscene(depleB,glbc)
  sauribbon(rundir,scenes=scenes,varqnts=depleBqnts,
              glbc=glbc,varname="depleB",console=FALSE,
              q90=Q90,intens=intensity,addleg=ribbonleg)
  # harvest rate tab-------------------------------------------------
  if (verbose) cat("Now doing the harvestR tab  \n")
  harvestR <- scenebyvar(dyn,byvar="harvestR",glb=glbc,projonly=TRUE)
  hqnts <- sauquantbyscene(harvestR,glbc)
  sauribbon(rundir,scenes=scenes,varqnts=hqnts,
            glbc=glbc,varname="harvestR",console=FALSE,
            q90=Q90,intens=intensity,addleg=ribbonleg)
  # phaseplots tab --------------------------------------------------
  if (verbose) cat("Now doing the phaseplots tab  \n")
  plotallphaseplots(rundir=rundir,dyn=dyn,prods,glb=glbc[[1]],scenes=scenes,width=9,
                    height=10,fnt=7,pntcex=1.5,zero=FALSE,
                    legloc="topright")
  makecompareoutput(rundir=rundir,glbc,scenes,scenarionames,postfixdir,
                    filesused=files[pickfiles],openfile=TRUE,verbose=FALSE)
  return(invisible(list(scenes=scenes,quantscen=quantscen,dyn=dyn,
                        prods=prods,scenprops=scenprops,devout=devout)))
} # end of do_comparison

