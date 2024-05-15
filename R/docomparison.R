

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
#' @param ... the ellipsis is here in case a jurisdiction specific function or
#'     functions is/are written to perform extra analyses, plots, and tables.
#'
#' @seealso{
#'    \link{scenebyvar}, \link{scenebyzone}. \link{RGB}
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
                          juris="",ribbonleg="topleft",...) {
  # rundir=rundir;postfixdir=postfixdir;outdir=outdir;files=files;pickfiles=c(1,7)
  #  verbose=TRUE; intensity=100; zero=FALSE; altscenes=c("targCE","targ150")
  #  juris="";ribbonleg="topleft"
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
  for (i in 1:nfile) { # i = 1
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
    prods[[i]] <- t(out$sauprod)
    scenes[i] <- out$ctrl$runlabel
    scores[[i]] <- out$outhcr
    zone[[i]] <- out$outzone
  }
  # end results extraction----------------------
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
  projequal <- vector(mode="logical",length=(length(nscenes)-1))
  for (i in 2:nscenes) {
    projequal[i-1] <- all.equal(glbc[[i-1]]$pyrs,glbc[[i]]$pyrs)
    if (!(all(projequal))) {
      label <- paste0("Number of projection years differs among scenario")
      warning(cat("Number of projection years differs among scenario"))
    }
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
    warning(cat("At least one important scenario property differs betweem scenarios \n"))
  }
  if (tolower(juris) == "tas") {
     cat("A specific jurisdiction fnction is required \n")
  }
  if (verbose) print(scenprops)
  if (verbose) cat("Now doing the comparisons  \n")
  setuphtml(rundir=rundir)
  filename <- "scenarioproperties.csv"
  addtable(scenprops,filen=filename,rundir=rundir,category="scenes",
           caption=paste0("Important scenario properties for comparability."))
  # dynamics tab-----------------------------------------
  quantscen <- comparedynamics(rundir=rundir,dyn,glbc,scenes)
  # productivity tab---------------------------------
  tabulateproductivity(rundir,prods,scenes)
  # Final Scores tab---------------------------------
  scenscore <- lapply(scores,"[[","finalsc")
  if (!is.null(scenscore[[1]])) {
    filename <- "compare_final_HSscores.png"
    meds <- comparefinalscores(rundir,scores,scenes,legloc="bottomright",
                               filen=filename,category="Scores")
    addplot(filen=filename,rundir=rundir,category="Scores",
            caption="The HS final scores for each sau.")
    tabulatefinalHSscores(rundir,meds,scenes,category="Scores")
  } else {
    warning(cat("TasHS should be first in comparisons, if used,",
                " no Score tab produced \n"))
  }
  # HSPM tab--------------------------------------------
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
  cdivmsy <- catchvsMSY(catch,glbc,prods,scenes)
  nscen <- length(scenes)
  for (i in 1:nscen) {
    plotsceneproj(rundir,cdivmsy[[i]],glbc[[i]],scenes[i],
                  filen=paste0("Catch_div_MSY_",scenes[i],".png"),
                  label="Catch / MSY",hline=1,Q=90)
  }
  # zone tab----------------------------------------
  outprod <- tabulatezoneprod(rundir,prods,scenes)
  plotzonedyn(rundir,scenes,zone,glbc,console=FALSE,q90=Q90,polys=TRUE,
              intens=intensity,hlines=list(catch=outprod[,"MSY"],
              spawnB=outprod[,"Bmsy"],harvestR=0,cpue=outprod[,"CEmsy"]))

  # ribbon plots by sau and dynamic variable
  # cpue <- scenebyvar(dyn=out$dyn,byvar="cpue",glb=out$glbc[[1]])
  # Catch ribbon tab -------------------------------------------
  glb <- glbc[[1]]
  catch <- scenebyvar(dyn,byvar="catch",glb=glb,projonly = TRUE)
  catqnts <- sauquantbyscene(catch,glb)
  nsau <- glbc[[1]]$nSAU
  for (sau in 1:nsau)
    sauribbon(rundir,scenes=scenes,sau=sau,varqnts=catqnts,
              glb=glb,varname="Catch",console=FALSE,
              q90=Q90,intens=intensity,addleg=ribbonleg)
  # cpue ribbon tab------------------------------------------
  cpue <- scenebyvar(dyn,byvar="cpue",glb=glb,projonly = TRUE)
  cpueqnts <- sauquantbyscene(cpue,glb)
  for (sau in 1:nsau)
    sauribbon(rundir,scenes=scenes,sau=sau,varqnts=cpueqnts,
              glb=glb,varname="cpue",console=FALSE,
              q90=Q90,intens=intensity,addleg=ribbonleg)
  # depletion spawnB ribbon tab--------------------------------------------
  deplsB <- scenebyvar(dyn,byvar="deplsB",glb=glb,projonly=TRUE)
  deplsBqnts <- sauquantbyscene(deplsB,glb)
  for (sau in 1:nsau)
    sauribbon(rundir,scenes=scenes,sau=sau,varqnts=deplsBqnts,
              glb=glb,varname="deplsB",console=FALSE,
              q90=Q90,intens=intensity,addleg=ribbonleg)
  # depletion exploitB ribbon tab--------------------------------------------
  depleB <- scenebyvar(dyn,byvar="depleB",glb=glb,projonly=TRUE)
  depleBqnts <- sauquantbyscene(depleB,glb)
  for (sau in 1:nsau)
    sauribbon(rundir,scenes=scenes,sau=sau,varqnts=depleBqnts,
              glb=glb,varname="depleB",console=FALSE,
              q90=Q90,intens=intensity,addleg=ribbonleg)
  # phaseplots tab --------------------------------------------------
  plotallphaseplots(rundir=rundir,dyn=dyn,prods,glb=glb,scenes=scenes,width=9,
                    height=10,fnt=7,pntcex=1.5,zero=FALSE,
                    legloc="topright")
  makecompareoutput(rundir=rundir,glbc,scenes,scenarionames,postfixdir,
                    filesused=files[pickfiles],openfile=TRUE,verbose=FALSE)
  return(invisible(list(scenes=scenes,ans=ans,quantscen=quantscen,dyn=dyn,
                        prods=prods,scenprops=scenprops)))
} # end of do_comparison

