
# hutils::listFunctions("C:/Users/User/Dropbox/A_code/aMSE/R/inputfiles.R")


#' Title
#' @param rundir the directory in which all files relating to a
#'     particular run are to be held.
#' @param controlfile default="control.csv", the filename of the control
#'     file present in rundir containing information regarding the run.
#' @param verbose Should progress comments be printed to console, default=TRUE#'
#'
#' @return nothing but it does plot a graph
#' @export
#'
#' @examples
#' print("wait on data sets")
checksizecompdata <- function(rundir,controlfile,verbose=TRUE) {
  zone1 <- readctrlfile(rundir,infile=controlfile,verbose=verbose)
  setuphtml(rundir)
  compdat <- zone1$condC$compdat$lfs
  if (is.null(compdat))
    stop(cat("No size-composition file found in copntrol file \n"))
  palfs <- zone1$condC$compdat$palfs
  saunames <- zone1$SAUnames
  nsau <- length(saunames)
  histyr <- zone1$condC$histyr
  setuphtml(rundir)
  addtable(palfs,"sizecomp_obs_by_year_by_SAU.csv",rundir,
           category="predictedcatchN",
           caption="Number of sizecomp observations by year and SAU.")
  for (plotsau in 1:nsau) {
    lfs <- preparesizecomp(compdat[,,plotsau],mincount=0)
    yrsize <- as.numeric(colnames(lfs))
    pickyr <- match(yrsize,histyr[,"year"])
    LML <- histyr[pickyr,"histLML"]
    plotsizecomp(rundir=rundir,incomp=lfs,SAU=saunames[plotsau],lml=LML,
                 catchN=NULL,start=NA,proportion=TRUE,
                 console=FALSE)
  }
    replist <- NULL
    projy <- 0
    runnotes <- c(controlfile,paste0("SAU = ",nsau))
    make_html(replist = replist,  rundir = rundir,
              controlfile=controlfile, datafile=zone1$ctrl$datafile,
              width = 500, openfile = TRUE,  runnotes = runnotes,
              verbose = verbose, packagename = "aMSE",  htmlname = "sizecomp")
    if (verbose) cat("finished  \n")
} # end of checksizecompdata

#' @title checkmsedata contains some tests of the niput MSE saudata file
#'
#' @param intxt the data file from readLines used in readsaudatafile
#' @param rundir the rundir for the scenario
#' @param verbose should test results be put to the console as well as filed,
#'     default=TRUE
#'
#' @seealso{
#'  \link{readsaudatafile}
#' }
#'
#' @return nothing but will write a file to rundir and may write to the console
#' @export
#'
#' @examples
#' print("wait on example data sets")
checkmsedata <- function(intxt,rundir,verbose=TRUE) { # intxt=indat;verbose=TRUE
  filen <- filenametopath(rundir,"input_data_tests.txt")
  cat("Input data test \n\n",file=filen,append=FALSE)
  endrun <- FALSE
  pickwta <- grep("Wtbtoa",intxt)[1]
  if (!is.na(pickwta)) {
    txt <- "datafile uses deprecated Wtbtoa input variable instead of Wta  \n"
    if (verbose) cat(txt)
    cat(txt,file=filen,append=TRUE)
    endrun=TRUE
  }
  pickL <- grep("lambda",intxt)[1]
  if (pickL == 0) {
    txt <- "datafile does not have a line containing a lambda value  \n"
    if (verbose) cat(txt)
    cat(txt,file=filen,append=TRUE)
    endrun=TRUE
  }
  pickL <- grep("qest",intxt)[1]
  if (pickL == 0) {
    txt <- "datafile does not have a line containing a qest value  \n"
    if (verbose) cat(txt)
    cat(txt,file=filen,append=TRUE)
    endrun=TRUE
  }
  # if (nrow(poprec) != numpop) {
  #   txt <- "Number of population recruitment proportions different from numpop  \n"
  #   if (verbose) cat(txt)
  #   cat(txt,file=filen,append=TRUE)
  #   endrun=TRUE
  # }
  if (endrun) stop(cat("See input_data_tests.txt for data errors \n"))
} # end of checkmsedata



#' @title ctrlfiletemplate generates a template input control file
#'
#' @description ctrlfiletemplate generates a standardized control file
#'     template. Generate this and then modify the contents to suit
#'     the system you are attempting to simulate. Defaults to 100
#'     replicates. There needs to be as many recdevs as there are conditioning
#'     catches (in the template = 58 rows. If they are all set to -1,
#'     the default, then they will have no effect and recdevs will be taken off
#'     the stock-recruitment curve with random deviates defined by withsigR
#'     found in the ctrl file. If devrec = 1.0 then the recdevs from 1980 to
#'     2016 will be set to 1.0 which means there will be no random variation
#'     and the system is ready to be conditioned on those recdevs to improve
#'     the predicted CPUE during the historical conditioning period. If
#'     devrec = 0.0, then the recdevs for 1980 - 2016 will be set to values
#'     that have already conditioned the model and provides a reasonable fit to
#'     the observed CPUE trends.
#'
#' @param indir directory in which to place the control.csv file
#' @param filename the name for the generated ctrlfile, a character
#'     string that defaults to controlsau.csv.
#' @param devrec what form should the recdevs take? valid values are -1, the
#'     default (= random deviates imposed dependent upon withsigR), 0.0, where
#'     a set of conditioned values between 1980 and 2016 are imposed that
#'     improve the fit to the observed CPUE as inserted by ctrlfiletemplate,
#'     or 1.0, which prepares the recdevs for a round of conditioning ready to
#'     improve the observed fit to cpue.
#'
#' @seealso{
#'   \link{datafiletemplate}, \link{readsaudatafile}, \link{readctrlfile}
#' }
#'
#' @return invisibly the fill path and name of the control file. More
#'     importantly, it write a control file template to that directory.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  yourdir <- tempdir()
#'  datafiletemplate(nSAU=8,yourdir,"saudata_test.csv")
#'  ctrlfiletemplate(yourdir,filename="testctrl.csv",devrec=-1)   #
#'  control <- readctrlfile(yourdir,"testctrl.csv")
#'  str(control,max.level=1)
#' }
ctrlfiletemplate <- function(indir,filename="controlsau.csv",devrec=-1) { # indir=rundir; filename="control2.csv"
   if (!(devrec %in% c(-1,0.0,1.0)))
      stop("devrec must be either -1, 0.0, or 1 in ctrlfiletemplate  \n")
   filename <- filenametopath(indir,filename)
   cat("DESCRIPTION \n",
       file=filename,append=FALSE)
   cat("Control file containing details of a particular run. Modify the  \n",
       file=filename,append=TRUE)
   cat("contents to suit your own situation. In particular. modify the  \n",
       file=filename,append=TRUE)
   cat("contents of this description to suit the scenario being executed \n",
       file=filename,append=TRUE)
   cat("\n",file=filename,append=TRUE)
   cat("START \n",file=filename,append=TRUE)
   cat("runlabel, a_meaningful_name, label for particular run \n",
       file=filename,append=TRUE)
   cat("datafile, saudata_test.csv, name of saudefs file \n",
       file=filename,append=TRUE)
   cat("bysau, 1, 1=TRUE and 0=FALSE  \n",file=filename,append=TRUE)
   cat("parfile, NA, name of file containing optimum parameters from sizemod \n",
       file=filename,append=TRUE)
   cat("\n\n",file=filename,append=TRUE)
   cat("zoneCOAST \n",file=filename,append=TRUE)
   cat("replicates,  100, number of replicates, usually 250, 500 or more  \n",
       file=filename, append=TRUE)
   cat("withsigR,  0.35, recruitment variability eg 0.5 \n",
       file=filename, append=TRUE)
   cat("withsigB,  0.1, process error on exploitable biomass \n",
       file=filename, append=TRUE)
   cat("withsigCE, 0.1, process error on cpue calculations  \n",
       file=filename, append=TRUE)
   cat("\n",file=filename,append=TRUE)
   cat("ZONE \n",file=filename,append=TRUE)
   cat("nSAU, 8, number of spatial management units eg 2 \n",
       file=filename,append=TRUE)
   cat("SAUpop, 3, 3, 5, 7, 9, 9, 12, 8, number of populations per SAU in sequence \n",
       file=filename,append=TRUE)
   cat("SAUnames, sau6, sau7, sau8, sau9, sau10, sau11, sau12, sau13, labels for each SAU \n",
       file=filename,append=TRUE)   # possibly deprecated
   cat("initdepl, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, initial depletion levels for each SAU \n",
       file=filename,append=TRUE)  # Only the control file initdepl is used
   cat("\n\n",file=filename,append=TRUE)
   cat("SIZE \n",file=filename,append=TRUE)
   cat("minc, 2, centre of minimum size class \n",file=filename,append=TRUE)
   cat("cw, 2, class width mm \n",file=filename,append=TRUE)
   cat("Nclass, 105, number of size classes \n",file=filename,append=TRUE)
   cat("\n",file=filename,append=TRUE)
   cat("RECRUIT \n",file=filename,append=TRUE)
   cat("larvdisp, 0.01, rate of larval dispersal eg 0.01 = 0.5 percent in either direction\n",
       file=filename,append=TRUE)
   cat("\n",file=filename,append=TRUE)
   cat("RANDOM \n",file=filename,append=TRUE)
   cat("randomseed, 3543304, for repeatability of population definitions set to 0 otherwise \n",
       file=filename,append=TRUE)
   cat("randomseedP, 0, for repeatability of projections set to 0 otherwise \n",
       file=filename,append=TRUE)
   cat("\n",file=filename,append=TRUE)
   cat("initLML, 140, initial LML for the unfished zone if no historical catches used \n",
       file=filename,append=TRUE)  # deprecated
   cat("\n",file=filename,append=TRUE)
   cat("PROJECT, 30, number of projection years for each simulation \n",
       file=filename,append=TRUE)
   cat("\n\n",file=filename,append=TRUE)
   cat("PROJLML, need the same number as there are projection years \n",file=filename,
       append=TRUE) # ensure there are Nyrs lines
   cat("2021, 145,  Legal Minimum Length (LML, MLL, MLS) e.g. 140 \n",
       file=filename,append=TRUE)
   for (i in 2:30) {
      yr <- 2021 + i - 1
      cat(as.character(yr),", 145 \n",file=filename,append=TRUE)
   }
   cat("\n",file=filename,append=TRUE)
   cat("CATCHES, 58, if > 1 then how many years in the histLML \n",
       file=filename,append=TRUE)
   cat("\n",file=filename,append=TRUE)
   cat("CondYears, LML, 6,7,8,9,10, 11, 12, 13 \n",file=filename,append=TRUE)
   cat("1963,127,0,0,0,0,0,0,0,0, \n",file=filename,append=TRUE)
   cat("1964,127,1,1,1,4,3, 5, 4,1, \n",file=filename,append=TRUE)
   cat("1965,127,2,3,4,17,15,21,19,5, \n",file=filename,append=TRUE)
   cat("1966,127,8,11,13,56,50, 72, 63,17, \n",file=filename,append=TRUE)
   cat("1967,127,22,31,34,149,133,191,169,46, \n",file=filename,append=TRUE)
   cat("1968,127,31,43,47,208,185,266,235,63, \n",file=filename,append=TRUE)
   cat("1969,127,24,33,36,160,142,204,180,49, \n",file=filename,append=TRUE)
   cat("1970,127,29,40,44,195,173,249,220,59, \n",file=filename,append=TRUE)
   cat("1971,127, 39,53,59,261,231,333,294,79, \n",file=filename,append=TRUE)
   cat("1972,127, 33,45,50,222,197,284,250,68, \n",file=filename,append=TRUE)
   cat("1973,127, 24,33,37,162,144,207,183,49, \n",file=filename,append=TRUE)
   cat("1974,127, 23,32,36,157,139,201,177,48, \n",file=filename,append=TRUE)
   cat("1975,127,21, 36, 42,126, 130, 191,143,43, \n",file=filename,append=TRUE)
   cat("1976,127,30,56,77,252,179,240,153,62, \n",file=filename,append=TRUE)
   cat("1977,127,19,24,22,123,98,153,189,39, \n",file=filename,append=TRUE)
   cat("1978,127,28,13,27,115,258,275,208,58, \n",file=filename,append=TRUE)
   cat("1979,127,31,19,23,172,166,269,325,63, \n",file=filename,append=TRUE)
   cat("1980,127,42,81,63,316,195,338,351,87, \n",file=filename,append=TRUE)
   cat("1981,127,48,88,87,444,260,417,246,99, \n",file=filename,append=TRUE)
   cat("1982,127,30,34,34,249,100,303,235,62, \n",file=filename,append=TRUE)
   cat("1983,127,38,102,58,199,175,430,242,78, \n",file=filename,append=TRUE)
   cat("1984,127,50,78,38,248,284,682,258,102, \n",file=filename,append=TRUE)
   cat("1985,127,36,99,23,246,140,479,155,74, \n",file=filename,append=TRUE)
   cat("1986,127,27,97,11,133,127,289,194,55, \n",file=filename,append=TRUE)
   cat("1987,132,31,84,44,251,82,339,195,64, \n",file=filename,append=TRUE)
   cat("1988,132,25,53,27,160,126,276,162,52, \n",file=filename,append=TRUE)
   cat("1989,132,21,49,46,120,109,212,145,44, \n",file=filename,append=TRUE)
   cat("1990,140,19,56,21,95,80,232,125,39, \n",file=filename,append=TRUE)
   cat("1991,140,20,54,30,102,106,219,140,42, \n",file=filename,append=TRUE)
   cat("1992,140,23,70,40,102,96,267,160,47, \n",file=filename,append=TRUE)
   cat("1993,140,21,67,40,110,65,197,177,42, \n",file=filename,append=TRUE)
   cat("1994,140,18,37,38,78,60,203,160,37, \n",file=filename,append=TRUE)
   cat("1995,140,17,33,17,44,69,186,185,34, \n",file=filename,append=TRUE)
   cat("1996,140,16,68,13,59,75,145,148,33, \n",file=filename,append=TRUE)
   cat("1997,140,24,75,28,140,66,224,227,49, \n",file=filename,append=TRUE)
   cat("1998,140,18,50,27,78,47,165,204,37, \n",file=filename,append=TRUE)
   cat("1999,140,23,60,24,115,58,220,251,47, \n",file=filename,append=TRUE)
   cat("2000,140,21,61,23,205,148,333,288,54, \n",file=filename,append=TRUE)
   cat("2001,140,49,32,15,186,157,321,304,43, \n",file=filename,append=TRUE)
   cat("2002,140,31,52,17,174,149,366,237,93, \n",file=filename,append=TRUE)
   cat("2003,140,34,104,27,142,246,352,232,67, \n",file=filename,append=TRUE)
   cat("2004,140,24,89,22,130,193,390,250,96, \n",file=filename,append=TRUE)
   cat("2005,140,26,110,26,92,149,389,311,65, \n",file=filename,append=TRUE)
   cat("2006,140,50,75,6,143,198,384,229,89, \n",file=filename,append=TRUE)
   cat("2007,140,34,39,18,178,231,354,267,68, \n",file=filename,append=TRUE)
   cat("2008,140,35,51,9,156,178,345,305,79, \n",file=filename,append=TRUE)
   cat("2009,140,47,107,51,155,110,244,327,77, \n",file=filename,append=TRUE)
   cat("2010,140,23, 110,37,159,158,245,278,68, \n",file=filename,append=TRUE)
   cat("2011,140,17,95,48,171,159,247,257,56, \n",file=filename,append=TRUE)
   cat("2012,140,59,97,19,172,146,273,267,44, \n",file=filename,append=TRUE)
   cat("2013,140,11,44,8,158,180,288,251,41, \n",file=filename,append=TRUE)
   cat("2014,140,34,44,5,98,142,220,249,37, \n",file=filename,append=TRUE)
   cat("2015,140,32,65,11,89,115,240,245,34, \n",file=filename,append=TRUE)
   cat("2016,140,19,31,12,62,77,170,299,34, \n",file=filename,append=TRUE)
   cat("2017,140,12,37,7, 56,78,195,274,42, \n",file=filename,append=TRUE)
   cat("2018,140,11,41,14,57,89,174,262,52, \n",file=filename,append=TRUE)
   cat("2019,140,16,53,5,65,81,179,251,53, \n",file=filename,append=TRUE)
   cat("2020,145, 7,27,6,39,58,129,227,50, \n",file=filename,append=TRUE)
   cat("\n",file=filename,append=TRUE)
   cat(" , 6,7,8,9,10,11,12,13  \n",file=filename,append=TRUE)
   cat("CEYRS, 29,  if >1 then number reflects number of historical CPUE records by SAU \n",
       file=filename,append=TRUE)
   cat("1992,,113.40507,94.24715,97.03372,99.43445,98.10875,100.19948,  \n",
       file=filename,append=TRUE)
   cat("1993,,116.80512,110.3889,109.37573,99.14548,107.88068,102.2077,  \n",
       file=filename,append=TRUE)
   cat("1994,,104.02704,114.1891,136.78298,122.78214,110.6589,101.31698,  \n",
       file=filename,append=TRUE)
   cat("1995,,114.25513,130.97739,134.70459,135.68894,122.81239,116.48724,  \n",
       file=filename,append=TRUE)
   cat("1996,,149.13849,189.5771,174.3451,152.3415,142.48358,134.99475,  \n",
       file=filename,append=TRUE)
   cat("1997,,159.97343,197.61029,168.36235,170.04146,152.4425,142.20638,  \n",
       file=filename,append=TRUE)
   cat("1998,,165.45067,193.0997,175.59561,187.56223,157.19827,141.65122,  \n",
       file=filename,append=TRUE)
   cat("1999,,173.87182,180.72641,193.4431,177.77114,170.55928,142.87788,  \n",
       file=filename,append=TRUE)
   cat("2000,173.36549,179.96736,166.54089,180.31256,196.45518,169.38026,143.01979,95.87332  \n",
       file=filename,append=TRUE)
   cat("2001,157.15512,175.30658,165.32838,167.92169,184.64006,152.99312,133.48884,93.99618  \n",
       file=filename,append=TRUE)
   cat("2002,137.60701,175.1555,197.57349,168.51901,187.69954,160.51576,123.51065,88.34786  \n",
       file=filename,append=TRUE)
   cat("2003,122.49438,159.56903,150.9974,156.1052,176.07989,145.71378,125.38706,88.87757  \n",
       file=filename,append=TRUE)
   cat("2004,115.2121,157.85257,141.69066,147.15008,166.85616,150.75278,110.18178,87.60967  \n",
       file=filename,append=TRUE)
   cat("2005,116.15036,152.5276,152.91814,155.24376,154.66732,142.1134,110.44022,89.25662  \n",
       file=filename,append=TRUE)
   cat("2006,117.17298,150.05658,175.7411,155.69089,168.4007,137.18411,113.4301,93.56854  \n",
       file=filename,append=TRUE)
   cat("2007,122.36628,165.94552,185.97047,157.3592,159.0966,130.29756,115.57321,94.061  \n",
       file=filename,append=TRUE)
   cat("2008,143.03856,181.61195,186.72795,151.78263,134.57655,119.60881,109.86311,99.9658  \n",
       file=filename,append=TRUE)
   cat("2009,132.29504,172.65738,173.13257,183.14568,137.51885,131.24645,120.07897,93.7539  \n",
       file=filename,append=TRUE)
   cat("2010,129.24964,143.20195,154.61341,160.46221,162.7066,132.48497,112.06231,101.15461  \n",
       file=filename,append=TRUE)
   cat("2011,127.0115,137.70191,157.92823,167.85705,153.21758,125.63879,117.32218,100.25653  \n",
       file=filename,append=TRUE)
   cat("2012,111.25489,124.41648,123.58717,138.91591,132.25876,121.88758,112.92296,95.33387  \n",
       file=filename,append=TRUE)
   cat("2013,103.80147,117.47827,149.88607,124.90246,121.07735,105.05426,99.62491,76.17599  \n",
       file=filename,append=TRUE)
   cat("2014,97.26094,102.64518,92.66122,103.28205,100.21447,97.79445,101.16735,73.88172  \n",
       file=filename,append=TRUE)
   cat("2015,81.42042,98.81462,105.34779,95.86434,86.34135,86.20905,96.06224,74.07834  \n",
       file=filename,append=TRUE)
   cat("2016,79.19767,100.96856,144.2149,100.76882,88.3244,90.96269,104.11514,84.90243  \n",
       file=filename,append=TRUE)
   cat("2017,90.25977,146.57429,166.5622,102.01255,96.66777,97.93013,105.64149,97.49921  \n",
       file=filename,append=TRUE)
   cat("2018,106.68431,143.5766,148.43001,133.7072,104.3823,95.01433,106.88331,104.08634  \n",
       file=filename,append=TRUE)
   cat("2019,92.22078,112.70921,107.11476,93.65297,91.89809,86.94485,90.03483,90.78764  \n",
       file=filename,append=TRUE)
   cat("2020,92.01,113.01,112.1,99.1,92.1,87.1,93.1,92.1  \n",
       file=filename,append=TRUE)
   cat("\n",file=filename,append=TRUE)
   cat("SIZECOMP, 1, if 1 then a single filename should follow, if 2 then two filenames, etc. \n",
       file=filename,append=TRUE)
   cat("lf_WZ90-20.csv \n",file=filename,append=TRUE)
   cat("\n",file=filename,append=TRUE)
   cat("\n",file=filename,append=TRUE)
   cat("RECDEV, 58  \n",file=filename,append=TRUE)
   cat("CondYears,6,7,8,9,10,11,12,13 \n",file=filename,append=TRUE)
   if (devrec < 0) {
      for (i in 1:58) {
         yr <- 1962 + i
         cat(as.character(yr),",-1,-1,-1,-1,-1,-1,-1,-1 \n",file=filename,append=TRUE)
      }
   }
   if (devrec == 1) {
      for (i in 1:17)
         cat(as.character((1962+i)),",-1,-1,-1,-1,-1,-1,-1,-1 \n",file=filename,append=TRUE)
      for (i in 18:54)
         cat(as.character((1962+i)),",1,1,1,1,1,1,1,1 \n",file=filename,append=TRUE)
      for (i in 55:58)
         cat(as.character((1962+i)),",-1,-1,-1,-1,-1,-1,-1,-1 \n",file=filename,append=TRUE)
   }
   if (devrec == 0){
      for (i in 1:19)
         cat(as.character((1962+i)),",-1,-1,-1,-1,-1,-1,-1,-1 \n",file=filename,append=TRUE)
     cat("1982,1,0.898011706,1.214738587,1.130112683,0.976241769,0.893551239,0.790051338,1  \n",
         file=filename,append=TRUE)
     cat("1983,1,0.924020966,1.061569181,1.016748193,0.97407056,0.855462095,0.807254058,1  \n",
         file=filename,append=TRUE)
     cat("1984,1,0.904397602,0.913885752,0.851196819,0.945096581,0.828747968,0.84114442,1  \n",
         file=filename,append=TRUE)
     cat("1985,1,0.887284703,0.850845565,0.834720081,0.906852318,0.865954227,0.854925667,1  \n",
         file=filename,append=TRUE)
     cat("1986,1,1.044487535,1.011602322,0.870319748,0.878467782,0.889352256,0.898698496,1  \n",
         file=filename,append=TRUE)
     cat("1987,1,1.019790172,1.248335334,1.244829755,0.893806461,1.238162698,1.028010868,1  \n",
         file=filename,append=TRUE)
     cat("1988,1,0.869125082,1.603827547,1.465344761,1.320242475,1.083471531,1.077481418,1  \n",
         file=filename,append=TRUE)
     cat("1989,1,1.131831212,1.503945427,1.088430947,1.36101713,1.309478202,1.946136712,1  \n",
         file=filename,append=TRUE)
     cat("1990,1.026883055,2.354756479,1.091624333,1.302493513,1.4294116,1.571634174,1.948869064,1.112995867  \n",
         file=filename,append=TRUE)
     cat("1991,1.025504393,1.592403174,0.94517318,1.12933916,1.480964047,1.305795389,1.370702218,1.169983501  \n",
         file=filename,append=TRUE)
     cat("1992,1.003477151,1.288926824,0.909627819,1.179802751,1.293363588,1.211025988,1.307728304,1.069162779  \n",
         file=filename,append=TRUE)
     cat("1993,0.950446781,1.392466919,0.946588265,1.384965705,1.087584767,1.454911166,1.455091766,1.009196818  \n",
         file=filename,append=TRUE)
     cat("1994,0.907463441,1.313074622,0.894908768,1.292202188,1.378920757,1.227064236,1.392855906,1.155102507  \n",
         file=filename,append=TRUE)
     cat("1995,0.858233289,1.268824542,0.604808622,1.012616194,1.148739557,1.011831269,1.013057965,1.286957479  \n",
         file=filename,append=TRUE)
     cat("1996,0.871778225,1.061975292,0.499616399,1.031537226,1.374171074,1.415594271,1.011898488,1.362088585  \n",
         file=filename,append=TRUE)
     cat("1997,0.974371042,1.102884484,0.703406095,1.023172864,1.057765181,1.10798134,0.948543997,1.620690891  \n",
         file=filename,append=TRUE)
     cat("1998,1.27313887,1.436630006,0.859824285,0.906421556,1.4294448,1.500686209,0.91166883,1.610706733  \n",
         file=filename,append=TRUE)
     cat("1999,1.709368379,1.510815833,0.977921567,0.862098634,1.216520334,1.146273685,1.403491749,1.338012311  \n",
         file=filename,append=TRUE)
     cat("2000,1.705211324,1.309104788,1.067001987,0.84889326,1.403527711,1.30775231,1.303527061,1.553556795  \n",
         file=filename,append=TRUE)
     cat("2001,1.314174499,1.268457964,0.870139431,1.303824827,1.246590614,1.195627718,1.302389947,1.365563275  \n",
         file=filename,append=TRUE)
     cat("2002,1.226713417,1.737888211,0.882613276,1.071529113,0.855268209,0.997172623,1.432902852,1.225570786  \n",
         file=filename,append=TRUE)
     cat("2003,0.943464143,1.066059689,0.944418721,1.143858121,1.225153199,1.271127418,1.674412001,1.028644968  \n",
         file=filename,append=TRUE)
     cat("2004,1.000904676,1.251436682,0.873038153,1.11663499,1.456359522,1.010730368,1.270680622,0.926752604  \n",
         file=filename,append=TRUE)
     cat("2005,1.199028024,1.419415417,0.657370317,1.115651184,1.136964872,0.991910748,1.619217133,0.66165756  \n",
         file=filename,append=TRUE)
     cat("2006,0.71316976,1.061463222,0.510491095,0.876820237,0.91626714,0.871127136,1.104387063,0.537637698  \n",
         file=filename,append=TRUE)
     cat("2007,0.606889001,0.754426065,0.420885208,0.848747182,0.715395043,0.803575214,1.137049936,0.586222133  \n",
         file=filename,append=TRUE)
     cat("2008,0.822076796,0.545765138,0.48530799,0.656962161,0.616700458,0.760318997,0.683776539,0.96392894  \n",
         file=filename,append=TRUE)
     cat("2009,1.071538265,0.6979601,0.614021871,0.4319214,0.6003162,0.529756987,0.783770727,1.638943019  \n",
         file=filename,append=TRUE)
     cat("2010,0.720795297,1.259884354,0.595582076,0.547748764,0.817399374,0.901720054,2.235830081,1.248255152  \n",
         file=filename,append=TRUE)
     cat("2011,0.588408755,1.463451263,0.570265262,0.660606001,0.878514687,1.064602073,1.260317538,1.0177701  \n",
         file=filename,append=TRUE)
     cat("2012,0.842392164,0.821685524,0.626249726,0.67382212,0.843745396,0.729183606,1.148005799,1.173249244  \n",
         file=filename,append=TRUE)
     cat("2013,0.994967427,0.75906534,0.865777016,0.508232173,0.748440695,0.725781167,1.13047107,1.089817397  \n",
         file=filename,append=TRUE)
     cat("2014,0.998185669,0.953031116,1.01835623,0.748138628,0.898743845,0.844394757,1.209967431,1.003727701  \n",
         file=filename,append=TRUE)
      for (i in 53:58)
         cat(as.character((1962+i)),",-1,-1,-1,-1,-1,-1,-1,-1 \n",file=filename,append=TRUE)
   }
   cat("\n",file=filename,append=TRUE)
   return(invisible(filename))
} # end of ctrlzonetemplate


#' @title datafiletemplate generates a template input datafile akin to M15h75
#'
#' @description datafiletemplate generates a standard input datafile
#'     to use as a template, go in and edit it appropriately to suit
#'     your own needs. It contains the probability distributions that
#'     are sampled to provide the necessary biological constants for
#'     each population. It also contains the proportional distribution of
#'     recruitment levels determined for Tasmania by examining the GPS
#'     data-logger data for the zone and allocating recruitment in
#'     proportion to the relative yield by area over the last 8 years.
#'
#' @param indir the directory into which to place the generated data file
#' @param filename the name for the generated datafile, a character
#'     string, defaults to saudata_test.csv, which is the default within
#'     the ctrlfiletemplate function
#'
#' @seealso{
#'   \link{ctrlfiletemplate}, \link{readsaudatafile}, \link{readctrlfile}
#' }
#'
#' @return a standard definition data file, to be read by readdatafile
#'     whose name and path is returned invisibly
#' @export
#'
#' @examples
#' \dontrun{
#'  yourdir <- tempdir()
#'  datafiletemplate(yourdir,"saudata_test.csv")
#'  constants <- readsaudatafile(yourdir,"saudata_test.csv")
#'  str(constants,max.level=1)
#'  print(constants[,1:10])
#' }
datafiletemplate <- function(indir,filename="saudata_test.csv") {
   genconst <- function(invect) { # a function to generate constant values
      nlab <- length(invect)
      invect <- as.character(invect)
      ans <- invect[1]
      if (nlab > 1) for (i in 2:nlab) ans <- paste(ans,invect[i],sep=",")
      return(ans)
   } # end of genconst
   filename <- filenametopath(indir,filename)
   nSAU <- 8  # this template generates data for an 8 sau zone
   cat("SAU definitions listing Probability density function parameters for each variable \n",
       file=filename,append=FALSE)
   cat("These are randomly aligned to each population except for the proportion of recruitment \n",
       file=filename,append=TRUE)
   cat("which is literally allocated down in the popREC section. Average recruitment has a tiny \n",
       file=filename,append=TRUE)
   cat("sd which is to aid in simplifying conditioning. \n",
       file=filename,append=TRUE)
   cat("\n\n",file=filename,append=TRUE)
   cat("SPATIAL \n",file=filename,append=TRUE)
   cat("nsau, 8, number of spatial management units  \n",file=filename,append=TRUE)
   cat("saupop, 3,3,5,7,9,9,12,8, number of populations per SAU in sequence  \n",
       file=filename,append=TRUE)
   cat("saunames, 6,7,8,9,10,11,12,13, labels for each SAU  \n",file=filename,append=TRUE)
   # cat("initdepl,",genconst(rep(0.7,nSAU)),", initial depletion levels for each SAU  \n",
   #     file=filename,append=TRUE)  # possibly deprecated this will be clarified
   cat("\n\n",file=filename,append=TRUE)
   cat("PDFs,  32   \n",file=filename, append=TRUE)
   cat("DLMax , 21.5060321746972,26.6005907815097,24.3520350566157,29.0972888599621,",
       "29.2904865501258,28.6097710536279,26.7206918605892,21.5736807435363 ,",
        "maximum growth increment \n",file=filename, append=TRUE)
   cat("sMaxDL ,",genconst(rep(1e-04,nSAU)),", variation of MaxDL \n",
       file=filename, append=TRUE)
   cat("L50 ,",genconst(rep(130.0,nSAU)),", Length at 50% MaxDL \n",
       file=filename, append=TRUE)
   cat("sL50 ,",genconst(rep(2e-04,nSAU)),", variation of L50 \n",
       file=filename, append=TRUE)
   cat("L50inc , 54.3495728192894,43.1978234068735,53.7880354443096,49.2417858122113,",
       "53.621543719204,51.5238879305982,48.7727579990218,54.2813866497762,",
       "L95 - L50 = delta =L50inc \n",file=filename, append=TRUE)
   cat("sL50inc ,",genconst(rep(1e-04,nSAU)),", variation of L50inc \n",
       file=filename, append=TRUE)
   cat("SigMax  ,",genconst(rep(3.3784,nSAU)),", max var around growth \n",
       file=filename, append=TRUE)
   cat("sSigMax ,",genconst(rep(1e-04,nSAU)),", var of SigMax \n",
       file=filename, append=TRUE)
   cat("LML ,",genconst(rep(132,nSAU)),", initial legal minimum length \n",
       file=filename, append=TRUE)  # deprecated. now set in the control file
   cat("Wtb ,3.162046,3.162046,3.162046,3.162046,3.292563,3.162046,3.183611,3.162046 ,",
       ", weight-at-length exponent \n", file=filename, append=TRUE)
   cat("sWtb ,",genconst(rep(0.000148,nSAU)),", var of Wtb \n",
       file=filename, append=TRUE)
   cat("Wta , 6.4224e-05,5.62e-05,6.4224e-05,6.4224e-05,3.388e-05,6.4224e-05,5.64e-05,",
       "5.62e-05, intercept of power curve  \n",
       file=filename,append=TRUE)
   cat("sWta ,",genconst(rep(1e-08,nSAU)),", \n",file=filename,append=TRUE)
   cat("Me ,",genconst(rep(0.15,nSAU)),", Nat Mort, 0.2 maxage=23, 0.1 maxage=46 \n",
       file=filename, append=TRUE)
   cat("sMe ,",genconst(rep(0.0003,nSAU)),", var of Me\n",
       file=filename, append=TRUE)
   cat("AvRec ,245104.49443123,466856.350762581,260150.843272428,1040501.51466161,",
       "838542.867453192,1764761.04011726,1514858.94140123,581059.118386339 , LnR0 \n",
       file=filename, append=TRUE)
   cat("sAvRec ,",genconst(rep(0.00001,nSAU)),", \n",file=filename, append=TRUE)
   cat("defsteep ,",genconst(rep(0.7,nSAU)),", Beverton-Holt steepness \n",
       file=filename, append=TRUE)
   cat("sdefsteep ,",genconst(rep(0.0002,nSAU)),", \n",file=filename,
       append=TRUE)
   cat("L50C ,",genconst(rep(126.4222,nSAU)),", length at 50% emergent \n",
       file=filename, append=TRUE)
   cat("sL50C ,",genconst(rep(5e-04,nSAU)),", \n",file=filename, append=TRUE)
   cat("deltaC ,",genconst(rep(10.0,nSAU)),", length at 95% emergent \n",
       file=filename, append=TRUE)
   cat("sdeltaC ,",genconst(rep(1e-04,nSAU)),", \n",file=filename, append=TRUE)
   cat("MaxCEpar,0.374212895,0.237660281,0.45,0.475,0.45,0.45,0.375,0.3,",
       "max cpue t-hr \n",file=filename, append=TRUE) # deprecated
   cat("sMaxCEvar,",genconst(rep(0.02,nSAU)),", \n",file=filename, append=TRUE)
   cat("selL50p ,",genconst(rep(0.0,nSAU)),", L50 of selectivity \n",
       file=filename, append=TRUE)
   cat("selL95p ,5.37464213465647,5.55601436678328,3.85858046229972,3.2701437202985,",
       "2.86997433675897,4.10479572671887,4.61890731485338,3.90711525125475 ,",
       " L95 of selectivity \n",file=filename, append=TRUE)
   cat("SaMa,  -22.371,-22.371,-22.2536,-22.1666,-24.0026,-15.2194,-21.2821,-22.1641 ,",
       " matuirity logistic a par \n",file=filename, append=TRUE)
   cat("L50Mat,98.8992042,98.8992042,116.395209,116.7892518,122.9011777,",
       "112.2374631,121.1964692,105.8963211 ,  L50 for maturity b = -1/L50\n",
       file=filename, append=TRUE)
   cat("sL50Mat,",genconst(rep(0.00015,nSAU)),", \n",file=filename, append=TRUE)
   cat("lambda,",genconst(rep(0.75,nSAU)),", \n",file=filename, append=TRUE)
   cat("qest,4.7580120450086,2.24922282881994,5.17035161801573,1.63308027497998,",
       "1.27670203205689,0.777980219036715,0.609440955039425,1.45560892679734 , \n",
       file=filename, append=TRUE)
   cat("\n\n",file=filename, append=TRUE)
   cat("propREC  \n",file=filename, append=TRUE)
   cat("SAU, pop,propR \n",file=filename, append=TRUE)
   cat("6,1,0.742, \n",file=filename, append=TRUE)
   cat("6,2,0.094, \n",file=filename, append=TRUE)
   cat("6,3,0.163, \n",file=filename, append=TRUE)
   cat("7,4,0.203, \n",file=filename, append=TRUE)
   cat("7,5,0.590, \n",file=filename, append=TRUE)
   cat("7,6,0.207, \n",file=filename, append=TRUE)
   cat("8,7,0.041, \n",file=filename, append=TRUE)
   cat("8,8,0.330, \n",file=filename, append=TRUE)
   cat("8,9,0.090, \n",file=filename, append=TRUE)
   cat("8,10,0.400, \n",file=filename, append=TRUE)
   cat("8,11,0.139, \n",file=filename, append=TRUE)
   cat("9,12,0.074, \n",file=filename, append=TRUE)
   cat("9,13,0.084, \n",file=filename, append=TRUE)
   cat("9,14,0.275, \n",file=filename, append=TRUE)
   cat("9,15,0.085, \n",file=filename, append=TRUE)
   cat("9,16,0.187, \n",file=filename, append=TRUE)
   cat("9,17,0.250, \n",file=filename, append=TRUE)
   cat("9,18,0.045, \n",file=filename, append=TRUE)
   cat("10,19,0.104, \n",file=filename, append=TRUE)
   cat("10,20,0.080, \n",file=filename, append=TRUE)
   cat("10,21,0.065, \n",file=filename, append=TRUE)
   cat("10,22,0.100, \n",file=filename, append=TRUE)
   cat("10,23,0.130, \n",file=filename, append=TRUE)
   cat("10,24,0.070, \n",file=filename, append=TRUE)
   cat("10,25,0.129, \n",file=filename, append=TRUE)
   cat("10,26,0.184, \n",file=filename, append=TRUE)
   cat("10,27,0.138, \n",file=filename, append=TRUE)
   cat("11,28,0.050, \n",file=filename, append=TRUE)
   cat("11,29,0.164, \n",file=filename, append=TRUE)
   cat("11,30,0.108, \n",file=filename, append=TRUE)
   cat("11,31,0.152, \n",file=filename, append=TRUE)
   cat("11,32,0.173, \n",file=filename, append=TRUE)
   cat("11,33,0.199, \n",file=filename, append=TRUE)
   cat("11,34,0.017, \n",file=filename, append=TRUE)
   cat("11,35,0.052, \n",file=filename, append=TRUE)
   cat("11,36,0.085, \n",file=filename, append=TRUE)
   cat("12,37,0.055, \n",file=filename, append=TRUE)
   cat("12,38,0.015, \n",file=filename, append=TRUE)
   cat("12,39,0.060, \n",file=filename, append=TRUE)
   cat("12,40,0.080, \n",file=filename, append=TRUE)
   cat("12,41,0.040, \n",file=filename, append=TRUE)
   cat("12,42,0.025, \n",file=filename, append=TRUE)
   cat("12,43,0.137, \n",file=filename, append=TRUE)
   cat("12,44,0.368, \n",file=filename, append=TRUE)
   cat("12,45,0.080, \n",file=filename, append=TRUE)
   cat("12,46,0.065, \n",file=filename, append=TRUE)
   cat("12,47,0.045, \n",file=filename, append=TRUE)
   cat("12,48,0.030, \n",file=filename, append=TRUE)
   cat("13,49,0.121, \n",file=filename, append=TRUE)
   cat("13,50,0.084, \n",file=filename, append=TRUE)
   cat("13,51,0.124, \n",file=filename, append=TRUE)
   cat("13,52,0.164, \n",file=filename, append=TRUE)
   cat("13,53,0.215, \n",file=filename, append=TRUE)
   cat("13,54,0.126, \n",file=filename, append=TRUE)
   cat("13,55,0.062, \n",file=filename, append=TRUE)
   cat("13,56,0.104, \n",file=filename, append=TRUE)
   return(invisible(filename))
}  # end of datafileTemplate

# Utility functions used within parseFile, not exported


#' @title readctrlfile checks rundir contains the required csv files
#'
#' @description readctrlfile checks rundir contains the required csv
#'     files including the named control file, which then contains
#'     the names of the region data file, and the population data
#'     file. The run stops if any are not present or are misnamed.
#'
#' @param rundir the directory in which all files relating to a
#'     particular run are to be held.
#' @param infile default="control.csv", the filename of the control
#'     file present in rundir containing information regarding the
#'     run.
#' @param verbose Should progress comments be printed to console, default=TRUE
#'
#' @return the control list for the run
#' @export
#'
#' @examples
#' \dontrun{
#' # this has sinoce been modified and needs updating
#' rundir <- tempdir()
#' ctrlzonetemplate(rundir)
#' datafiletemplate(6,rundir,filename="zone1sau2pop6.csv")
#' ctrl <- readctrlfile(rundir)
#' ctrl
#' }
readctrlfile <- function(rundir,infile="control.csv",verbose=TRUE) {
   # rundir=rundir; infile=controlfile; verbose=verbose
   filenames <- dir(rundir)
   if (length(grep(infile,filenames)) != 1)
      stop(cat(infile," not found in ",rundir," \n"))
   filename <- filenametopath(rundir,infile)
   indat <- readLines(filename)   # reads the whole file as character strings
   begin <- grep("START",indat) + 1
   runlabel <- getStr(indat[begin],1)
   datafile <- getStr(indat[begin+1],1)
   bysau <- getsingleNum("bysau",indat)
   parfile <- getStr(indat[begin+3],1)
   optpars=NULL
   parsin=FALSE
   if ((!is.na(parfile)) & (nchar(parfile) > 4)) {
     paramfile <- filenametopath(rundir,parfile)
     optpars <- exp(read.csv(paramfile,header=TRUE,row.names=1))
     parsin <- TRUE
   }
   batch <- getsingleNum("batch",indat)
   reps <- getsingleNum("replicates",indat)
   withsigR <- getsingleNum("withsigR",indat)
   withsigB <- getsingleNum("withsigB",indat)
   withsigCE <- getsingleNum("withsigCE",indat)
   hyrs=40 # minimum to set up equilibrium; should this be altered?
   filenames2 <- dir(rundir)
   if (length(grep(datafile,filenames2)) != 1)
      stop("population data file not found \n")
   if (verbose) cat("All required files appear to be present \n")
   # Now read zone data
   nSAU <-  getsingleNum("nSAU",indat) # number of spatial management units
   begin <- grep("SAUpop",indat)
   SAUpop <-  getConst(indat[begin],nSAU) # how many populations per SAU
   numpop <- sum(SAUpop)
   SAUnames <- getStr(indat[begin+1],nSAU)
   initdepl <- getConst(indat[begin+2],nSAU)
   minc <-  getsingleNum("minc",indat) # minimum size class
   cw    <- getsingleNum("cw",indat) # class width
   Nclass <- getsingleNum("Nclass",indat) # number of classes
   midpts <- seq(minc,minc+((Nclass-1)*cw),2)
   larvdisp <- getsingleNum("larvdisp",indat)
   randomseed <- getsingleNum("randomseed",indat)
   randomseedP <- getsingleNum("randomseedP",indat)
   initLML <- getsingleNum("initLML",indat)
   projyrs <- getsingleNum("PROJECT",indat)
   rows <- c("runlabel","datafile","bysau","replicates","withsigR","withsigB",
             "withsigCE","nSAU","SAUpop","SAUnames","initdepl","minc","cw",
             "Nclass","larvdisp","randomseed","randomseedP","initLML",
             "Projections")
   controldf <- as.data.frame(matrix(0,nrow=length(rows),ncol=nSAU,
                                     dimnames=list(rows,SAUnames)))
   controldf[,1] <- c(runlabel,datafile,bysau,reps,withsigR,withsigB,
                      withsigCE,nSAU,NA,NA,NA,minc,cw,
                      Nclass,larvdisp,randomseed,randomseedP,initLML,
                      projyrs)
   controldf["SAUpop",] <- SAUpop
   controldf["SAUnames",] <- SAUnames
   controldf["initdepl",] <- initdepl
   projLML <- NULL
   pyrnames <- NULL
   HS <- NULL
   histCatch <- NULL
   histyr <- NULL
   hyrnames <- NULL
   histCE <- NULL
   yearCE <- NULL
   compdat=NULL
   if (projyrs > 0) {
      projLML <- numeric(projyrs)
      from <- grep("PROJLML",indat)
      for (i in 1:projyrs) {
         from <- from + 1
         projLML[i] <- getConst(indat[from],1)
      }
   } # end of projyrs if test
   catches <- getsingleNum("CATCHES",indat)
   if (catches > 0) {
      if (catches > hyrs) hyrs <- catches # fix historical years duration
      begin <- grep("CondYears",indat)[1]
      histCatch <- matrix(0,nrow=catches,ncol=nSAU)
      colnames(histCatch) <- SAUnames
      histyr <- matrix(0,nrow=hyrs,ncol=2)
      colnames(histyr) <- c("year","histLML")
      for (i in 1:hyrs) {
         begin <- begin + 1
         asnum <- as.numeric(unlist(strsplit(indat[begin],",")))
         histyr[i,] <- asnum[1:2]
         histCatch[i,] <- asnum[3:(nSAU+2)]
      }
      rownames(histCatch) <- histyr[,1]
      rownames(histyr) <- histyr[,1]
   } # end of catches loop
   yrce <- getsingleNum("CEYRS",indat)
   if (yrce > 0) {
      # yrce <- hyrs - startce + 1
      begin <- grep("CEYRS",indat)
      histCE <- matrix(NA,nrow=yrce,ncol=nSAU)
      yearCE <- numeric(yrce)
      colnames(histCE) <- SAUnames
      for (i in 1:yrce) {
          begin <- begin + 1
          cenum <- as.numeric(unlist(strsplit(indat[begin],",")))
          yearCE[i] <- cenum[1]
          histCE[i,] <- cenum[2:(nSAU+1)]
       }
      rownames(histCE) <- yearCE
      hyrnames <- as.numeric(histyr[,1])
      firstyear <- tail(hyrnames,1) + 1
      if (projyrs > 0) pyrnames <- firstyear:(firstyear + projyrs - 1)
   } # end of if(yrce == 0)
   sizecomp <- getsingleNum("SIZECOMP",indat)
   if (sizecomp > 0) {
      lffiles <- NULL
      locsizecomp <- grep("SIZECOMP",indat)
      if (sizecomp > 1) {
         compdat <- vector("list",sizecomp)
         for (i in 1:sizecomp) {
            lffilename <- removeEmpty(unlist(strsplit(indat[locsizecomp+i],",")))
            lffiles <- c(lffiles,lffilename)
            compdat[[i]] <- getLFdata(rundir,lffilename)
         }
       } else {
         lffilename <- removeEmpty(unlist(strsplit(indat[locsizecomp+1],",")))
         lffiles <- c(lffiles,lffilename)
         compdat <- getLFdata(rundir,lffilename)
      }
   }  # end of sizecomp loop
   recdevs <- matrix(-1,nrow=hyrs,ncol=nSAU,dimnames=list(hyrnames,SAUnames))
   if (parsin) {
     firstrec <- 6
     endrec <- length(optpars[,1])
     recd <- as.matrix(optpars[(firstrec:endrec),])
     recyrs <- as.numeric(gsub("d","",rownames(recd)))
     pickrows <- match(recyrs,hyrnames)
     recdevs[pickrows,] <- recd
   } else {
     rdevs <- getsingleNum("RECDEV",indat)
     if (rdevs > 0) {
        if (rdevs != hyrs)  #rdevs <- hyrs-4
           stop("rows of control file recdevs not equal to conditioning years \n")
        begin <- grep("RECDEV",indat) + 1
        for (i in 1:rdevs) {
           begin <- begin + 1
           devs <- as.numeric(unlist(strsplit(indat[begin],",")))
           recdevs[i,] <- devs[2:(nSAU+1)]
        }
     }
   }# end of recdev loop
   # make output objects
   condC <- list(histCatch=histCatch,histyr=histyr,
                 histCE=histCE,yearCE=yearCE,initdepl=initdepl,
                 compdat=compdat,recdevs=recdevs,parsin=parsin,optpars=optpars,
                 sizecomp=sizecomp,lffiles=lffiles,poprec=NULL)
   projC <- list(projLML=projLML,projyrs=projyrs,
                 Sel=NULL,SelWt=NULL,histCE=histCE)
   outctrl <- list(runlabel,datafile,infile,reps,randomseed,randomseedP,
                   withsigR,withsigB,withsigCE,catches,projyrs,bysau,rundir)
   names(outctrl) <- c("runlabel","datafile","controlfile","reps","randseed",
                       "randseedP","withsigR","withsigB","withsigCE",
                       "catches","projection","bysau","rundir")
   globals <- list(numpop=numpop, nSAU=nSAU, midpts=midpts,Nclass=Nclass,
                   reps=reps,hyrs=hyrs,pyrs=projyrs,hyrnames=hyrnames,
                   pyrnames=pyrnames,saunames=SAUnames,SAUpop=SAUpop,
                   larvdisp=larvdisp)
   totans <- list(SAUnames,SAUpop,minc,cw,larvdisp,randomseed,
                  initLML,condC,projC,globals,outctrl,catches,projyrs)
   names(totans) <- c("SAUnames","SAUpop","minc","cw","larvdisp","randomseed",
                     "initLML","condC","projC","globals","ctrl",
                     "catches","projyrs")
   return(totans)
} # end of readctrlzone

#' @title readdatafile reads in a matrix of data defining each population
#'
#' @description readdatafile expects a matrix of probability density
#'     function definitions that are used to define the populations
#'     used in the simulation. These constitute the definition of
#'     popdefs.
#'
#' @param numpop the total number of populations across the zone
#' @param indir directory in which to find the date file
#' @param infile character string with filename of the data file
#'
#' @return a matrix of values defining the PDFs used to define the
#'     properties of each population. The contents of popdefs
#' @export
#'
#' @examples
#' \dontrun{
#' data(zone1)
#' glb <- zone1$globals
#' glb
#' data(constants)
#' constants
#' ctrlfile <- "control.csv"
#' ctrl <- readctrlfile(glb$numpop,rundir,ctrlfile)
#' reg1 <- readzonefile(rundir,ctrl$zonefile)
#' popdefs <- readdatafile(reg1$globals,rundir,ctrl$datafile)
#' print(popdefs)
#' }
readdatafile <- function(numpop,indir,infile) {  # indir=rundir;infile="zone1sau2pop6.csv";numpop=6
   filename <- filenametopath(indir,infile)
   indat <- readLines(filename)   # reads the whole file as character strings
   begin <- grep("PDFs",indat)
   npar <- getConst(indat[begin],1)
   rows <- c("popnum","SAU","DLMax","sMaxDL","L50","sL50","L50inc","sL50inc","SigMax",
             "sSigMax","LML","Wtb","sWtb","Wtbtoa","sWtbtoa","Me","sMe",
             "AvRec","sAvRec","defsteep","sdefsteep","L50C","sL50C",
             "deltaC","sdeltaC","MaxCEpars","sMaxCEpars","selL50p",
             "selL95p","SaMa","L50Mat","sL50Mat")
   ans <- matrix(0,nrow=length(rows),ncol=numpop)
   begin <- begin + 1
   for (i in 1:npar) {
      ans[i,] <- getConst(indat[begin],numpop)
      begin <- begin + 1
   } # completed filling ans matrix
   rownames(ans) <- rows
   colnames(ans) <- ans["popnum",]
   return(ans)
} # end of readdatafile


#' @title readsaudatafile generates the constants matrix from sau data
#'
#' @description readsaudatafile uses data described at the SAU level to make
#'     the constants file, which is then used to generate the popdefs matrix
#'     containing the specific values. This change from explicitly defining
#'     each population has been implemented to simplify the conditioning of
#'     each operating model.
#'
#' @param rundir the directory in which the data file is to be found. This will
#'     usually be the rundir for the scenario run
#' @param infile the name of the specific datafile used.
#' @param optpar the optimum parameters from sizemod used to replace inputs for
#'     the main parameters estimated.
#'
#' @return a list of the constants matrix with values for each population and
#'     the original matrix of sau values from readsaudatafile
#' @export
#'
#' @examples
#' print("wait on suitable data sets")
#' # rundir=rundir; infile=ctrl$datafile;optpar=opar
readsaudatafile <- function(rundir,infile,optpar=NULL) {
   filename <- filenametopath(rundir,infile)
   indat <- readLines(filename)   # reads the whole file as character strings
   nsau <- getsingleNum("nsau",indat)
   saupop <- getConst(indat[grep("saupop",indat)],nsau)
   numpop <- sum(saupop)
   checkmsedata(intxt=indat,rundir=rundir,verbose=TRUE)
   txt <- indat[grep("saunames",indat)]
   saunames <- unlist(strsplit(txt,","))[2:(nsau+1)]
   rows <- c("DLMax","sMaxDL","L50","sL50","L50inc","sL50inc","SigMax",
             "sSigMax","LML","Wtb","sWtb","Wta","sWta","Me","sMe",
             "AvRec","sAvRec","defsteep","sdefsteep","L50C","sL50C",
             "deltaC","sdeltaC","MaxCEpars","sMaxCEpars","selL50p",
             "selL95p","SaMa","L50Mat","sL50Mat","lambda","qest")
 #  numrow <- length(rows)
   npar <- length(rows)
   ans <- matrix(0,nrow=npar,ncol=nsau)
   begin <- grep("PDFs",indat) + 1
   for (i in 1:npar) {
      ans[i,] <- getConst(indat[begin],nsau)
      begin <- begin + 1
   } # completed filling ans matrix
   rownames(ans) <- rows
   colnames(ans) <- saunames
   if (!is.null(optpar)) {
      sourcerows <- c("LnR0","MaxDL","L95","qest","seldelta")
      targetrow <- c("AvRec","DLMax","L50inc","qest","selL95p")
      numvar <- length(sourcerows)
      for (i in 1:numvar) {
        replace <- optpar[sourcerows[i],]
        if (sourcerows[i] == "L95")
            replace <- optpar[sourcerows[i],] - ans["L50",]
        ans[targetrow[i],] <- replace
      }
   }
   poprec <- matrix(0,nrow=numpop,ncol=3,
                    dimnames=list(1:numpop,c("sau","pop","prec")))
   begin <- grep("propREC",indat) + 2
   for (i in 1:numpop) {
      poprec[i,] <- getConst(indat[begin],nb=3,index=1)
      begin <- begin + 1
   }
   pop <- 1  # needed for the following counting
   sauindex <- numeric(numpop); names(sauindex) <- poprec[,"sau"]
   for (i in 1:nsau) {
      npop <- saupop[i]
      for (j in 1:npop) {
         sauindex[pop] <- i
         pop <- pop + 1
      }
   }
   outrows <- c("SAU",rows)
   numrow <- length(outrows)
   consts <- matrix(0,nrow=numrow,ncol=numpop,
                    dimnames=list(outrows,poprec[,"pop"]))
   consts["SAU",] <- poprec[,"sau"]
   for (index in 1:npar) { # index=17
      vect <- ans[rows[index],]
      if (rows[index] == "AvRec") {
         consts[rows[index],] <- log(vect[sauindex] * poprec[,"prec"])
      } else {
         consts[rows[index],] <- vect[sauindex]
      }
   }
   return(list(constants=consts,saudat=ans,poprec=poprec))
} # end of readsaudatafile

#' @title replaceVar replaces values of a variable in the input datafile
#'
#' @description replaceVar replaces the values of a variable in the input
#'     datafile. TAKE CARE, it overwrites the original file! However, it
#'     also saves the original file by adding an '_old' to the filename.
#'     Obviously if the function is used multiple times the original
#'     'original' file will be over-written on the second use.
#'     Alternatively, one could readin the datafile using readdateFile and
#'     then directly alter the condDat, e.g., in the case of a four block
#'     zone in which the AvRec vale is changed: condDat$constants["AvRec",]
#'     <- c(12.5,12.2,12.4,12.1)
#' @param infile the name, and path,  of the data file to be altered
#' @param invar the text name of the variable to be changed
#' @param newval the new value with which to replace the current values
#' @return The function over-writes the original file but saves the original
#'     by adding an '_old' to the end of the original filename.
#' @export
#' @examples
#' \dontrun{
#' filename <- datafiletemplate(numblock=1,filename="oneblock.csv")
#' replace(filename,"AvRec",15.75)
#' condDat <- readdatafile(filename)
#' print(round(condDat$constants,4))
#' }
replaceVar <- function(infile,invar,newval) {
  indat <- readLines(infile)   # reads the whole file as character strings
  filenamelen <- nchar(infile)
  oldfilename <- paste0(substr(infile,1,(filenamelen-4)),"_old.csv")
  write(indat,file=oldfilename) #,row.names=FALSE)
  begin <- grep("BLOCKNAMES",indat)
  numblk <- getConst(indat[begin],1) # number of blocks
  begin <- grep(invar,indat)
  if (length(begin) > 0) {
    begin <- begin[1]
  } else {
    stop("Error in 'replace' unknown variable used \n")
  }
  ans <- newval
  if (numblk > 1) for (i in 2:numblk) ans <- paste(ans,newval,sep=",")
  ans <- paste(invar,ans,",",sep=",")
  indat[begin] <- ans
  write(indat,file=infile) #,row.names=FALSE)
}  # end of replaceVar




#' @title rewritecontrolfile generates a revised control file
#'
#' @description rewritecontrolfile is used to generate a new control file after
#'     the option of reading in parameters from sizemod has been used. The
#'     rewritten file contains the revised parameter values.
#'
#' @param indir the directory into which the new file should be written. This
#'     would usually be the same as rundir, the scenario directory
#' @param zone1 a large object inside the output object from makeequilzone that
#'     contains SAUnames, SAUpop, minc, cw, larvdisp, randomseed, initLML,
#'     condC, projC, globals, ctrl, catches, and projyrs. Obviously there is
#'     some redundency.
#' @param controlfile name of the original control file
#'
#' @return nothing but it does write a file called 'oldcontrolfilename_new.csv'
#'     into rundir
#' @export
#'
#' @examples
#' \dontrun{
#'   data(zone1)
#'   tmpdir <- tempdir()
#'   rewritecontrolfile(tmpdir,zone1,"controlS21.csv")
#'   dir(tmpdir)
#' }
rewritecontrolfile <- function(indir,zone1,controlfile) { # indir=rundir; zone1=zone$zone1
  ctrl <- zone1$ctrl
  glb <- zone1$globals
  condC <- zone1$condC
  filen <- gsub(".csv","_new.csv",controlfile)
  filename <- filenametopath(indir,filen)
  cat("DESCRIPTION \n",file=filename,append=FALSE)
  cat("Control file containing details of a particular run. Modify the  \n",
      file=filename,append=TRUE)
  cat("contents and this description to suit your own situation.  \n",
      file=filename,append=TRUE)
  cat("\n\n",file=filename,append=TRUE)
  cat("START \n",file=filename,append=TRUE)
  cat("runlabel,",ctrl$runlabel,", label for particular run \n",
      file=filename,append=TRUE)
  cat("datafile,",ctrl$datafile,", name of saudefs file \n",
      file=filename,append=TRUE)
  cat("bysau,",ctrl$bysau,", 1=TRUE and 0=FALSE  \n",file=filename,append=TRUE)
  cat("parfile,NA, name of file containing optimum parameters from sizemod \n",
      file=filename,append=TRUE)
  cat("\n\n",file=filename,append=TRUE)
  cat("zoneCOAST \n",file=filename,append=TRUE)
  cat("replicates,",glb$reps,", number of replicates, usually 250, 500 or more  \n",
      file=filename, append=TRUE)
  cat("withsigR,",ctrl$withsigR,", recruitment variability eg 0.5 \n",
      file=filename, append=TRUE)
  cat("withsigB,",ctrl$withsigB,", process error on exploitable biomass \n",
      file=filename, append=TRUE)
  cat("withsigCE,",ctrl$withsigCE,", process error on cpue calculations  \n",
      file=filename, append=TRUE)
  cat("\n\n",file=filename,append=TRUE)
  cat("ZONE \n",file=filename,append=TRUE)
  cat("nSAU,",glb$nSAU,", number of spatial management units eg 2 \n",
      file=filename,append=TRUE)
  cat("SAUpop,",paste0(zone1$SAUpop,collapse=','),", number of populations per SAU in sequence \n",
      file=filename,append=TRUE)
  cat("SAUnames,",paste0(zone1$SAUnames,collapse=','),", labels for each SAU \n",
      file=filename,append=TRUE)   # possibly deprecated
  cat("initdepl,",paste0(condC$initdepl,collapse=','),", initial depletion levels for each SAU \n",
      file=filename,append=TRUE)  # Only the control file initdepl is used
  cat("\n\n",file=filename,append=TRUE)
  cat("SIZE \n",file=filename,append=TRUE)
  cat("minc,",zone1$minc,", centre of minimum size class \n",file=filename,append=TRUE)
  cat("cw,",zone1$cw,", class width mm \n",file=filename,append=TRUE)
  cat("Nclass,",glb$Nclass,", number of size classes \n",file=filename,append=TRUE)
  cat("\n",file=filename,append=TRUE)
  cat("RECRUIT \n",file=filename,append=TRUE)
  cat("larvdisp,",glb$larvdisp,", rate of larval dispersal eg 0.01 = 0.5 percent in either direction\n",
      file=filename,append=TRUE)
  cat("\n",file=filename,append=TRUE)
  cat("RANDOM \n",file=filename,append=TRUE)
  cat("randomseed,",ctrl$randseed,", for repeatability of population definitions set to 0 otherwise \n",
      file=filename,append=TRUE)
  cat("randomseedP,",ctrl$randseedP,", for repeatability of projections set to 0 otherwise \n",
      file=filename,append=TRUE)
  cat("\n",file=filename,append=TRUE)
  cat("initLML,",zone1$initLML,", the initial LML if no historical catches present \n",
      file=filename,append=TRUE)  # deprecated
  cat("\n",file=filename,append=TRUE)
  cat("PROJECT,",glb$pyrs,", number of projection years for each simulation \n",
      file=filename,append=TRUE)
  cat("\n\n",file=filename,append=TRUE)
  cat("PROJLML, need the same number as there are projection years \n",file=filename,
      append=TRUE) # ensure there are Nyrs lines
  pyrnames <- glb$pyrnames
  projLML <- zone1$projC$projLML
  for (i in 1:glb$pyrs) {
    cat(pyrnames[i],",",projLML[i],", \n",
        file=filename,append=TRUE)
  }
  cat("\n\n",file=filename,append=TRUE)
  cat("CATCHES,",glb$hyrs,", if > 1 then how many years in the histLML \n\n",
      file=filename,append=TRUE)
  cat("CondYears, LML,",paste0(zone1$SAUnames,collapse=','),", \n",
      file=filename,append=TRUE)
  catches <- condC$histCatch
  histyr <- condC$histyr
  nyr <- nrow(catches)
  for (i in 1:nyr) {   #  i=1
    cat(paste0(histyr[i,],collapse=','),",",paste0(catches[i,],collapse=','),
        ", \n",file=filename,append=TRUE)
  }
  cat("\n\n",file=filename,append=TRUE)
  cat(" ,",paste0(zone1$SAUnames,collapse=','),", \n",
      file=filename,append=TRUE)
  cpue <- condC$histCE
  cat("CEYRS ,",nrow(cpue),", \n",file=filename,append=TRUE)
  yrs <- as.numeric(rownames(cpue))
  nyr <- length(yrs)
  for (i in 1:nyr) { # i=1
    txt <- paste0(yrs[i],",",paste0(cpue[i,],collapse=','),",,")
    txt <- gsub("NA","",txt)
    cat(txt,", \n",file=filename,append=TRUE)
  }
  cat("\n\n",file=filename,append=TRUE)
  cat("SIZECOMP,", condC$sizecomp ,", 1 = 1 filename, 2 = 2 filenames, etc. \n",
      file=filename,append=TRUE)
  for (i in 1:condC$sizecomp)
    cat(condC$lffiles[i]," \n",file=filename,append=TRUE)
  cat("\n\n",file=filename,append=TRUE)
  recdevs <- condC$recdevs
  yrs <- as.numeric(rownames(recdevs))
  nyrs <- nrow(recdevs)
  cat("RECDEV,",nyrs," \n",file=filename,append=TRUE)
  cat("CondYears,",paste0(zone1$SAUnames,collapse=',')," \n",
      file=filename,append=TRUE)
  for (i in 1:nyrs) {
    cat(yrs[i],",",paste0(recdevs[i,],collapse=','),", \n",
        file=filename,append=TRUE)
  }
  cat("\n",file=filename,append=TRUE)
} # end of rewritecontrolfile




#' Title
#'
#' @param indir the directory into which the new file should be written. This
#'     would usually be the same as rundir, the scenario directory
#' @param zone1 a large object inside the output object from makeequilzone that
#'     contains SAUnames, SAUpop, minc, cw, larvdisp, randomseed, initLML,
#'     condC, projC, globals, ctrl, catches, and projyrs. Obviously there is
#'     some redundency.
#'
#' @return nothing but it does write a file called 'lffilename.csv'
#'     into rundir, which will overwrite the same if already there.
#' @export
#'
#' @examples
#' \dontrun{
#'   data(zone1)
#'   tmpdir <- tempdir()
#'   rewritecompdata(tmpdir,zone1)
#'   dir(tmpdir)
#' }
rewritecompdata <- function(indir,zone1) {
  lfs <- zone1$condC$compdat$lfs
  label <- dimnames(lfs)
  lengths <- as.numeric(label[[1]])
  nlen <- length(lengths)
  years <- as.numeric(label[[2]])
  nyr <- length(years)
  sau <- label[[3]]
  ans <- as.data.frame(matrix(0,nrow=(2*nlen),ncol=(nyr + 2)))
  colnames(ans) <- c("length","sau",years)
  ans[,"length"] <- c(lengths,lengths)
  ans[,"sau"] <- c(rep(sau[1],nlen),rep(sau[2],nlen))
  ans[1:nlen,3:(nyr+2)] <- lfs[,,1]
  ans[(nlen+1):(2*nlen),3:(nyr+2)] <- lfs[,,2]
  filename <- filenametopath(indir,zone1$condC$lffiles)
  write.table(ans,file=filename,sep=",",row.names=FALSE)
}


#' @title rewrwitedatafile generates a revised saudatafile
#'
#' @description rewritedatafile is used to generate a new saudata file after
#'     the option of reading in parameters from sizemod has been used. The
#'     rewritten file contains the revised parameter values.
#'
#' @param indir the directory into which the new file should be written. This
#'     would usually be the same as rundir, the scenario directory
#' @param zone1 a large object inside the output object from makeequilzone that
#'     contains SAUnames, SAUpop, minc, cw, larvdisp, randomseed, initLML,
#'     condC, projC, globals, ctrl, catches, and projyrs. Obviously there is
#'     some redundency.
#' @param saudat the main contents of the saudata file
#'
#' @return nothing but it does write a file called 'oldsaudatafile_new.csv' into
#'     the indir
#' @export
#'
#' @examples
#' \dontrun{
#'   data(zone1)
#'   data(saudat)
#'   tmpdir <- tempdir()
#'   rewritedatafile(tmpdir,zone1,saudat)
#'   dir(tmpdir)
#' }
rewritedatafile <- function(indir,zone1,saudat) {
  olddatfile <- zone1$ctrl$datafile
  filen <-gsub(".csv","_new.csv",olddatfile)
  glb <- zone1$globals
  filename <- filenametopath(indir,filen)
  cat("SAU definitions of Probability density function parameters for each variable \n",
      file=filename,append=FALSE)
  cat("Randomly asigned to each population except the proportion of recruitment \n",
      file=filename,append=TRUE)
  cat("which is literally allocated down in the popREC section. Average recruitment \n",
      file=filename,append=TRUE)
  cat("has a tiny sd which is to aid in simplifying conditioning. \n",
      file=filename,append=TRUE)
  cat("\n\n",file=filename,append=TRUE)
  cat("SPATIAL \n",file=filename,append=TRUE)
  cat("nsau,",glb$nSAU,", number of spatial management units  \n",file=filename,append=TRUE)
  cat("saupop,",paste0(glb$SAUpop,collapse=','),", number of populations per SAU in sequence  \n",
      file=filename,append=TRUE)
  cat("saunames,",paste0(glb$SAUnames,collapse=','),", labels for each SAU  \n",file=filename,append=TRUE)
  cat("\n\n",file=filename,append=TRUE)
  nvars <- nrow(saudat)
  vars <- rownames(saudat)
  cat("PDFs,",nvars,", \n",file=filename,append=TRUE)
  for (i in 1:nvars) {
    cat(vars[i],",",paste0(saudat[i,],collapse=','),", \n",file=filename,append=TRUE)
  }
  cat("\n\n",file=filename,append=TRUE)
  poprec <- zone1$condC$poprec
  npop <- nrow(poprec)
  cat("propREC,   \n",file=filename,append=TRUE)
  cat("SAU, pop, propR,  \n",file=filename,append=TRUE)
  for (i in 1:npop) {
    cat(paste0(poprec[i,],collapse=','),", \n",file=filename,append=TRUE)
  }
  cat("\n\n",file=filename,append=TRUE)
} # end of rewritesaudatafle.csv



#' @title writecompdata writes a size-composition file to indir
#'
#' @param indir the directory into which the new file should be written. This
#'     would usually be the same as rundir, the scenario directory
#' @param lfs the 3-dimensional array of length-composition data to be written
#'     to lf_WZ90-20.csv.
#' @param filen the filename of the size-comp data, default='lf_WZ90_20.csv'
#'
#' @return nothing but it does write a file called 'lf_WZ90-20.csv'
#'     into indir=rundir, which will overwrite the same if already there.
#' @export
#'
#' @examples
#' \dontrun{
#'   data(lfs)
#'   tmpdir <- tempdir()
#'   writecompdata(tmpdir,lfs)
#'   dir(tmpdir)
#' }
writecompdata <- function(indir,lfs,filen="lf_WZ90-20.csv") {
  label <- dimnames(lfs)
  lengths <- as.numeric(label[[1]])
  nlen <- length(lengths)
  years <- as.numeric(label[[2]])
  nyr <- length(years)
  sau <- label[[3]]
  nsau <- length(sau)
  ans <- as.data.frame(matrix(0,nrow=(nsau*nlen),ncol=(nyr + 2)))
  colnames(ans) <- c("length","sau",years)
  ans[,"length"] <- c(rep(lengths,nsau))
  tmp <- NULL
  for (i in 1:nsau) tmp <- c(tmp,rep(sau[i],nlen))
  ans[,"sau"] <- tmp
  ans[1:nlen,3:(nyr+2)] <- lfs[,,1]
  for (i in 2:nsau) {
    ans[((i-1)*nlen+1):(i*nlen),3:(nyr+2)] <- lfs[,,i]
  }
  filename <- filenametopath(indir,filen)
  write.table(ans,file=filename,sep=",",row.names=FALSE)
} # end of writecompdata







