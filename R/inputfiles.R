
# rutilsMH::listFunctions("C:/Users/User/Dropbox/A_code/aMSE/R/inputfiles.R")


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
   cat("runlabel, zoneX, label for particular run \n",
       file=filename,append=TRUE)
   cat("datafile, saudata_test.csv, name of saudefs file \n",
       file=filename,append=TRUE)
   cat("bysau, 1, 1=TRUE and 0=FALSE  \n",file=filename,append=TRUE)
   cat("\n",file=filename,append=TRUE)
   cat("zoneCOAST \n",file=filename,append=TRUE)
 #  cat("batch,  0, deprecated as batch jobs are run differently \n",
 #      file=filename, append=TRUE)
   cat("replicates,  100, number of replicates, usually 500 or more  \n",
       file=filename, append=TRUE)
   cat("withsigR,  0.3, recruitment variability eg 0.3 \n",
       file=filename, append=TRUE)
   cat("withsigB,  0.1, process error on exploitable biomass \n",
       file=filename, append=TRUE)
   cat("withsigCE, 0.15, process error on cpue calculations  \n",
       file=filename, append=TRUE)
   cat("\n",file=filename,append=TRUE)
   cat("ZONE \n",file=filename,append=TRUE)
   cat("nSAU, 8, number of spatial management units eg 2 \n",
       file=filename,append=TRUE)
   cat("SAUpop, 3, 3, 5, 7, 9, 9, 12, 8, number of populations per SAU in sequence \n",
       file=filename,append=TRUE)
   cat("SAUnames, sau6, sau7, sau8, sau9, sau10, sau11, sau12, sau13, labels for each SAU \n",
       file=filename,append=TRUE)   # possibly deprecated
   cat("initdepl, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, initial depletion levels for each SAU \n",
       file=filename,append=TRUE)  # possibly deprecated
   cat("\n",file=filename,append=TRUE)
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
   cat("initLML, 140, the initial LML for generating the unfished zone if no historical catches present \n",
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
   cat("SIZECOMP, 0, if 1 then a single filename should follow, if 2 then two filenames, etc. \n",
       file=filename,append=TRUE)
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
      for (i in 1:17)
         cat(as.character((1962+i)),",-1,-1,-1,-1,-1,-1,-1,-1 \n",file=filename,append=TRUE)
      cat("1980,1.20727,1.09742,0.93953,1.01431,0.6429,0.91117,0.50361,1.19751 \n",file=filename,append=TRUE)
      cat("1981,1.08316,0.14083,0.50044,1.54715,0.59112,0.50305,0.06646,0.64524 \n",file=filename,append=TRUE)
      cat("1982,0.922,0.91519,1.13053,0.99442,0.7901,0.9954,0.80744,1.06755 \n",file=filename,append=TRUE)
      cat("1983,1.01998,0.93356,1.12931,1.11423,0.88758,0.92337,0.78474,1.3606 \n",file=filename,append=TRUE)
      cat("1984,1.06304,0.91256,0.95035,0.96208,1.01281,0.91243,0.75447,1.12396 \n",file=filename,append=TRUE)
      cat("1985,1.05086,0.68669,0.94454,1.11627,0.98552,0.87839,0.75148,1.08797 \n",file=filename,append=TRUE)
      cat("1986,1.00913,0.38365,0.84443,1.02884,0.99212,1.05791,0.71612,1.0487 \n",file=filename,append=TRUE)
      cat("1987,1.16901,0.92909,0.85184,1.03493,0.97751,0.91665,0.64895,0.97533 \n",file=filename,append=TRUE)
      cat("1988,0.89468,0.86044,1.24648,1.06569,0.97573,1.15223,0.83098,0.78318 \n",file=filename,append=TRUE)
      cat("1989,0.8737,1.04136,1.29412,1.08022,0.90066,1.01103,1.02249,0.95434 \n",file=filename,append=TRUE)
      cat("1990,1.08596,1.34419,1.42049,0.82212,1.02519,1.07576,1.09688,0.93282 \n",file=filename,append=TRUE)
      cat("1991,1.06348,1.56493,2.22142,1.92953,1.61953,1.49053,1.52866,0.9503 \n",file=filename,append=TRUE)
      cat("1992,0.91822,1.26505,1.02269,0.9003,1.13051,1.04972,1.24638,0.90663 \n",file=filename,append=TRUE)
      cat("1993,0.88964,1.52089,0.86533,0.96544,1.09779,1.10389,1.13422,0.91746 \n",file=filename,append=TRUE)
      cat("1994,0.99061,1.07266,0.81647,1.09396,0.93548,1.13104,1.2707,0.65045 \n",file=filename,append=TRUE)
      cat("1995,0.88451,0.98911,0.39698,0.99352,1.07403,1.07916,1.03953,0.93436 \n",file=filename,append=TRUE)
      cat("1996,0.78078,0.93078,0.90728,1.06832,1.08105,1.04196,1.08,1.06868 \n",file=filename,append=TRUE)
      cat("1997,0.57762,0.99118,0.89432,1.13435,0.96783,1.0463,1.00321,0.83011 \n",file=filename,append=TRUE)
      cat("1998,0.80189,1.19164,0.06814,0.46693,1.54344,1.50513,0.94818,1.54093 \n",file=filename,append=TRUE)
      cat("1999,1.04149,1.08524,0.95322,0.91927,0.91088,1.09505,0.94709,1.40757 \n",file=filename,append=TRUE)
      cat("2000,1.48862,1.50403,1.17549,0.96855,1.15971,1.11267,0.80134,1.30331 \n",file=filename,append=TRUE)
      cat("2001,1.47781,1.14656,1.05932,1.00169,1.19329,1.09973,1.12002,1.57254 \n",file=filename,append=TRUE)
      cat("2002,1.46786,1.19021,1.09307,1.01745,1.03101,0.93205,1.11123,1.5447 \n",file=filename,append=TRUE)
      cat("2003,1.34601,1.33293,0.05729,1.09981,1.67906,1.44746,1.5489,1.29404 \n",file=filename,append=TRUE)
      cat("2004,0.97519,1.06773,0.94356,1.67039,0.7622,1.00882,1.18776,1.16923 \n",file=filename,append=TRUE)
      cat("2005,0.91323,1.09527,0.87703,1.05038,1.25863,0.99086,1.2271,1.36836 \n",file=filename,append=TRUE)
      cat("2006,0.71417,0.84869,0.93993,0.68625,1.13591,0.92534,1.16808,1.39069 \n",file=filename,append=TRUE)
      cat("2007,0.7424,0.61702,0.93948,0.84653,1.06592,0.98738,0.81591,0.23896 \n",file=filename,append=TRUE)
      cat("2008,0.57009,0.46971,0.05829,0.37503,0.11891,0.16579,1.17357,0.10712 \n",file=filename,append=TRUE)
      cat("2009,0.90396,0.96523,0.11664,0.54098,0.84044,0.92302,1.0334,0.6854 \n",file=filename,append=TRUE)
      cat("2010,0.59246,1.08064,0.25012,0.49833,0.63672,0.53775,0.70679,1.04342 \n",file=filename,append=TRUE)
      cat("2011,0.9251,1.08402,2.16096,0.4875,0.92798,0.96639,1.05925,1.21706 \n",file=filename,append=TRUE)
      cat("2012,1.02053,1.05273,0.28883,0.83263,0.90625,1.0863,1.28118,1.48546 \n",file=filename,append=TRUE)
      cat("2013,0.16033,0.60397,0.05743,0.33133,0.17506,0.16141,1.38297,0.78258 \n",file=filename,append=TRUE)
      cat("2015,1.00332,1.05593,0.43213,0.55848,0.98658,0.87003,0.99218,0.92184 \n",file=filename,append=TRUE)
      cat("2015,1.22426,1.63342,0.98686,0.87168,1.03108,0.92531,1.25532,1.3359 \n",file=filename,append=TRUE)
      cat("2016,1.22574,1.30988,1.17181,1.53304,1.16856,1.04363,1.4575,1.18581 \n",file=filename,append=TRUE)
      for (i in 55:58)
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
#' @param nSAU the number of SAU in the zone
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
#'  datafiletemplate(nSAU=8,yourdir,"saudata_test.csv")
#'  constants <- readsaudatafile(yourdir,"saudata_test.csv")
#'  str(constants,max.level=1)
#'  print(constants[,1:10])
#' }
datafiletemplate <- function(nSAU,indir,filename="saudata_test.csv") {
   genconst <- function(invect) {
      nlab <- length(invect)
      invect <- as.character(invect)
      ans <- invect[1]
      if (nlab > 1) for (i in 2:nlab) ans <- paste(ans,invect[i],sep=",")
      return(ans)
   }
   filename <- filenametopath(indir,filename)
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
   cat("initdepl,",genconst(rep(0.7,nSAU)),", initial depletion levels for each SAU  \n",
       file=filename,append=TRUE)  # possibly deprecated this will be clarified
   cat("\n",file=filename,append=TRUE)
   cat("PDFs,  30   \n",file=filename, append=TRUE)
   cat("MaxDL ,26,30,30,34,34,34,34,29, maximum growth increment \n",
       file=filename, append=TRUE)
   cat("sMaxDL ,",genconst(rep(2.0,nSAU)),", variation of MaxDL \n",
       file=filename, append=TRUE)
   cat("L50 ,",genconst(rep(135.276,nSAU)),", Length at 50% MaxDL \n",
       file=filename, append=TRUE)
   cat("sL50 ,",genconst(rep(2,nSAU)),", variation of L50 \n",
       file=filename, append=TRUE)
   cat("L50inc ,",genconst(rep(39,nSAU)),", L95 - L50 = delta =L50inc \n",
       file=filename, append=TRUE)
   cat("sL50inc ,",genconst(rep(1.5,nSAU)),", variation of L50inc \n",
       file=filename, append=TRUE)
   cat("SigMax  ,",genconst(rep(3.3784,nSAU)),", max var around growth \n",
       file=filename, append=TRUE)
   cat("sSigMax ,",genconst(rep(0.1,nSAU)),", var of SigMax \n",
       file=filename, append=TRUE)
   cat("LML ,",genconst(rep(127,nSAU)),", initial legal minimum length \n",
       file=filename, append=TRUE)  # deprecated. now set in the control file
   cat("Wtb ,",genconst(rep(3.161963,nSAU)),", weight-at-length exponent \n",
       file=filename, append=TRUE)
   cat("sWtb ,",genconst(rep(0.148,nSAU)),", var of Wtb \n",
       file=filename, append=TRUE)
   cat("Wtbtoa ,",genconst(rep(962.8098,nSAU)),
       ", intercept of power curve between Wtb and Wta \n",file=filename,
       append=TRUE)
   cat("sWtbtoa ,",genconst(rep(-14.3526,nSAU)),
       ", exponent of power curve between Wtb and Wta \n",file=filename,
       append=TRUE)
   cat("Me ,",genconst(rep(0.15,nSAU)),", Nat Mort, 0.2 maxage=23, 0.1 maxage=46 \n",
       file=filename, append=TRUE)
   cat("sMe ,",genconst(rep(0.003,nSAU)),", var of Me\n",
       file=filename, append=TRUE)
   cat("AvRec ,217500,410000,185000,905000,760000,1440000,1350000,355000,  LnR0 \n",
       file=filename, append=TRUE)
   cat("sAvRec ,",genconst(rep(0.0001,nSAU)),", \n",file=filename, append=TRUE)
   cat("defsteep ,",genconst(rep(0.75,nSAU)),", Beverton-Holt steepness \n",
       file=filename, append=TRUE)
   cat("sdefsteep ,",genconst(rep(0.02,nSAU)),", \n",file=filename,
       append=TRUE)
   cat("L50C ,",genconst(rep(126.4222,nSAU)),", length at 50% emergent \n",
       file=filename, append=TRUE)
   cat("sL50C ,",genconst(rep(0.5,nSAU)),", \n",file=filename, append=TRUE)
   cat("deltaC ,",genconst(rep(10.0,nSAU)),", length at 95% emergent \n",
       file=filename, append=TRUE)
   cat("sdeltaC ,",genconst(rep(0.1,nSAU)),", \n",file=filename, append=TRUE)
   cat("MaxCEpar,0.4,0.425,0.45,0.475,0.45,0.45,0.375,0.3, max cpue t-hr \n",
       file=filename, append=TRUE)
   cat("sMaxCEvar,",genconst(rep(0.02,nSAU)),", \n",file=filename, append=TRUE)
   cat("selL50p ,",genconst(rep(0.0,nSAU)),", L50 of selectivity \n",
       file=filename, append=TRUE)
   cat("selL95p ,",genconst(rep(1.5,nSAU)),", L95 of selectivity \n",
       file=filename, append=TRUE)
   cat("SaMa,",genconst(rep(-16,nSAU)),", matuirity logistic a par \n",
       file=filename, append=TRUE)
   cat("L50Mat,104,114,114,123,121,120,117,114, L50 for maturity b = -1/L50\n",
       file=filename, append=TRUE)
   cat("sL50Mat,",genconst(rep(2.5,nSAU)),", \n",file=filename, append=TRUE)
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


#' @title readctrlfile checks datadir contains the required csv files
#'
#' @description readctrlfile checks datadir contains the required csv
#'     files including the named control file, which then contains
#'     the names of the region data file, and the population data
#'     file. The run stops if any are not present or are misnamed.
#'
#' @param rundir the directory in which all files relating to a
#'     particular run are to be held. If datadir != rundir then the main data
#'     files are kept in datadir
#' @param infile default="control.csv", the filename of the control
#'     file present in datadir containing information regarding the
#'     run.
#' @param datadir default = rundir. If datadir != rundir then the data files
#'     for a series of scenarios are kept in this directory.
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
readctrlfile <- function(rundir,infile="control.csv",datadir=rundir,verbose=TRUE) {
   # rundir=rundir; infile="controlM15h75.csv"; datadir=datadir; verbose=verbose
   filenames <- dir(rundir)
   if (length(grep(infile,filenames)) != 1)
      stop(cat(infile," not found in ",rundir," \n"))
   filename <- filenametopath(rundir,infile)
   indat <- readLines(filename)   # reads the whole file as character strings
   begin <- grep("START",indat) + 1
   runlabel <- getStr(indat[begin],1)
   datafile <- getStr(indat[begin+1],1)
   bysau <- getsingleNum("bysau",indat)
   batch <- getsingleNum("batch",indat)
   reps <- getsingleNum("replicates",indat)
   withsigR <- getsingleNum("withsigR",indat)
   withsigB <- getsingleNum("withsigB",indat)
   withsigCE <- getsingleNum("withsigCE",indat)
   hyrs=40 # minimum to set up equilibrium; should this be altered?
   filenames2 <- dir(datadir)
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
   projLML <- NULL
   HS <- NULL
   histCatch <- NULL
   histyr <- NULL
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
      begin <- grep("CondYears",indat)
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
   startce <- getsingleNum("CEYRS",indat)
   yrce <- hyrs - startce + 1
   if (yrce == 0) {
      warning("CPUE calibration has no data")
   } else {
      begin <- grep("CEYRS",indat)
      histCE <- matrix(NA,nrow=yrce,ncol=nSAU)
      yearCE <- numeric(yrce) # of same length as nSAU
      colnames(histCE) <- SAUnames
      for (i in 1:yrce) {
         begin <- begin + 1
         cenum <- as.numeric(unlist(strsplit(indat[begin],",")))
         yearCE[i] <- cenum[1]
         histCE[i,] <- cenum[2:(nSAU+1)]
      }
      rownames(histCE) <- yearCE
   } # end of if(yrce == 0)
   hyrnames <- as.numeric(histyr[,1])
   firstyear <- tail(hyrnames,1) + 1
   pyrnames <- firstyear:(firstyear + projyrs - 1) # projection year names
   sizecomp <- getsingleNum("SIZECOMP",indat)
   if (sizecomp > 0) {
      lffiles <- NULL
      locsizecomp <- grep("SIZECOMP",indat)
      if (sizecomp > 1) {
         compdat <- vector("list",sizecomp)
         for (i in 1:sizecomp) {
            lffilename <- removeEmpty(unlist(strsplit(indat[locsizecomp+i],",")))
            compdat[[i]] <- getLFdata(datadir,lffilename)
         }
       } else {
         lffilename <- removeEmpty(unlist(strsplit(indat[locsizecomp+1],",")))
         compdat <- getLFdata(datadir,lffilename)
      }
   }  # end of sizecomp loop
   recdevs <- matrix(-1,nrow=hyrs,ncol=nSAU,dimnames=list(hyrnames,SAUnames))
   rdevs <- getsingleNum("RECDEV",indat)
   if (rdevs > 0) {
      if (rdevs != hyrs) rdevs <- hyrs-4
      #   stop("rows of recdevs not equal to conditioning years \n")
      begin <- grep("RECDEV",indat) + 1
      for (i in 1:rdevs) {
         begin <- begin + 1
         devs <- as.numeric(unlist(strsplit(indat[begin],",")))
         recdevs[i,] <- devs[2:(nSAU+1)]
      }
   } # end of recdev loop
   # make output objects
   condC <- list(histCatch=histCatch,histyr=histyr,
                 histCE=histCE,yearCE=yearCE,initdepl=initdepl,
                 compdat=compdat,recdevs=recdevs,Sel=NULL,SelWt=NULL)
   projC <- list(projLML=projLML,projyrs=projyrs,
                 Sel=NULL,SelWt=NULL,histCE=histCE)
   outctrl <- list(runlabel,datafile,reps,randomseed,randomseedP,
                   withsigR,withsigB,withsigCE,catches,projyrs,bysau)
   names(outctrl) <- c("runlabel","datafile","reps","randseed",
                       "randseedP","withsigR","withsigB","withsigCE",
                       "catches","projection","bysau")
   globals <- list(numpop=numpop, nSAU=nSAU, midpts=midpts,Nclass=Nclass,
                   reps=reps,hyrs=hyrs,pyrs=projyrs,hyrnames=hyrnames,
                   pyrnames=pyrnames,saunames=SAUnames,
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
#' ctrl <- readctrlfile(glb$numpop,datadir,ctrlfile)
#' reg1 <- readzonefile(datadir,ctrl$zonefile)
#' popdefs <- readdatafile(reg1$globals,datadir,ctrl$datafile)
#' print(popdefs)
#' }
readdatafile <- function(numpop,indir,infile) {  # indir=datadir;infile="zone1sau2pop6.csv";numpop=6
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
#' @param datadir the directory in which the data file is to be found. This will
#'     usually be the rundir for the scenario run
#' @param infile the name of the specific datafile used.
#'
#' @return the constants matrix with values for each population
#' @export
#'
#' @examples
#' print("wait on suitable data sets")
readsaudatafile <- function(datadir,infile) {  # rundir=rundir; infile=ctrl$datafile
   filename <- filenametopath(datadir,infile)
   indat <- readLines(filename)   # reads the whole file as character strings
   nsau <- getsingleNum("nsau",indat)
   saupop <- getConst(indat[grep("saupop",indat)],nsau)
   numpop <- sum(saupop)
   txt <- indat[grep("SAUnames",indat)]
   saunames <- unlist(strsplit(txt,","))[2:(nsau+1)]
#   saunames <- getConst(indat[grep("saunames",indat)],nsau)
   initdepl <- getConst(indat[grep("initdepl",indat)],nsau)
   begin <- grep("PDFs",indat)
   npar <- getConst(indat[begin],1)
   rows <- c("DLMax","sMaxDL","L50","sL50","L50inc","sL50inc","SigMax",
             "sSigMax","LML","Wtb","sWtb","Wtbtoa","sWtbtoa","Me","sMe",
             "AvRec","sAvRec","defsteep","sdefsteep","L50C","sL50C",
             "deltaC","sdeltaC","MaxCEpars","sMaxCEpars","selL50p",
             "selL95p","SaMa","L50Mat","sL50Mat")
   numrow <- length(rows)
   ans <- matrix(0,nrow=numrow,ncol=nsau)
   begin <- begin + 1
   for (i in 1:npar) {
      ans[i,] <- getConst(indat[begin],nsau)
      begin <- begin + 1
   } # completed filling ans matrix
   rownames(ans) <- rows
   colnames(ans) <- saunames
   poprec <- matrix(0,nrow=numpop,ncol=3,
                    dimnames=list(1:numpop,c("sau","pop","prec")))
   begin <- grep("propREC",indat) + 2
   for (i in 1:numpop) {
      poprec[i,] <- getConst(indat[begin],nb=3,index=1)
      begin <- begin + 1
   }
   pop <- 1
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
   for (index in 1:npar) { # index=16
      vect <- ans[rows[index],]
      if (rows[index] == "AvRec") {
         consts[rows[index],] <- log(vect[sauindex] * poprec[,"prec"])
      } else {
         consts[rows[index],] <- vect[sauindex]
      }
   }
   return(consts)
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


