
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LATEST UPDATE

-   2021\_07-04 aMSE 0.0.0.1500 Another big change. Functions added to
    aid conditioning, and aspects of the projections streamlined. No
    errors, warnings, or notes

-   2021-07-04 aMSE 0.0.0.1900 More big changes. I have removed the
    numbers-at-size large matrices from the zoneDP object into a NAS
    object. zoneDP is now only 12Mb instead of 450Mb and so can now be
    stored conveniently for each run. In addition, fixed recruitment
    deviates have now been implemented (although the default is to not
    use them - all are set to -1). But, expreimentally, in Tasmania,
    setting the deviate in 1991 = 2.0 meaning that the base level of
    recruitment off the recruitment curve is doubled, improves the
    relationship between predicted cpue and observed cpue so this holds
    great promise for improved conditioning.

# aMSE

<!-- badges: start -->
<!-- badges: end -->

This packages up a new Abalone Management Strategy Evaluation framework.
It has a novel structure that is based around the spatial ideas of
populations within spatial management units (SAUs), within a zone:

-   zone - highest geogrpahical level. Is simply the totality of the
    spatial management units.

-   SAU - spatial assessment units. In Tasmania these would currently be
    the classical statistical blocks.

-   population - literally a population. These are the active components
    within the simulation framework. The dynamics of the simulation are
    based around the populations, although, with positive larval
    dispersal (the default) there is some dependency of neighbouring
    populations on each other.

This conceptual structure has implications for how the simulations are
conditioned during the management strategy testing. While it is still
the case that the general properties of each SAU are used to set the
scene for each SAU’s component populations, the advent of sufficient GPS
logger data now allows a more formal definition of the productivity of
each population. Currently, the modelling assumes that each population
is linearly arranged around the coastline, and this will allow, where
available, specific catch/yield data to be used to define the bounds of
each population within an SAU.

An important change from previous designs is that larval dispersal is
implemented explicitly rather than implicitly being held constant. This
alters the dynamics so that analytical equilibrium methods no longer
work and we need to resort to iterative approaches to equilibrium if the
larval dispersal rate is &gt; 0.0.

Of course, adding such a component increases the number of options that
may need to be explored, but currently we envisage including a very low,
a middle range, and a relatively high level of larval dispersal, to
examine its implications.

## Installation

Once the MSE framework is made public you can install the development
version from [GitHub](https://github.com/haddonm/aMSE) with:

``` r
if (!require(devtools)){install.packages("devtools")} 

devtools::install_github("https://github.com/haddonm/aMSE",build_vignettes = TRUE)
```

Alternatively, while the development version remains private, you can
generate a branch that you can work on by cloning the repository, which,
again, can be done very simply within RStudio. Open the New Project
option in the project dialog at the top right of the RStudio screen and
selection Version Control, then use ‘<https://github.com/haddonm/aMSE>’
in the top box, identify where you want the new directory put, and press
return. ALternatively, you could download the zip file from inside the
‘code’ button and establish an R project from that.

It would be a good idea to read Hadley Wickham’s draft chapter on Git
and GitHub at <https://r-pkgs.org/index.html>.

## Example

This is an example which illustrates the generation of an initial
equilibrium, which then goes on to apply an early version of the
Tasmania MCDA using only 50 replicates (that bit takes about 15 seconds
on my computer, hence 2.5 minutes for 1000). It uses built in data-sets
but usually you would read in a control file, which would contain the
name of the biological datafile describing each population.

``` r
# a constant TAC example
starttime <- (Sys.time())
library(aMSE)
library(rutilsMH)
library(makehtml)
library(knitr)
# OBVIOUSLY you should modify the rundir and datadir to suit your own computer
if (dir.exists("c:/Users/User/DropBox")) {
  ddir <- "c:/Users/User/DropBox/A_codeUse/"
} else {
  ddir <- "c:/Users/Malcolm/DropBox/A_codeUse/"
}
doproject <- TRUE  # change to FALSE if only conditioning is required
verbose <- TRUE
rundir <- paste0(ddir,"aMSEUse/scenarios/HS652510")
datadir <- paste0(ddir,"aMSEUse/scenarios/tasdata")
alldirExists(rundir,datadir,verbose=TRUE)
#> rundir,  c:/Users/User/DropBox/A_codeUse/aMSEUse/scenarios/HS652510 :  exists  
#> datadir,  c:/Users/User/DropBox/A_codeUse/aMSEUse/scenarios/tasdata :  exists
# equilibrium zone -------------------------------------------------------------
# You now need to ensure that there is, at least, a control.csv, and a 
# constantsdata.csv file in the data directory plus some other data .csv files
# depending on how conditioned you want the model to be. Templates for the
# correct format can be produced using ctrlfiletemplate(), datafiletemplate().
# 
# Of course, usually one would use data files, control.csv and a zone.csv, which
# is listed as the datafile within the control.csv. These must be stored in 
# rundir. There are example control and data files in the DropBox folder:
# C:\Users\User\Dropbox\National abalone MSE\aMSE_files. Copy the control2.csv
# and zonewest.csv into you rundir. Then, assuming yo have the very latest
# version of aMSE = 5300, the following code should work. Recently both aMSE and
# makehtml have been altered so that all references to resdir have been changed
# to rundir, so both packages will need updating. See their respective GitHub
# readme pages for details.

# TasmanianHS.R should be in rundir, but during development is in data-raw
source(paste0(datadir,"/TasmanianHS.R"))
controlfile <- "controlsau.csv"
# run the scenario --------------------- Obviously unhash this to make it work
#out <- do_MSE(rundir,controlfile,datadir,hsargs=hsargs,
#              hcrfun=mcdahcr,sampleCE=tasCPUE,sampleFIS=tasFIS,sampleNaS=tasNaS,
#              getdata=tasdata,calcpopC=calcexpectpopC,varyrs=7,startyr=32,
#              cleanslate=TRUE,verbose=TRUE,doproject=TRUE,ndiagprojs=3)
# make results webpage ---------------------------------------------------------
# replist <- list(starttime=as.character(out$starttime),
#                 endtime=as.character(out$projtime))
# glb <- out$glb
# runnotes <- paste0(out$ctrl$runlabel,":  RunTime = ",out$tottime,
#                    "  replicates = ",glb$reps,",   years projected = ",glb$pyrs,
#                    "  Populations = ",glb$numpop," and SAU = ",glb$nSAU,
#                    "  Randomseed for conditioning = ",out$ctrl$randseed)
# Unhash the make_html command to generate the results webspage
# make_html(
#   replist = replist,
#   rundir = rundir,
#   width = 500,
#   openfile = TRUE,
#   runnotes = NULL,
#   verbose = FALSE,
#   packagename = "aMSE",
#   htmlname = "hs652510"
# )
```

After running the whole, even if you do not generate the results
webpage, you could also try:

``` r
#str(out$zoneDP,max.level=1)
```

To see the structure of the dynamic object generated by the projections.

See the vignette Running\_aMSE.Rmd for a more detailed example.

See the New.Rd for recent developments
