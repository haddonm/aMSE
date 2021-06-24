
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LATEST UPDATE

  - 2021-06-24 aMSE 0.0.0.2100 Big jump in number as some larger
    changes. I have removed a bunch of deprecated functions (which are
    now in ‘deprecated.R’ in the data-raw directory). I have tidied many
    other functions, and have modified listall files to aid in the
    auto-documentation of the package.

# aMSE

<!-- badges: start -->

<!-- badges: end -->

This packages up a new Abalone Management Strategy Evaluation framework.
It has a novel structure that is based around the spatial ideas of
populations within spatial management units (SAUs), within a zone:

  - zone - highest geogrpahical level. Is simply the totality of the
    spatial management units.

  - SAU - spatial assessment units. In Tasmania these would currently be
    the classical statistical blocks.

  - population - literally a population. These are the active components
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
larval dispersal rate is \> 0.0.

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
#> Warning: package 'knitr' was built under R version 4.0.5
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
# run the scenario ----------------------------------------------------
#source(paste0(ddir,"aMSEUse/scenarios/MSE_Source.R"))
out <- do_MSE(rundir,controlfile,datadir,hsargs=hsargs,
              hcrfun=mcdahcr,sampleCE=tasCPUE,sampleFIS=tasFIS,sampleNaS=tasNaS,
              getdata=tasdata,calcpopC=calcexpectpopC,varyrs=7,startyr=32,
              cleanslate=TRUE,verbose=TRUE,doproject=TRUE,ndiagprojs=3)
#> All required files appear to be present 
#> Files read, now making zone 
#> Now estimating population productivity 
#> matureB Stable 
#> exploitB Stable 
#> recruitment varies 
#> spawning depletion Stable 
#> Time difference of 41.31817 secs
#> Conditioning on the Fishery data
#> Doing the projections
# make results webpage ---------------------------------------------------------
replist <- list(starttime=as.character(out$starttime),
                endtime=as.character(out$projtime))
glb <- out$glb
runnotes <- paste0(out$ctrl$runlabel,":  RunTime = ",out$tottime,
                   "  replicates = ",glb$reps,",   years projected = ",glb$pyrs,
                   "  Populations = ",glb$numpop," and SAU = ",glb$nSAU,
                   "  Randomseed for conditioning = ",out$ctrl$randseed)

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
str(out$zoneDP,max.level=1)
#> List of 17
#>  $ SAU     : num [1:56] 6 6 6 7 7 7 8 8 8 8 ...
#>  $ matureB : num [1:30, 1:56, 1:100] 139 139 143 145 153 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ exploitB: num [1:30, 1:56, 1:100] 130 113 105 108 111 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ midyexpB: num [1:30, 1:56, 1:100] 152 134 127 133 133 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ catch   : num [1:30, 1:56, 1:100] 11.7 12.3 13.3 16.2 13.3 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ acatch  : num [1:30, 1:8, 1:100] 13.6 16.3 18 18.8 18.8 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ harvestR: num [1:30, 1:56, 1:100] 0.0897 0.1094 0.126 0.15 0.1198 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ cpue    : num [1:30, 1:56, 1:100] 174 152 143 148 150 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ cesau   : num [1:30, 1:8, 1:100] 170 150 145 148 149 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ catsau  : num [1:30, 1:8, 1:100] 16 16.9 18.4 23.1 18.7 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ recruit : num [1:30, 1:56, 1:100] 167429 195829 168013 103647 142860 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ deplsB  : num [1:30, 1:56, 1:100] 0.433 0.433 0.445 0.45 0.474 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ depleB  : num [1:30, 1:56, 1:100] 0.449 0.389 0.364 0.372 0.384 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ catchN  : num [1:105, 1:30, 1:56, 1:100] 2.35e-131 1.18e-127 4.39e-124 1.18e-120 2.30e-117 ...
#>   ..- attr(*, "dimnames")=List of 4
#>  $ Nt      : num [1:105, 1:30, 1:56, 1:100] 1.67e+05 2.21e-06 9.08e-05 2.70e-03 5.83e-02 ...
#>   ..- attr(*, "dimnames")=List of 4
#>  $ NumNe   : num [1:105, 1:30, 1:56, 1:100] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..- attr(*, "dimnames")=List of 4
#>  $ TAC     : num [1:30, 1:100] 703 517 417 425 433 ...
#>   ..- attr(*, "dimnames")=List of 2
```

To see the structure of the dynamic object generated by the projections.

See the vignette Running\_aMSE.Rmd for a more detailed example.

See the New.Rd for recent developments
