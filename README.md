
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LATEST UPDATE

-   2021-05-07 0.0.0.2500 Some relatively large changes this time. I
    have pulled most of the code used to run the MSE into a separate
    source file ‘MSE\_source.R’. You can find all required files in
    abalone MSE\_files. I have moved TasmanianHS.R into tasdata as it is
    currently common to different uses. I will endeavour to generalize
    the contents of MSE\_source.R to turn it into a function that can be
    included in the package. Then any changes made there will flow
    naturally to other instances of each scenario file. Also, now there
    is a minimu requirement of 2 SAU but each SAU can have a minimum of
    1 population.

-   2021-04-20 0.0.0.2900 Continuing process of modifying the core
    functions to allow for the numbers-at-size prior to fishing and the
    midyear-exploitB (midyexpB) to be included in teh zoneD object.
    doprojection still requires fursther modification to account for
    applying to chosen HCR and HS to the completed conditioned zoneDD to
    provide the first year of expected catches in the projections (which
    are to be defined by teh conditoned data prior to projections
    beginning).

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
Tasmania MCDA using only 100 replicates (that bit takes about 15 seconds
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
alldirExists(rundir,datadir)
#> rundir,  c:/Users/Malcolm/DropBox/A_codeUse/aMSEUse/scenarios/HS652510 :  exists  
#> datadir,  c:/Users/Malcolm/DropBox/A_codeUse/aMSEUse/scenarios/tasdata :  exists
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
source(paste0(ddir,"aMSEUse/scenarios/MSE_Source.R"))
#> All required files appear to be present 
#> Files read, now making zone 
#> Now estimating population productivity 
#> matureB Stable 
#> exploitB Stable 
#> recruitment Stable 
#> spawning depletion Stable 
#> Time difference of 21.79685 secs
#> Conditioning on the Fishery data
#> Doing the projections

# make results webpage ---------------------------------------------------------
replist <- list(starttime=as.character(starttime),endtime=as.character(projtime))

# Unhash the make_html command to generate the results webspage
# make_html(
#   replist = replist,
#   rundir = rundir,
#   width = 500,
#   openfile = TRUE,
#   runnotes = NULL,
#   verbose = FALSE,
#   packagename = "aMSE",
#   htmlname = "aMSE"
# )
```

After running the whole, even if you do not generate the results
webpage, you could also try:

``` r
str(zoneDP,max.level=1)
#> List of 17
#>  $ SAU     : num [1:28] 6 6 6 7 7 7 8 8 9 9 ...
#>  $ matureB : num [1:30, 1:28, 1:50] 29.5 31.9 33.3 32.9 31.7 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ exploitB: num [1:30, 1:28, 1:50] 27.3 23.4 25.2 27.1 26.5 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ midyexpB: num [1:30, 1:28, 1:50] 32.2 28.9 31.1 32.9 32.4 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ catch   : num [1:30, 1:28, 1:50] 2.8 3.62 3.88 3.68 3.83 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ acatch  : num [1:30, 1:8, 1:50] 13.6 15 15 15.7 16.5 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ harvestR: num [1:30, 1:28, 1:50] 0.103 0.154 0.154 0.136 0.144 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ cpue    : num [1:30, 1:28, 1:50] 194 171 184 196 192 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ cesau   : num [1:30, 1:8, 1:50] 148 127 133 151 172 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ catsau  : num [1:30, 1:8, 1:50] 16 14.8 15.3 13.9 15.7 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ recruit : num [1:30, 1:28, 1:50] 21244 61309 70590 28067 24583 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ deplsB  : num [1:30, 1:28, 1:50] 0.42 0.454 0.473 0.468 0.451 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ depleB  : num [1:30, 1:28, 1:50] 0.421 0.362 0.389 0.418 0.409 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ catchN  : num [1:105, 1:30, 1:28, 1:50] 3.76e-132 1.90e-128 7.03e-125 1.89e-121 3.69e-118 ...
#>   ..- attr(*, "dimnames")=List of 4
#>  $ Nt      : num [1:105, 1:30, 1:28, 1:50] 2.12e+04 3.14e-07 1.29e-05 3.82e-04 8.25e-03 ...
#>   ..- attr(*, "dimnames")=List of 4
#>  $ NumNe   : num [1:105, 1:30, 1:28, 1:50] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..- attr(*, "dimnames")=List of 4
#>  $ TAC     : num [1:30, 1:50] 703 538 411 409 420 ...
#>   ..- attr(*, "dimnames")=List of 2
```

To see the structure of the dynamic object generated by the projections.

See the vignette Running\_aMSE.Rmd for a more detailed example.

See the New.Rd for recent developments
