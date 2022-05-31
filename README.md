
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LATEST UPDATE

-   2022-05-31 aMSE 0.0.4 Added functions that assist in conditioning
    the operating model by optimizing the fit to CPUE and
    size-composition data after opimizing the AvRec values for each SAU.
    No EWN.

-   2021-11-29 aMSE 0.0.1 Removed datadir. Started implementing
    hyper-stability while keeping the interface the same, so that the
    results from sizemod can be used.

-   2021-11-23 aMSE 0.0.0.1 Added functions to simplify the optimization
    of AvRec relative to the CPUE data after initial parameter
    estimation for each SAU using the sizemod package. Also added a
    comparison of size-composition data to the do\_condition function.
    Next step will be to include the optimization of the recdevs using
    both CPUE and size-comp data.

-   2021-11-21 aMSE 0.0.0.020 So as to use any observed size-composition
    data in the conditioning, a number of functions have been added to
    convert the predicted numbers-at-size in the catch by population
    into NAS x SAU. These can now be plotted against the observed NAS
    for each SAU and this can be either to the console or to a new tab
    in the webpage called ‘predictedcatchN’. This work derived from the
    progress made with the sizemod package used to fit size-based
    assessment models to individual SAU. Because the MSE sub-divides
    each SAU into a number of populations each with somewhat different
    MaxDL, L95, etc, even after fitting the models, when using that data
    it is then necessary to refit the AvRec values, which, because of
    the sub-division into populations each having different biology can
    lead to the MSE AvRec being either smaller or larger than the
    estimated R0 for a SAU. Similar arguments apply to teh recruitment
    residuals but generally the pattern provides a reasonable fit.

-   2021-10-17 aMSE 0.0.0.010 As cleanslate, within setuphtml() has been
    deprecated in the package makehtml, I have removed reference to it
    throughout aMSE. Now, to clean a directory of aMSE files use
    cleanrundir.

# aMSE

<!-- badges: start -->
<!-- badges: end -->

This packages up a new Abalone Management Strategy Evaluation framework.
It has a novel structure that is based around the spatial ideas of
populations within spatial management units (SAUs), within a zone:

-   zone - highest geographical level. Is simply the totality of the
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
options("show.signif.stars"=FALSE,
        "stringsAsFactors"=FALSE,
        "max.print"=50000,
        "width"=240)
# declare libraries ------------------------------------------------------------
suppressPackageStartupMessages({
  library(aMSE)
  library(TasHS)
  library(hutils)
  library(hplot)
  library(makehtml)
  library(knitr)
})
#> Warning: package 'knitr' was built under R version 4.1.3
dbdir <- getDBdir()
# Obviously you should modify the rundir and datadir to suit your own setup
prefixdir <- paste0(dbdir,"A_CodeUse/aMSEUse/scenarios/")
postfixdir <- "M15h75"
rundir <- paste0(prefixdir,postfixdir)
datadir <- rundir
alldirExists(rundir,datadir,verbose=TRUE)
#> rundir,  C:/Users/Malcolm/Dropbox/A_CodeUse/aMSEUse/scenarios/M15h75 :  exists
controlfile <- "controlM15h75.csv"
# equilibrium zone -------------------------------------------------------------
# You now need to ensure that there is, at least, a control.csv, and a 
# constantsdata.csv file in the data directory plus some other data .csv files
# depending on how conditioned you want the model to be. Templates for the
# correct format can be produced using ctrlfiletemplate(), datafiletemplate().
# 
# Of course, usually one would use data files, control.csv and a saudata.csv, which
# is listed as the datafile within the control.csv. These must be stored in 
# rundir. There are example files in the DropBox folder:
# .\..\Dropbox\National abalone MSE\aMSE_files\scenarios\examplerun. Read the 
# 'Using_aMSE' section in the Documentation .docx file in the documentation 
# directry  Then, assuming you have the very latest version of aMSE = 500, 
# the code in the documentation file (or in the 'run_aMSE.R' file within the
# examplerun subdirectory should work. 

# HarvestStrategy.R should also be in rundir
source(paste0(rundir,"/TasHS1_Tas.R"))
# run the scenario --------------------- Obviously unhash this to make it work
#
# out <- do_MSE(rundir,controlfile,datadir,hsargs=hsargs,hcrfun=mcdahcr,
#               sampleCE=tasCPUE,sampleFIS=tasFIS,sampleNaS=tasNaS,
#               getdata=tasdata,calcpopC=calcexpectpopC,makeouthcr=makeouthcr,
#               varyrs=7,startyr=48,verbose=TRUE,
#               ndiagprojs=4,savesauout=TRUE)
# 
# makeoutput(out,rundir,datadir,postfixdir,controlfile,openfile=TRUE,verbose=FALSE)
```

After running the whole, even if you do not generate the results
webpage, you could also try:

``` r
#str(out$zoneDP,max.level=1)
#str1(out$hcrout)
```

To see the structure of the dynamic object generated by the projections.

See the vignette Running\_aMSE.Rmd for a more detailed example.

See the New.Rd for recent developments
