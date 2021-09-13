
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LATEST UPDATE

-   2021-09-13 aMSE 0.0.0.075 Implemented the option of having an
    initial depletion as well as application of historical catches.

-   2021-09-10 aMSE 0.0.0.100 Many small, but significant, changes.
    Modified getdata which interacts with the \#\#\#HS.R file and added
    the zoneDP$NAS data input. Added a summary recdev plot to condition
    tab.

-   2021-09-03 aMSE 0.0.0.200 Lots of tidying of help pages, added
    numbersatsizeSAU. Added a fixed TasmanianHS.R file to the data-raw
    directory.

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
library(aMSE)
library(TasHS)
library(rutilsMH)
library(makehtml)
library(knitr)
# Obviously you should modify the rundir and datadir to suit your own setup
prefixdir <- "C:/Users/user/Dropbox/A_CodeUse/aMSEUse/scenarios/"
postfixdir <- "M15h75"
rundir <- paste0(prefixdir,postfixdir)
datadir <- rundir
alldirExists(rundir,datadir,verbose=TRUE)
#> rundir,  C:/Users/user/Dropbox/A_CodeUse/aMSEUse/scenarios/M15h75 :  exists
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
#               varyrs=7,startyr=48,cleanslate=TRUE,verbose=TRUE,
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
