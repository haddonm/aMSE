
<!-- README.md is generated from README.Rmd. Please edit that file -->

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

You can install the development version from
[GitHub](https://github.com/haddonm/aMSE) with:

``` r
if (!require(devtools)){install.packages("devtools")} 

devtools::install_github("https://github.com/haddonm/aMSE",build_vignettes = TRUE)
```

Alternatively, you can generate a branch that you can work on by cloning
the repository, which, again, can be done very simply within RStudio.
Open the New Project option in the project dialog at the top right of
the RStudio screen and selection Version Control, then use
‘<https://github.com/haddonm/aMSE>’ in the top box, identify where you
want the new directory put, and press return.

It would be a good idea to read Hadley Wickham’s draft chapter on Git
and GitHub at <https://r-pkgs.org/index.html>.

## Example

This is a basic example which illustrates the generation of an initial
equilibrium. It uses built in data-sets but usually you would generate
zoneC (constants) and zoneD (the dynamic components) from the
conditioning data:

``` r
# a minimal example
starttime <- as.character(Sys.time())
library(aMSE)
library(rutilsMH)
library(makehtml)
# Obviously you should modify the resdir to suit your own computer
resdir <- "./../../A_code/aMSEUse/out/testrun"
dirExists(resdir,make=TRUE,verbose=TRUE)
#> ./../../A_code/aMSEUse/out/testrun :  exists
# The following function calls will ensure that there is a control.csv,
# a zone1sm\au2pop6.csv
# and region1.csv file in the data directory
ctrlfiletemplate(resdir)
zonefiletemplate(resdir)
datafiletemplate(6,resdir,filename="zone1sau2pop6.csv")

ctrl <- checkresdir(resdir)
#> All required files appear to be present
runname <- ctrl$runlabel
zone1 <- readzonefile(resdir,ctrl$zonefile)
glb <- zone1$globals     # glb without the movement matrix
constants <- readdatafile(glb$numpop,resdir,ctrl$datafile)

out <- setupzone(constants,zone1) # make operating model
zoneC <- out$zoneC
zoneD <- out$zoneD
glb <- out$glb        # glb now has the movement matrix
product <- out$product     # important bits usually saved in resdir
          # did the larval dispersal level disturb the equilibrium?
regDe <- testequil(zoneC,zoneD,glb)
#> [1] matureB Stable
#> [1] exploitB Stable
#> [1] recruitment Stable
#> [1] spawning depletion Stable

resfile <- setuphtml(resdir,runname)# prepare to save and log results

plotproductivity(resdir,runname,product,glb)
biology_plots(resdir, runname, glb, zoneC)
numbersatsize(resdir, runname, glb, zoneD)

endtime <- as.character(Sys.time())

reportlist <- list(runname=runname,
                   starttime=starttime,endtime=endtime,
                   zoneC=zoneC, zoneD=zoneD, product=product,
                   glb=glb,constants=constants)

runnotes <- "This is a bare-bones example."
# If you unhash the make_html call it will generate a local website 
# inside resdir and open it so you can see the results so far.
# make_html(replist=reportlist,resdir=resdir,width=500,
#          openfile=TRUE,runnotes=runnotes,verbose=FALSE,
#          packagename="aMSE",analysis="aMSE")
```

See the vignette Running\_aMSE.Rmd for a more detailed example.

See the New.Rd for recent developments
