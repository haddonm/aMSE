
<!-- README.md is generated from README.Rmd. Please edit that file -->

# aMSE

<!-- badges: start -->

<!-- badges: end -->

This packages up a new Abalone Management Strategy Evaluation framework.
It has a novel structure that is based around the spatial ideas of
populations within spatial management units (SMUs), within a region:

  - region - highest geogrpahical level. Is simply the totality of the
    spatial management units.

  - SMU - spatial management units. In Tasmania these would currently be
    the classical statistical blocks.

  - population - literally a population. These are the active components
    within the simulation framework. The dynamics of the simulation are
    based around the populations, although, with positive larval
    dispersal (the default) there is some dependency of neighbouring
    populations on each other.

This conceptual structure has implications for how the simulations are
conditioned during the management strategy testing. While it is still
the case that the general properties of each SMU are used to set the
scene for each SMU’s component populations, the advent of sufficient GPS
logger data now allows a more formal definition of the productivity of
each population. Currently, the modelling assumes that each population
is linearly arranged around the coastline, and this will allow, where
available, specific catch/yield data to be used to define the bounds of
each population within an SMU.

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

It would be a good idea to read Hadley Wickham’s draft chpater on Git
and GitHub at <https://r-pkgs.org/index.html>.

## Example

This is a basic example which illustrates the generation of an initial
equilibrium. It uses built in data-sets but usually you would generate
regionC (constants) and regionD (the dynamic components) from the
conditioning data:

``` r
# a minimal example
starttime <- as.character(Sys.time())
library(aMSE)
# Obviously you should modify the resdir to suit your own computer
resdir <- "./../../rcode2/aMSEUse/out/testrun"
dirExists(resdir,make=TRUE,verbose=TRUE)
#> ./../../rcode2/aMSEUse/out/testrun :  exists
# You now need to ensure that there is a control.csv, reg1smu2pop6.csv
# and region1.csv file in the data directory
ctrlfiletemplate(resdir)  # puts a template ctrlfile into resdir
regionfiletemplate(resdir)  # puts a template region file into resdir
datafiletemplate(6,resdir,filename="reg1smu2pop6.csv") # etc
ctrl <- checkresdir(resdir) # checks the data files are present
runname <- ctrl$runlabel    # and reads the control file, if present
region1 <- readregionfile(resdir,ctrl$regionfile)
glb <- region1$globals
constants <- readdatafile(glb$numpop,resdir,ctrl$datafile)

out <- setupregion(constants, glb, region1) # make operating model 
regionC <- out$regionC     # constant bits
regionD <- out$regionD     # dynamics bits
product <- out$product     # all usually saved in resdir, not here
  # did the larval dispersal level disturb the equilibrium?
regDe <- testequil(regionC,regionD,glb) 
#> [1] matureB Stable
#> [1] exploitB Stable
#> [1] recruitment Stable
#> [1] spawning depletion Stable

resfile <- setuphtml(resdir,runname)# prepare to save and log results

plotproductivity(resdir,runname,product,glb)
biology_plots(resdir, runname, glb, regionC)
numbersatsize(resdir, runname, glb, regionD)

endtime <- as.character(Sys.time())

reportlist <- list(runname=runname,
                   starttime=starttime,endtime=endtime,
                   regionC=regionC, regionD=regionD, product=product,
                   glb=glb,constants=constants
)

runnotes <- "This is a bare-bones example."
# If you unhash the make_html component it will generate a local 
# website inside resdir and open it so you can see the results so far.
#
# make_html(replist=reportlist,resdir=resdir,width=500,
#          openfile=TRUE,runnotes=runnotes,verbose=FALSE)
```

See the vignette Running\_aMSE.Rmd for a more detailed example.

See the New.Rd for recent developments
