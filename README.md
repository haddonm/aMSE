
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LATEST UPDATE

  - 2021-04-20 0.0.0.2900 Continuing process of modifying the core
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
  ddir <- "c:/Users/User/DropBox/A_code/"
} else {
  ddir <- "c:/Users/Malcolm/DropBox/A_code/"
}
rundir <- paste0(ddir,"aMSEUse/scenarios/HS652510")
datadir <- paste0(ddir,"aMSEUse/scenarios/tasdata")
alldirExists(rundir,datadir)
#> rundir,  c:/Users/User/DropBox/A_code/aMSEUse/scenarios/HS652510 :  exists  
#> datadir,  c:/Users/User/DropBox/A_code/aMSEUse/scenarios/tasdata :  exists
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
source(paste0(rundir,"/TasmanianHS.R"))
# generate equilibrium zone ----------------------------------------------------
starttime2 <- (Sys.time())
zone <- makeequilzone(rundir,"controlsau.csv",datadir=datadir,cleanslate = TRUE)
#> All required files appear to be present 
#> Files read, now making zone 
#> Now estimating population productivity 
#> matureB Stable 
#> exploitB Stable 
#> recruitment Stable 
#> spawning depletion Stable
equiltime <- (Sys.time()); print(equiltime - starttime2)
#> Time difference of 20.96232 secs
# declare main objects ---------------------------------------------------------
glb <- zone$glb
ctrl <- zone$ctrl
zone1 <- zone$zone1
projC <- zone1$projC
condC <- zone1$condC
zoneC <- zone$zoneC
zoneD <- zone$zoneD
product <- zone$product
# add some tables and plots into rundir
biology_plots(rundir, glb, zoneC)
plotproductivity(rundir,product,glb)
numbersatsize(rundir, glb, zoneD)
# Condition on Fishery ---------------------------------------------------------
zoneDD <- dohistoricC(zoneD,zoneC,glob=glb,condC,sigR=1e-08,sigB=1e-08)

propD <- t(getzoneprops(zoneC,zoneDD,glb,year=47))
addtable(round(propD,4),"propertyDD.csv",rundir,category="zoneDD",caption=
           "Properties of zoneD after conditioning on historical catches.")
addtable(round(t(zoneDD$harvestR[40:47,]),4),"final_harvestR.csv",rundir,
         category="zoneDD",caption="Last eight years of harvest rate.")
popdefs <- getlistvar(zone$zoneC,"popdef")
addtable(round(t(popdefs),3),"popdefs.csv",rundir,category="zoneDD",caption=
           "Population specific definitions")
# Prepare projections ----------------------------------------------------------
plotconditioning(zoneDD,glb,zoneC,condC$histCE,1973:2019,rundir)

#explicit definition of indata rather than using the tasdata function
indata <- list(arrce=condC$histCE,yearnames=rownames(condC$histCE),
               acatches=tail(condC$histCatch,1),fis=NULL,nas=NULL)
cmcda <- mcdahcr(indata,hsargs=hsargs,saunames=glb$saunames)
pms <- cmcda$pms
multTAC <- cmcda$multTAC
out <- prepareprojection(projC=projC,condC=condC,zoneC=zoneC,glb=glb,
                            multTAC=multTAC,zoneDD=zoneDD,ctrl=ctrl,varyrs=7,
                            lastsigR = ctrl$withsigR)
# prepareprojection needs to be modified so that is uses expected catches rather
# than calculating them internally (which will differ by jurisdiction)
zoneDP <- out$zoneDP
projC <- out$projC
zoneCP <- out$zoneCP
zoneDDR <- out$zoneDDR

zoneDP <- doprojections(ctrl,zoneDP,zoneCP,condC$histCE,glb,mcdahcr,hsargs,
                        tasCPUE,tasFIS,tasNaS)
# save more plots to rundir ----------------------------------------------------
histCE <- condC$histCE
B0 <- getvar(zoneC,"B0")
ExB0 <- getvar(zoneC,"ExB0")
sauout <- sauplots(zoneDP,zoneDDR,glb,rundir,B0,ExB0,startyr=25,addCI=TRUE,
                  histCE=histCE)
diagnosticsproj(sauout$zonePsau,glb,rundir,nrep=3)
outzone <- poptozone(zoneDP,glb,
                     B0=sum(getvar(zoneC,"B0")),
                     ExB0=sum(getvar(zoneC,"ExB0")))
plotZone(outzone,rundir,CIprobs=c(0.05,0.5,0.95),addfile=TRUE)
projtime <- Sys.time()
print(projtime - starttime)
#> Time difference of 58.96815 secs
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
#>  $ matureB : num [1:30, 1:28, 1:100] 27.9 30.1 30.8 29.7 30.1 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ exploitB: num [1:30, 1:28, 1:100] 25.8 21.9 23 24 23.8 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ midyexpB: num [1:30, 1:28, 1:100] 30.5 27.3 29.3 30.5 29.3 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ catch   : num [1:30, 1:28, 1:100] 2.68 3.76 4.49 4.57 3.63 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ acatch  : num [1:30, 1:8, 1:100] 13.6 15 15 15 15.7 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ harvestR: num [1:30, 1:28, 1:100] 0.104 0.172 0.195 0.19 0.153 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ cpue    : num [1:30, 1:28, 1:100] 186 162 172 180 175 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ cesau   : num [1:30, 1:8, 1:100] 142 121 125 140 159 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ catsau  : num [1:30, 1:8, 1:100] 14.2 15.4 17.7 17.3 15 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ recruit : num [1:30, 1:28, 1:100] 62962 18512 24351 16443 19163 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ deplsB  : num [1:30, 1:28, 1:100] 0.402 0.433 0.443 0.427 0.432 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ depleB  : num [1:30, 1:28, 1:100] 0.402 0.342 0.359 0.375 0.371 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ catchN  : num [1:105, 1:30, 1:28, 1:100] 3.73e-132 1.88e-128 6.97e-125 1.87e-121 3.65e-118 ...
#>   ..- attr(*, "dimnames")=List of 4
#>  $ Nt      : num [1:105, 1:30, 1:28, 1:100] 6.30e+04 3.07e-07 1.26e-05 3.74e-04 8.08e-03 ...
#>   ..- attr(*, "dimnames")=List of 4
#>  $ NumNe   : num [1:105, 1:30, 1:28, 1:100] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..- attr(*, "dimnames")=List of 4
#>  $ TAC     : num [1:30, 1:100] 703 550 497 496 498 ...
#>   ..- attr(*, "dimnames")=List of 2
```

To see the structure of the dynamic object generated by the projections.

See the vignette Running\_aMSE.Rmd for a more detailed example.

See the New.Rd for recent developments
