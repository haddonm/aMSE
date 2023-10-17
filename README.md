
<!-- README.md is generated from README.Rmd. Please edit that file -->

# aMSE

<!-- badges: start -->
<!-- badges: end -->

This packages up a new Abalone Management Strategy Evaluation framework.
It has a novel structure that is based around the spatial ideas of
populations within spatial management units (SAUs), within a zone:

- zone - highest geographical level. Is simply the totality of the
  spatial management units.

- SAU - spatial assessment units. In Tasmania these would currently be
  the classical statistical blocks.

- population - literally a population. These are the active components
  within the simulation framework. The dynamics of the simulation are
  based around the populations, although, with positive larval dispersal
  (the default) there is some dependency of neighbouring populations on
  each other.

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

# LATEST UPDATES

- 2023-10-17 aMSE 0.1.13 added poplevelplot() plots to the output. So
  far this includes catches and cpue. Population level plots can give
  insight into the effectiveness of the fleet dynamics and whether each
  population is being exploited.

- 2023-08-01 aMSE 0.1.10 Added scoreplots and plotmultflags into the
  arguments for do_MSE rather than requiring the functions be called
  plotfinalscores() and plotmultandflags(), so any name can be given to
  them.

- 2023-07-26 aMSE 0.1.9 Added tasphaseplots to the kobe plots for
  comparison.

- 2023-07-24 aMSE 0.1.8 Revised use of TasHS and of comparescenarios,
  meta-rule flagging now working correctly.

- 2023-07-20 aMSE 0.1.7 Implemented the use of incremental updates of HS
  scores in a jurisdiction generic manner. Added a meta-rule flagging
  system.

- 2023-07-10 aMSE 0.1.6 Modified TasHS and aMSE so that instead of using
  hsargs$startCE and hsargs$endCE to define a vector of reference years
  when using the constantrefhcr in Tasmania, one should now use a
  hsargs\$refperiodCE which should be a vector of years describing the
  constant reference period to be used (assuming one is using the
  constantrefhcr). Within TasHS there is also a checkhsargs() function
  which provides a brief check of the hsargs being used (strictly that
  is still under development).

- 2023-06-11 aMSE 0.1.5 Modified the condition plot and fixed up a few
  places where numbers were hard wired into the software. Introduced a
  new option to save output following projections (just in case of
  experiments with subsequent analyses crash - I’m looking at you SA).
  Now one can set interimout =c(outdir,postfixdir) and that should save
  a tempout into your outdir directory. See the help for do_MSE.

- 2023-06-04 aMSE 0.1.4 Small modifications to generalize across
  jurisdictions and allow for the consthcr constant catch harvest
  strategy; also removed irrelavent messages when plotting the phase
  plot

- 2023-06-03 aMSE 0.1.3 Modified to allow the mcdahcr and constrefhcr
  HCR to run successfully with the conditioned SA setup.

- 2023-04-17 aMSE 0.1.2 Amended save_hsargs so it can handle data.frame
  within hsargs.

- 2023-04-12 aMSE 0.1.1 Started generalizing the outputs from aMSE by
  beginning to separate functions that relate to specific jurisdiction
  outputs into a separate R package (in Tasmania it is currently called
  aMSEExtra).

- 2023-04-04 aMSE 0.1.0 Added curryear to hcrfun so all HS functions
  should now include curryear inthier arguments. Within doprojecitons()
  it is set to year inside the loop. This was needed to allow the SA HS
  to use its fleet dynamics model and is not needed in Tasmania, but is
  required in the argument list.

- 2023-03-31 aMSE 0.0.25 Modified calcpopC througout to allow for SA’s
  harvest strategy and fleet dynamics.

- 2022-12-14 aMSE 0.0.24 Lots of minor edits and modifications aimed at
  improving the comparison of scenarios.

- 2022-11-11 aMSE 0.0.23 Added a new tab = recruits, containing plots of
  the stock recruitment relationships at an sau and at a population
  within sau level. Added a table to the Tables tab.

- 2022-11-10 aMSE 0.0.22 now have the capacity to introduce productivity
  events in specific projection years. This can now be done at sau
  specific levels. One defines the survivorship of recruits in each
  selected year and also the survivorship of all post-settlement
  animals. No changes are needed to one’s control file if this option is
  not being used.

- 2022-11-02 aMSE 0.0.21 added changes to how do_MSE responds to hsargs.

- 2022-10-31 aMSE 0.0.20 added do_comp_outputs, scenebyvar, scenebyzone,
  projectiononly, and sauquantbyscene to assist with comparing
  scenarios.

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
dbdir <- getDBdir()
# Obviously you should modify the rundir and datadir to suit your own setup
prefixdir <- paste0(dbdir,"A_CodeUse/aMSEUse/scenarios/")
postfixdir <- "M15h75"
rundir <- paste0(prefixdir,postfixdir)
confirmdir(rundir)
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

See the vignette Running_aMSE.Rmd for a more detailed example.

See the New.Rd for recent developments
