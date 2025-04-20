
<!-- README.md is generated from README.Rmd. Please edit that file -->

# aMSE

<!-- badges: start -->
<!-- badges: end -->

This is the codebase for an Abalone Management Strategy Evaluation
framework, which we call **aMSE**. By ‘we’ the intent is Malcolm Haddon
and Craig Mundy of University of Tasmania, IMAS-FA, Cathy Dichmont of
CDC-Consultancy, and Owen Burnell of SARDI, SA.

The MSE has a novel structure that is based around the spatial ideas of
populations (or areas of persistent production or apps, reflecting how
they are defined), within spatial assessment units (SAUs), within a
zone:

- zone - highest geographical level. Is simply the totality of the
  spatial assessment units, or more simply a quota zone.

- SAU - spatial assessment units. In Tasmania these would currently be
  the classical statistical blocks within each quota zone.

- population/app - a population, herein defined as an area of persistent
  production. These are the active components within the simulation
  framework. The dynamics of the simulation are based around the
  populations, although, with positive larval dispersal (a low level is
  the default) there is some, although very limited, dependency of
  neighbouring populations on each other.

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

**aMSE** and the other packages used, all come with a version number.
This is because they are all under different levels of on-going
development. For example, **makehtml** hasn’t changed since September
2024, and that was only a minor change to the *addtext* function. Other
packages changes more often. **aMSE** had some additions included on 20
April 2025. In fact, **aMSE** may be expected to continue being
developed as it is a useful tool in the on-going examination and
improvement of existing harvest strategies and the code is added to each
time something novel is required. However, such additions are designed
so as not to disturb the operation and functionality of what is already
available. This paragraph is simply here to warn that changes may occur
but to point out that these should not detract from everyday use.

## Installation

You can install the development version from
[GitHub](https://github.com/haddonm/aMSE) with:

``` r
if (!require(devtools)){install.packages("devtools")} 

devtools::install_github("https://github.com/haddonm/aMSE")

# similarly for codeutils, hplot, EGHS, and makehtml, which are also required
# knitr you should get from CRAN
```

The working of the **aMSE** software and how to use it has been
documented during the development and this is available as a GitBook at:

<https://haddonm.github.io/aMSEGuide/>

That can either be read online or a PDF of the whole downloaded (look
for the small Acrobat symbol at top left and click that to download the
latest published version). AMSEGuide can be considered as a ‘living’
document in that it may receive updates and additions as extra
functionality is added to the **aMSE** codebase.

The version number of the first publicly available version is 1.0.0,
with a date of 20-04-2025. This will alter if changes are made to the
software.

## Example

This is an example which illustrates the generation of an initial
equilibrium, which then goes on to apply an early version of the
Tasmania Multi-Criterion Decision Analysis (MCDA) using only 50
replicates (that bit takes about 15 seconds on my computer, hence 2.5
minutes for 1000). It uses built in data-sets but usually you would read
in a control file, which would contain the name of the biological data
file describing each population.

``` r
# a constant TAC example
starttime <- (Sys.time())
options("show.signif.stars"=FALSE,
        "stringsAsFactors"=FALSE,
        "max.print"=50000,
        "width"=240)
# declare libraries ------------------------------------------------------------
suppressPackageStartupMessages({
  library(aMSE)       # from github.com/haddonm/aMSE MSE codebase
  library(EGHS)       # from github.com/haddonm/EGHS example Tas harvest strategy
  library(codeutils)  # from github.com/haddonm/codeutils a set of code utilities
  library(hplot)      # from github.com/haddonm/hplot a set of plotting routines
  library(makehtml)   # from github.com/haddonm/makehtml used to display results
  library(knitr)      # from CRAN for generating tables
})
dbdir <- getDBdir()   # utility to define the dropbox directory
# Obviously you should modify the rundir and datadir to suit your own setup
prefixdir <- paste0(dbdir,"A_CodeUse/aMSEUse/scenarios/")
postfixdir <- "M15h75"      # name your scenario, here natural mortality=0.15, steepness=0.75 
rundir <- paste0(prefixdir,postfixdir) # files and results directory
confirmdir(rundir)
controlfile <- "controlM15h75.csv"    # name of the control file
# equilibrium zone -------------------------------------------------------------
# You now need to ensure that there is, at least, a control??.csv, and a 
# saudata??.csv file in the data directory plus some other data .csv files
# depending on how conditioned you want the model to be. Templates for the
# correct format can be produced using ctrlfiletemplate(), datafiletemplate().
# 
# Of course, usually one would use data files, control.csv and a saudata.csv, which
# is listed as the datafile within the control.csv. These must be stored in 
# rundir. Read the Using_aMSE' section in the Documentation .docx file in the 
# documentation. To see examples of code to run the MSE for a real question.

# Define the arguments used in the harvest strategy
hsargs <- list(mult=0.1,
               wid = 4,
               targqnt = 0.55,
               maxtarg = c(150,150,150,150,150,150,150,150),
               pmwts = c(0.65,0.25,0.1),
               hcr = c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2),
               hcrm3 = c(0.25,0.75,0.8,0.85,0.9,1,1.1,1.15,1.2,1.25),
               startCE = 1992,
               endCE = 2020,
               metRunder = 2,
               metRover = 2,
               decrement=1,
               pmwtSwitch = 4,
               stablewts = c(0.4, 0.5, 0.1),
               hcrname="constantrefhcr",
               printmat=NULL)

hsargs <- checkhsargs(hsargs)

deleteyrs <- 0  # none removed in the western zone

# One would need to unhash this to get it to run (obviously)

#out <- do_MSE(rundir,controlfile,hsargs=hsargs,hcrfun=constantrefhcr,
#              sampleCE=tasCPUE,sampleFIS=tasFIS,sampleNaS=tasNaS,
#              getdata=tasdata,calcpopC=calcexpectpopC,makeouthcr=makeouthcr,
#              fleetdyn=NULL,scoreplot=plotfinalscores,
#              plotmultflags=plotmultandflags,interimout="",
#              varyrs=7,startyr=38,verbose=TRUE,ndiagprojs=4,cutcatchN=56,
#              matureL=c(40,170),wtatL=c(50,200),mincount=120,
#              includeNAS = FALSE,depensate=0,kobeRP=c(0.4,0.2,0.15),
#              nasInterval=1,minsizecomp=c(100,135),uplimH=0.375,incH=0.005,
#              deleteyrs=deleteyrs,selectyr=0)

#makeoutput(out,rundir,postfixdir,controlfile,hsfile="TasHS Package",
#           openfile=TRUE,verbose=FALSE)

# save(out,file=paste0(outdir,postfixdir,".RData"))
```

After running the whole, even if you do not generate the results
webpage, you could also try:

``` r
#str(out$zoneDP,max.level=1)
#str1(out$hcrout)
```

To see the structure of the dynamic object generated by the projections.

See the *aMSEGuide* for detailed examples. This living document, which
is also open to on-going development, provides many more details of the
operation and structure of the software components used, but also many
more examples of scenarios runs and the comparison of such scenarios,
which after all, is the point of the software.

The **EGHS** package is a snapshot of the Tasmanian constant reference
period harvest strategy as of December 2024. It provides a working HS
for examples and for any other user to get started. Similar example HS
are not currently available for South Australia or Victoria. The example
HS provides the capacity to explore the effects of changing details
within the various HS arguments used in Tasmania (as expressed in the
*hsargs* object). As with all these details, see the *aMSEGuide* for
details.

### Acknowledgments

All the coding and documentation for this Management Strategy Evaluation
(MSE) framework was part of the FRDC project 2019-118, with funding
provided by the Fisheries Research Development Corporation, as well as
by the University of Tasmania and the South Australian Research and
Development Institute (SARDI). The occurrence of the Covid 19 pandemic
greatly complicated the development and progress of this project. The
eventual complexity of coding for the MSE, and the coding of each
jurisdiction’s harvest strategy took far longer than initial estimates.
This led to significant additional contributions from Malcolm Haddon and
Cathy Dichmont.

Malcolm Haddon takes responsibility for the software in the ‘haddonm’
GitHub based packages used here (and hence any remaining bugs; although
for TasHS any remaining bugs are shared with Craig Mundy!). Colleagues
in the FRDC MSE project, Cathy Dichmont, Craig Mundy, and Owen Burnell
were vital in checking that the MSE code operated as intended, asking
the right questions and making the right requests and suggestions. By
providing the necessary scientifically critical context for the work we
managed to succeed with what proved to be a very difficult task. Writing
the R code to capture the workings of the harvest strategies in each
jurisdiction was an especially complicated task. For Tasmania this was
completed by Malcolm Haddon and Craig Mundy, while that for South
Australia was a combination of Johnathon Smart, Cathy Dichmont, and Owen
Burnell. Cathy Dichmont wrote the R code for the Victorian Harvest
Strategy by modifying code written by Jo Potts from MRAG.
