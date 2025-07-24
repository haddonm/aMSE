
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

## Constant Reference Years Harvest Strategy

This is an example which illustrates the generation of an initial
equilibrium, which then goes on to apply an early version of the
Tasmania Multi-Criterion Decision Analysis (MCDA) using only 50
replicates (that bit takes about 15 seconds on my computer, hence 2.5
minutes for 1000). It uses built in data-sets but usually you would read
in a control file, which would contain the name of the biological data
file describing each population.

``` r
# The following is copied from aMSEGuide Chapter 3
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
# The next four code lines define my directory structure. 
# OBVIOUSLY, modify the rundir definition to suit your own setup!!!
dbdir <- getDBdir()   # utility to identify the dropbox directory
prefixdir <- pathtopath(dbdir,"A_codeR/aMSEGuide/runs/")
postfixdir <- "EG"     # a descriptive name for the rundir 
rundir <- pathtopath(prefixdir,postfixdir) # define rundir
startime <- Sys.time() # to document the time taken
verbose <- TRUE        # send messages to the console about the run progress
controlfile <- paste0("control",postfixdir,".csv") # match control file name
outdir <- "C:/aMSE_scenarios/EG/"   # storage on a non-cloud hard-drive
confirmdir(rundir,ask=FALSE)   # make rundir if it does not exist
confirmdir(outdir,ask=FALSE)   # to be interactive needs ask = TRUE
# You now need to ensure that there is, at least, a control??.csv, and a 
# saudata??.csv file in the data directory plus some other data .csv files
# Templates for the correct format can be produced using ctrlfiletemplate()
# and datafiletemplate().
# 
# Of course, usually one would use data files, control.csv and a saudata.csv, 
# which is listed as the datafile within the control.csv. These must be stored 
# in rundir. Read the Using_aMSE' section in aMSEGuide.
controlfile <- "controlEG.csv" # default example file name, change as needed
# Now make the controlfile. Read the help file using ?ctrlfiletemplate. This
# explains what the devrec argument is all about.

ctrlfiletemplate(indir=rundir,filename=controlfile,devrec=0)
# Within controlfile, the default data file name is 'saudatapostfixdir.csv'. If
# you want to call it something else (maybe aloysius.csv or machynlleth.csv?)
# then edit the seventh line of the controlfile that holds the definition
# of the name of the SAU data file, and change it in this code.

datafiletemplate(indir=rundir,filename="saudataEG.csv")
# The program also needs some length composition data. For this we are going 
# to use one of the inbuilt data-sets called 'lfs' (see its help). The default
# filen is already written into the controlfile. Change that as necessary.

data(lfs)
writecompdata(indir=rundir,lfs,filen="lf_WZ90-20.csv")

# Define the arguments used in the constantrefhcr harvest strategy
hsargs <- list(mult=0.1, #expansion factor for cpue range when calc the targqnt
               wid = 4, # number of years in the grad4 PM
               targqnt = 0.55, # quantile defining the cpue target
               maxtarg = c(150,150,150,150,150,150,150,150), # max cpue Target
               pmwts = c(0.65,0.25,0.1),  # relative weights of PMs
               hcr = c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2), # hcr mults
               hcrm3 = c(0.25,0.75,0.8,0.85,0.9,1,1.1,1.2,1.25,1.3),# mr 3
               startCE = 2000, # used in constant reference period HS
               endCE = 2019,   # used in constant reference period HS
               metRunder = 0,  # should the metarules be used. 0 = No
               metRover = 0,   # use metarules 0 = No
               decrement=1, # use fishery data up to the end of the time series
               pmwtSwitch = 0, # n years after reaching the targCE to replace
               stablewts = c(0.8, 0.15, 0.05), # pmwts with stablewts mr3
               hcrname="constantrefhcr",     # the name of the HCR used
               printmat=NULL) # something needed for some HS, not TAS

checkhsargs(hsargs) # lists chosen options before run, only for TAS 
# if you want to check your own set then include a suitable check function 
# within your jurisdiction's HS source file or package (see EGHS as example).

deleteyrs <- 0  # no years of data removed from the zone

# ?do_MSE provides more detailed descriptions of each function argument
# see EGHS package for help with the HS functions
# do_MSE function can take a few minutes to run so be patient
out <- do_MSE(rundir,controlfile, # already known
              hsargs=hsargs,    # defined as global object
              hcrfun=constantrefhcr,   # the main HS function
              sampleCE=tasCPUE, # processes cpue data
              sampleFIS=tasFIS, # processes FIS data (if any)
              sampleNaS=tasNaS, # processes Numbers-at-Size data
              getdata=tasdata,  # extracts the data from the zoneDP object
              calcpopC=calcexpectpopC, #distributes catches to populations
              makeouthcr=makeouthcr, # generates updateable HS output object
              fleetdyn=NULL,   # only used in SA and VIC so far
              scoreplot=plotfinalscores, # a function to plot HS scores
              plotmultflags=plotmultandflags, # function to plot multipliers
              interimout="", # save the results after projecting,
              varyrs=7, # years prior to projections for random recdevs
              startyr=48, # in plots of projections what year to start
              verbose=TRUE, # send progress reports to the console
              ndiagprojs=4, # individual trajectories in DiagProj plots
              cutcatchN = 56,  # cut size-at-catch array < class = 112mm
              matureL = c(70, 200), # length range for maturity plots
              wtatL = c(50, 200),   # length range for weight-at-length plots
              mincount = 120,   # minimum sample size in size-composition data
              includeNAS = FALSE, #include the numbers-at-size in saved output?
              depensate=0, #will depensation occur? see do_MSE help for details
              kobeRP=c(0.4,0.2,0.15), # ref points for a kobe-like phase plot
              nasInterval=5,  # year interval to use when plotting pred NaS 
              minsizecomp=c(100,135), #min size for pred sizes, min for catches
              uplimH=0.35,incH=0.005, #H range when estimating productivity 
              deleteyrs=0)    # all length comp years used 

makeoutput(out,rundir,postfixdir,controlfile,
           hsfile="EGHS Package",openfile=TRUE,verbose=FALSE)

# unhash this to save the output to your outdir
# save(out,file=paste0(outdir,postfixdir,".RData"))
```

After running the whole, even if you do not generate the results
webpage, to see the structure of the dynamic object generated by the
projections you could try:

``` r
str(out$zoneDP,max.level=1)
str1(out$hcrout)
```

### Constant Catch Harvest Strategy

Within EGHS there are two harvest control rules. The code above uses the
constant reference years harvest strategy ‘constantrefhcr’. There is
also the ‘consthr’ which simply implements a constant unvarying
aspirational catch in each SAU each year. To run that you should first
copy the required data files into a new directory. As before, setup the
R code to match whatever directory structure you intend to use:

``` r
suppressPackageStartupMessages({
  library(aMSE)       # from github.com/haddonm/aMSE MSE codebase
  library(EGHS)       # from github.com/haddonm/EGHS example Tas harvest strategy
  library(codeutils)  # from github.com/haddonm/codeutils a set of code utilities
  library(hplot)      # from github.com/haddonm/hplot a set of plotting routines
  library(makehtml)   # from github.com/haddonm/makehtml used to display results
  library(knitr)      # from CRAN for generating tables
})
# ONCE AGAIN YOU WILL NEED TO SETUP YOUR OWN DIRECTORY STRUCTURE  TO DEFINE
# THE prefixdir
dropdir <-getDBdir()   
prefixdir <- pathtopath(dropdir,"/A_codeR/aMSETas/scenarios/")
filelist=c("controlEG.csv","lf_WZ90-20.csv","saudataEG.csv")
aMSE::copyto(prefixdir=prefixdir,fromdir="EG",todir="EGconst",
              filelist=filelist)
# this copies the files across and renames them appropriately.
# NO CHANGES ARE NEEDED IN THE FILES BECAUSE THE CHANGES OCCUR IN hsargs
```

The change in the harvest strategy is implemented in the hsargs rather
than in the code, although it is necessary to let the do_MSE() function
know of the change:

``` r
# amend the run code above with the following
postfixdir <- "EGconst"
verbose <- TRUE
rundir <- pathtopath(prefixdir,postfixdir)
controlfile <- paste0("control",postfixdir,".csv")
outdir <- "C:/aMSE_scenarios/EG/"
confirmdir(rundir,ask=FALSE)
confirmdir(outdir,ask=FALSE)

# current Tas HS
hsargs <- list(mult=0.025, # cpue range expansion factor when calc the targqnt
               wid = 4, # number of years in the grad4 PM
               targqnt = 0.55, # quantile defining the cpue target
               maxtarg = c(150,150,150,150,150,150,150,150), # max cpue Target
               acatch =  c(22,48,26,125,105,210,175,47),#ONLY used for consthcr
               pmwts = c(0.65,0.25,0.10),  # relative weights of PMs default
               hcr = c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2), 
               hcrm3 = c(0.25,0.75,0.8,0.85,0.9,1,1.1,1.2,1.25,1.3),
               startCE = 2000, # used in constant reference period HS
               endCE = 2019,   # used in constant reference period HS
               metRunder = 2,  # MR2 if CE below targ multTAC = 1
               metRover = 2,   # MR1 CE>targ for 2+ yrs before multTAC increases
               decrement=1, # use fishery data up to the end of the time series
               pmwtSwitch = 4, # number of years after reaching the targCE to
               stablewts = c(0.65, 0.25, 0.1), # replace pmwts with stablewts
               hcrname="consthcr", # name of HCR: constantrefhcr / consthcr
               printmat=NULL)

hsargs <- checkhsargs(hsargs)


# NOTE THE CHANGE FROM tasdata TO constdata AND makeouthcr TO makeoutconst 
# IN THIE THIRD LINE
out <- do_MSE(rundir,controlfile,hsargs=hsargs,hcrfun=consthcr,  # constantrefhcr,     #
              sampleCE=tasCPUE,sampleFIS=tasFIS,sampleNaS=tasNaS,
              getdata=constdata,calcpopC=calcexpectpopC,makeouthcr=makeoutconst,
              fleetdyn=NULL,scoreplot=plotfinalscores,
              plotmultflags=plotmultandflags,interimout="",
              varyrs=7,startyr=38,verbose=TRUE,ndiagprojs=4,cutcatchN=56,
              matureL=c(40,170),wtatL=c(50,200),mincount=120,
              includeNAS = FALSE,depensate=0,kobeRP=c(0.4,0.2,0.15),
              nasInterval=1,minsizecomp=c(100,135),uplimH=0.35,incH=0.005,
              deleteyrs=0)

makeoutput(out,rundir,postfixdir,controlfile,hsfile="EGHS Package",
           openfile=TRUE,verbose=FALSE)
```

### Final Words

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
