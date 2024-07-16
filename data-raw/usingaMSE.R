# R-code from the file  04_Using_aMSE.qmd

# ```{r}
#| echo: false
source("_common.R")
# ```



# ```{r }
#| label: fig-Rstudio-Install
#| echo: false
#| fig-cap: |
#|    The Install Packages interface from RStudio. The dropdown arrow is used to select 'Package Archive File (.zip; .tar.gz)', and then you browse to point the box at the required source file.
#| fig-alt: |
#|   Image of the Install Packages interface from RStudio.
#| out-width: 50%
filen <- pathtopath(basedir,"/figures/install_tar.gz_file.png")
knitr::include_graphics(filen,dpi=300)
# ```



# ```{r}
#| echo: true
#| label: setup
options("show.signif.stars"=FALSE, # some R options that I find can help
        "stringsAsFactors"=FALSE,  # now the default in R4
        "max.print"=50000,
        "width"=240)
suppressPackageStartupMessages({  # declare libraries -------------
  # this is the minimum, others can be added if desired
  library(aMSE)
  library(TasHS)     # obviously only if using the Tasmanian HS
  library(codeutils)
  library(hplot)
  library(makehtml)
  library(knitr)
})
# OBVIOUSLY, modify the rundir definition to suit your own setup!!!
dropdir <- getDBdir()
basedir <- pathtopath(dropdir,"A_Code/aMSEGuide/")
prefixdir <- pathtopath(basedir,"/runs/")
postfixdir <- "EG"     # a descriptive name for the rundir
rundir <- pathtopath(prefixdir,postfixdir) # define rundir
startime <- Sys.time() # to document the time taken
verbose <- TRUE        # send messages to the console about the run progress
controlfile <- paste0("control",postfixdir,".csv") # match control file name
outdir <- "C:/aMSE_scenarios/EG/"   # storage on a non-cloud hard-drive
confirmdir(rundir,ask=FALSE)   # make it if it does not exist
confirmdir(outdir,ask=FALSE)   # to be interactive needs ask = TRUE
# ```



# ```{r}
#| echo: true
#| label: make_data
controlfile <- "controlEG.csv" # name this CSV file whatever you wish
# Now make the controlfile. Read the help file using  ?ctrlfiletemplate. This
# explains what the devrec argument is all about.
ctrlfiletemplate(indir=rundir,filename=controlfile,devrec=0)
# Within controlfile, the default data file name is 'saudatapostfixdir.csv'. If
# you want to call it something else (maybe aloysius.csv or machynlleth.csv?)
# then edit the seventh line of the controlfile that holds the definition
# of the name of the SAU data file, and change it in this code.
datafiletemplate(indir=rundir,filename="saudataEG.csv")
# The program also needs some length composition data. For this we are going
# to use one of the inbuilt data-sets called 'lfs' (see its help). The default
# filen is already written into the controlfile.
data(lfs)
writecompdata(indir=rundir,lfs,filen="lf_WZ90-20.csv")
# dir(rundir)  # listing the contents of rundir can be useful
# ```



# ```{r}
#| echo: true
#| label: read_HS_functions_and_constants
# In Tasmania, the HS is now included in a package called TasHS, and this
# contains all the required functions.

# But, if the TasHS package is used instead of TasHS1_Tas.R, we still need
# the global object hsargs, which contains the settings used by the the
# Tasmanian HS. For details of hsargs, see the documentation for ?TasHS.

hsargs <- list(mult=0.1, #expansion factor for cpue range when calc the targqnt
               wid = 4, # number of years in the grad4 PM
               targqnt = 0.55, # quantile defining the cpue target
               maxtarg = c(150,150,150,150,150,150,150,150), # max cpue Target
               pmwts = c(0.65,0.25,0.1),  # relative weights of PMs
               hcr = c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2), # hcr mults
               hcrm3 = c(0.25,0.75,0.8,0.85,0.9,1,1.1,1.15,1.2,1.25),# metar 3
               startCE = 2000, # used in constant reference period HS
               endCE = 2011,   # used in constant reference period HS
               metRunder = 2,  # should the metarules be used. 0 = No
               metRover = 2,   # use metarules 0 = No
               decrement=1, # use fishery data up to the end of the time series
               pmwtSwitch = 4, # n years after reaching the targCE to replace
               stablewts = c(0.4, 0.5, 0.1), # replace pmwts with stablewts mr3
               hcrname="mcdahcr",     # the name of the HCR used
               printmat=NULL) # something needed for some HS
# ```



# ```{r}
#| echo: true
#| eval: false
tasFIS <- function(x) { # currently no FIS data is used in TAS
  return(NULL)          # though this may change
}
# ```



# ```{r}
#| echo: true
#| label: do conditioning
prodout <- FALSE  # no estimation of production properties saves time but
out <- do_condition(rundir,controlfile,  # gives no details re MSY
                    calcpopC=calcexpectpopC,
                    verbose = TRUE,
                    doproduct = prodout,
                    dohistoric=TRUE,
                    matureL = c(70, 200), # length range for maturity plots
                    wtatL = c(80, 200), # lengths for weight-at-length plots
                    mincount=120,  # minimum obs for including length comps
                    uplimH=0.35,   # not used because doproduct=FALSE
                    incH=0.005,
                    deleteyrs=0,  # all length comp years used
                    prodpops=NULL) # no individual pop producivity plotted

makeoutput(out,rundir,postfixdir,controlfile,hsfile="TasHS Package",
           doproject=prodout,openfile=TRUE,verbose=FALSE)
# ```



# ```{r}
#| echo: true
#| label: doing_MSE
checkhsargs(hsargs)  # remind user of selected options before
                     # just ot be sure!
# ?do_MSE provides more detailed descriptions of each function argument
# see TasHS package for help with the HS functions
# do_MSE function can take a few minutes to run so be patient
out <- do_MSE(rundir,controlfile, # already known
              hsargs=hsargs,    # defined as global object
              hcrfun=mcdahcr,   # the main HS function
              sampleCE=tasCPUE, # processes cpue data
              sampleFIS=tasFIS, # processes FIS data (see above)
              sampleNaS=tasNaS, # processes Numbers-at-Size data
              getdata=tasdata,  # extracts the data from zoneDP objects
              calcpopC=calcexpectpopC, #distributes catches to populations
              makeouthcr=makeouthcr, # generates updateable HS object
              fleetdyn=NULL,   # only used in SA so far
              scoreplot=plotfinalscores, # a function to plot HS scores
              plotmultflags=plotmultandflags, # fun to plot multipliers
              interimout="", # save the results after projecting,
              varyrs=7, # years prior to projections for random recdevs
              startyr=48, # in plots of projections what year to start
              verbose=TRUE, # send progress reports to the console
              ndiagprojs=4, # individual trajectories in DiagProj plots
              cutcatchN = 56,  # reduce size-at-catch array below this class
              matureL = c(70, 200), # length range for maturity plots
              wtatL = c(50, 200),   # length range for weight-at-length plots
              mincount = 120,    # minimum sample size in size-composition data
              includeNAS = FALSE, # include the numbers-at-size in saved output?
              depensate=0, # will depensation occur? see do_MSE help for details
              kobeRP=c(0.4,0.2,0.15), # ref points for a kobe-like phase plot
              nasInterval=5,  # year interval to use when plotting pred NaS
              minsizecomp=c(100,135), #min size for pred sizes, min for catches
              uplimH=0.4,incH=0.005, #H range when estimating productivity
              fissettings=NULL,    # no FIS in Tasmania
              fisindex=NULL,
              deleteyrs=0)    # all length comp years used
# ```



# ```{r}
#| echo: true
str1(out)  # function from codeutils = str(out,max.levels=1)
# ```



# ```{r}
#| label: make_MSE_web_page
#| echo: true
makeoutput(out,rundir,postfixdir,controlfile,
           hsfile="TasHS Package",openfile=TRUE,verbose=FALSE)
# ```



# ```{r}
#| label: dont save the output object
#| eval: false
save(out,file=pathtopath(outdir,paste0(postfixdir,".RData")))
# ```



# ```{r }
#| label: fig-make-full-website
#| echo: false
#| fig-cap: |
#|   The home page of the internal web-site generated from the results of the do_MSE R function from aMSE.  Of course, your own will differ from this in the details but all the tabs should be there.
#| fig-alt: |
#|   Image of he home-page of the internal web-site generated from the results of the do_MSE R function from aMSE.
#| out-width: 100%
filen <- pathtopath(basedir,"/figures/homepage_aMSERun.png")
knitr::include_graphics(filen,dpi=270)
# ```



# ```{r }
#| label: fig-back-arrow
#| echo: false
#| fig-cap: |
#|   A view of a browser page where a particular figure has been enlarged by clicking on it (here most of it is obscured to the right). Use the highlighted arrow to return to the ordinary scale view.
#| fig-alt: |
#|   When clicking on a figure it enlarges. This image highlights the backarrow to use to return the figure to its usual size.
#| out-width: 800%
filen <- pathtopath(basedir,"/figures/topleft_return_arrow.png")
knitr::include_graphics(filen,dpi=270)
# ```



# ```{r }
#| label: fig-compareCPUE
#| echo: false
#| fig-cap: |
#|   The top plot from the condition page of the internal web-page generated from the results of the do_condition and do_MSE R functions from aMSE. The values beside the SAU names are the simple sum-of-squared differences between the observed and predicted values.
#| fig-alt: |
#|   Image of how well the observed CPUE is matched by the conditioned model.
#| out-width: 40%
filen <- pathtopath(rundir,"/compareCPUE.png")
knitr::include_graphics(filen)
# ```



# ```{r }
#| label: fig-04-sau10-conditioned
#| echo: false
#| fig-cap: |
#|   The SAU plot illustrating the dynamics of the combined populations within SAU 10. This, and plots for all other SAU are to be found in the 'condition' tab.
#| fig-alt: |
#|   An image illustrating the quality of the conditioning for each sau in the simulated zone.
#| out-width: 80%
filen <- pathtopath(rundir,"/SAU10_conditioned.png")
knitr::include_graphics(filen,dpi=270)
# ```



