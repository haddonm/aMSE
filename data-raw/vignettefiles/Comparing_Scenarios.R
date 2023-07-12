## ---- include = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE,
  message = FALSE,
  warning = FALSE
)

options(knitr.kable.NA = "",
        knitr.table.format = "pandoc")

options("show.signif.stars"=FALSE,
        "stringsAsFactors"=FALSE,
        "max.print"=50000,
        "width"=240)

suppressPackageStartupMessages({
library(aMSE)
library(TasHS)
library(codeutils)
library(hplot)
library(knitr)
library(makehtml)
library(captioner)
})

tab_nums <- captioner(prefix = "Table")
fig_nums <- captioner(prefix = "Figure")

dbdir <- getDBdir()
ddir <- paste0(dbdir,"/A_CodeUse/aMSEDoc/") 
prefixdir <- paste0(dbdir,"/A_CodeUse/aMSEUse/scenarios/")



## ----echo=TRUE, eval=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  # A subset of commands to run a given scenario
#  postfixdir <- "BC"                   # BC stands for basecase
#  rundir <- filenametopath(prefixdir,postfixdir)
#  outdir <- "C:/aMSE_scenarios/BC/"
#  confirmdir(rundir,ask=FALSE)     # check to see if it already exists
#  confirmdir(outdir,ask=FALSE)     # and if not create it without asking.
#  
#  hsargs <- list(mult=0.1, # expansion factor for cpue range when calc the targqnt
#                 wid = 4, # number of years in the grad4 PM
#                 targqnt = 0.55, # quantile defining the cpue target
#                 maxtarg = c(150,150,150,150,150,150,150,150), # max cpue Target
#                 pmwts = c(0.65,0.25,0.1),  # relative weights of PMs
#                 hcr = c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2), # hcr multipliers
#                 startCE = 2000, # used in constant reference period HS
#                 endCE = 2019,   # used in constant reference period HS
#                 metRunder = 1,  # should the metarules be used. o =
#                 metRover = 1,   # use metarules
#                 decrement=1, # use fishery data up to the end of the time series
#                 pmwtSwitch = 0, # number of years after reaching the targCE to
#                 stablewts = c(0.4, 0.5, 0.1), # replace pmwts with stablewts
#                 hcrname="mcdahcr",     # the name of the HCR used
#                 printmat=NULL)
#  
#  
#  out <- do_MSE(rundir,controlfile,hsargs=hsargs,hcrfun=mcdahcr,
#                sampleCE=tasCPUE,sampleFIS=tasFIS,sampleNaS=tasNaS,
#                getdata=tasdata,calcpopC=calcexpectpopC,makeouthcr=makeouthcr,
#                hcrscoreoutputs=extractandplotscores,HSPMs=getcpueHS,
#                fleetdyn=NULL,interimout=c(outdir,postfixdir),
#                varyrs=7,startyr=38,verbose=TRUE,ndiagprojs=4,cutcatchN=56,
#                matureL=c(40,170),wtatL=c(50,200),mincount=120,
#                includeNAS = FALSE,depensate=0,kobeRP=c(0.4,0.2,0.15))
#  
#  makeoutput(out,rundir,postfixdir,controlfile,hsfile="TasHS Package",
#             openfile=TRUE,verbose=FALSE)
#  
#  save(out,file=paste0(outdir,postfixdir,".RData"))  # the results are saved

