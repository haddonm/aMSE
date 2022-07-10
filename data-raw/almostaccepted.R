




options("show.signif.stars"=FALSE,"stringsAsFactors"=FALSE,
        "max.print"=50000,"width"=240)
suppressPackageStartupMessages({
  library(aMSE)
  library(TasHS)
  library(hutils)
  library(hplot)
  library(makehtml)
  library(knitr)
})
dropdir <- getDBdir()
prefixdir <- paste0(dropdir,"A_codeUse/aMSEUse/scenarios/")

startime <- Sys.time()
postfixdir <- "S21"
verbose <- TRUE
rundir <- filenametopath(prefixdir,postfixdir)
controlfile <- paste0("control",postfixdir,".csv")
confirmdir(rundir)


data(zone1)
data(saudat)


rewritecontrolfile(rundir,zone1,controlfile=controlfile)
rewritedatafile(rundir,zone1,saudat)
rewritecompdata(rundir,zone1)

# need to rename the control and datafiles to remove the '_new' postfix.

hsargs <- list(mult=0.1,
               wid = 4,
               targqnt = 0.55,
               maxtarg = c(150,150,150,150,150,150,150,150),
               pmwts = c(0.65,0.25,0.1),
               hcr = c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2),
               startCE = 1992)



out <- do_condition(rundir,controlfile,
                    calcpopC=calcexpectpopC,
                    verbose = TRUE,
                    doproduct = TRUE,
                    dohistoric=TRUE,
                    mincount=120)

aMSE::makeoutput(out,rundir,postfixdir,controlfile,hsfile="TasHS Package",
                 doproject=FALSE,openfile=TRUE,verbose=FALSE)




adjustavrec(rundir,out$glb,out$ctrl,calcpopC=calcexpectpopC,verbose=TRUE)



























