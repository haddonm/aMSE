# Used in readme.Rmd to illustrate the running of aMSE

# a minimal example
starttime <- as.character(Sys.time())
library(aMSE)
# Obviously you should modify the resdir to suit your own computer
resdir <- "./../../rcode2/aMSEUse/out/testrun"
dirExists(resdir,make=TRUE,verbose=TRUE)
# You now need to ensure that there is a control.csv, reg1smu2pop6.csv
# and region1.csv file in the data directory
ctrlfiletemplate(resdir)
regionfiletemplate(resdir)
datafiletemplate(6,resdir,filename="reg1smu2pop6.csv")
ctrl <- checkresdir(resdir)
runname <- ctrl$runlabel
region1 <- readregionfile(resdir,ctrl$regionfile)
glb <- region1$globals
constants <- readdatafile(glb$numpop,resdir,ctrl$datafile)

out <- setupregion(constants, glb, region1) # make operating model
regionC <- out$regionC
regionD <- out$regionD
product <- out$product     # important bits usually saved in resdir
          # did the larval dispersal level disturb the equilibrium?
regDe <- testequil(regionC,regionD,glb)

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
# If you unhash this component it will generate a local website inside
# resdir and open it so you can see the results so far.
# make_html(replist=reportlist,resdir=resdir,width=500,
#          openfile=TRUE,runnotes=runnotes,verbose=FALSE)

