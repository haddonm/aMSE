# Used in readme.Rmd to illustrate the running of aMSE

# a minimal example
starttime <- as.character(Sys.time())
library(aMSE)
library(rutilsMH)
library(makehtml)
# Obviously you should modify the resdir to suit your own computer
resdir <- "./../../A_code/aMSEUse/conddata/generic2"
dirExists(resdir,make=TRUE,verbose=TRUE)
# You now need to ensure that there is a control.csv, zone1sm\au2pop6.csv
# and region1.csv file in the data directory
zone1 <- readctrlzone(resdir,infile="control.csv")
ctrl <- zone1$ctrl
glb <- zone1$globals     # glb without the movement matrix
constants <- readdatafile(glb$numpop,resdir,ctrl$datafile)

#zone1$initLML <- 140
out <- setupzone(constants,zone1) # make operating model
zoneC <- out$zoneC
zoneD <- out$zoneD
glb <- out$glb        # glb now has the movement matrix
product <- out$product     # important bits usually saved in resdir
# did the larval dispersal level disturb the equilibrium?



regDe <- testequil(zoneC,zoneD,glb)

resfile <- setuphtml(resdir)# prepare to save and log results

plotproductivity(resdir,product,glb)
biology_plots(resdir, glb, zoneC)
numbersatsize(resdir, glb, zoneD)

endtime <- as.character(Sys.time())

reportlist <- list(starttime=starttime,endtime=endtime,
                   zoneC=zoneC, zoneD=zoneD, product=product,
                   glb=glb,constants=constants)

runnotes <- "This is a bare-bones example."
# If you unhash this component it will generate a local website inside
# resdir and open it so you can see the results so far.
 make_html(replist=reportlist,resdir=resdir,width=500,
          openfile=TRUE,runnotes=runnotes,verbose=FALSE,
          packagename = "aMSE",htmlname="testrun")

