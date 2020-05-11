# Used in readme.Rmd to illustrate the running of aMSE

# a minimal example
starttime <- as.character(Sys.time())
library(aMSE)
# Obviously you should modify the rundir to suit your own computer
rundir <- "./../../rcode2/aMSEUse/run2"
outdir <- setupdirs(rundir)
datadir <- outdir$datadir
resdir <- outdir$resdir

data(ctrl)
data(constants)
data(region1)
glb <- region1$globals
runname <- ctrl$runlabel

out <- setupregion(constants, glb, region1)
regionC <- out$regionC
regionD <- out$regionD
product <- out$product

regDe <- testequil(regionC,regionD,glb)

resfile <- setuphtml(resdir,runname)

plotproductivity(resdir,runname,product,glb)
biology_plots(resdir, runname, glb, regionC, regionD, product)
numbersatsize(resdir, runname, glb, regionD)

endtime <- as.character(Sys.time())


reportlist <- list(runname=runname,
  starttime=starttime,endtime=endtime,
  regionC=regionC,
  regionD=regionD,
  product=product,
  glb=glb,constants=constants
)

runnotes <- "This is a bare-bones example."
# If you unhash this component it will generate a local website inside
# resdir and open it so you can see the results so far.
# make_html(replist=reportlist,rundir=rundir,width=500,
#          openfile=TRUE,runnotes=runnotes,verbose=FALSE)

