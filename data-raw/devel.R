# file reading ------------------------------------------------------
starttime <- as.character(Sys.time())
library(rutilsMH)
library(aMSE)
library(microbenchmark)


# read data files ----------------------------------------------------

 sourcedir <- "./../../rcode2/aMSE/data-raw"
 #source(filenametopath(sourcedir,"sourcer.R"))

 rundir <- "./../../rcode2/aMSEUse/run2"
 outdir <- setupdirs(rundir)
 datadir <- outdir$datadir
 resdir <- outdir$resdir

 # You now need to ensure that there is a control.csv, reg1smu2pop6.csv
 # and region1.csv file in the data directory

 ctrl <- readctrlfile(datadir,infile="control.csv")
 runname <- ctrl$runlabel
 region1 <- readregionfile(datadir,ctrl$regionfile)
 glb <- region1$globals
 constants <- readdatafile(glb$numpop,datadir,ctrl$datafile)

 # source(filenametopath(sourcedir,"sourcer.R"))
 out <- setupregion(constants, glb, region1)
 regionC <- out$regionC
 regionD <- out$regionD
 product <- out$product

# save the primary ojects to the resdir
resfile <- setuphtml(resdir,runname)

save(ctrl,file=filenametopath(resdir,"ctrl.RData"))
save(glb,file=filenametopath(resdir,"glb.RData"))
save(product,file=filenametopath(resdir,"product.RData"))
save(regionC,file=filenametopath(resdir,"regionC.RData"))
save(regionD,file=filenametopath(resdir,"regionD.RData"))

# testing the equilibrium --------------------------------------------
# This runs the dynamics for Nyrs with zero harvest rate, a constant
# larval dispersal rate (from globals), and negligible recruitment
# variation and expects the calculated components to remain constant
# within three decimal places.
 regDe <- testequil(regionC,regionD,glb)
 #  str(regDe)

# characterize productivity and unfished biology ---------------------
plotproductivity(resdir,runname,product,glb)
biology_plots(resdir, runname, glb, regionC, regionD, product)
numbersatsize(resdir, runname, glb, regionC, regionD, product)
# end characterize biology -------------------------------------------
# store the initial properties
unfishprops <- getregionprops(regC=regionC,regD=regionD,glb=glb,year=1)
filename <- filenametopath(resdir,"unfishprops.csv")
write.table(round(unfishprops,4),file = filename,sep=",")
#  or use tmp <- read.csv(file=filename,header=TRUE,row.names=1)
caption <- paste0("The unfished equilibrium properties of the ",
              "populations and region before any initial depletion.")
logfilename(filename,resfile=resfile,"Tables",caption=caption)

ctrl$initdepl <-  0.40

if (ctrl$initdepl < 1.0) {
  regionDD <- dodepletion(regionC, regionD, glb, depl=ctrl$initdepl, product)
  # store the initial properties after depletion
  initprops <- getregionprops(regC=regionC,regD=regionDD,glb=glb,year=1)
  filename <- filenametopath(resdir,"initprops.csv")
  write.table(round(initprops,4),file = filename,sep=",")
   #  or use tmp <- read.csv(file=filename,header=TRUE,row.names=1)
  caption <- paste0("The equilibrium properties of the populations ",
                "and region after the initial depletion.")
  logfilename(filename,resfile=resfile,"Tables",caption=caption)
  save(regionDD,file=filenametopath(resdir,"regionDD.RData"))
}

 endtime <- as.character(Sys.time())


 reportlist <- list(
         runname=runname,
         starttime=starttime,
         endtime=endtime,
         regionC=regionC,
         regionD=regionD,
         product=product,
         glb=glb,
         constants=constants
 )
 str(reportlist,max.level = 1)

 runnotes <- paste0("The results presented here relate to the included data-sets testregC, ",
                    "testregD, and product. They are for a region made up of 2 SMU and 6 population. ",
                    "These results are currently under development and there are many more needed yet.")

#  source(filenametopath(sourcedir,"sourcer.R"))
 make_html(replist=reportlist,rundir=rundir,width=500,
           openfile=TRUE,runnotes=runnotes,verbose=FALSE)



# end of run ---------------------------------------------------------





 str(regionC[[1]])
 str(regionC[[4]])
 str(product)

 # Summarize MSY and related statistics ------------------------------
 approxMSY <- findmsy(product)
 approxMSY

 # Some summaries ----------------------------------------------------
 sapply(regionC,"[[","MSY")           # msy by population
 sum(sapply(regionC,"[[","MSY"))      # total msy
 sum(sapply(regionC,"[[","effB0"))    # total effective B0
 sum(sapply(regionC,"[[","B0"))       # total B0
 sum(sapply(regionC,"[[","effExB0"))  # total effective exploitable B0
 sum(sapply(regionC,"[[","ExB0"))     # total exploitable B0



 filen <- paste0("product_",runname,".RData")
 save(product,file=filenametopath(tabledir,filen))

 filen <- paste0("msytable_",runname,".csv")
 xval <- findmsy(product)
 write.csv(xval,file = filenametopath(tabledir,filen))

 filen <- filenametopath(tabledir,paste0("F1table_",runname,".csv"))
 xval <- findF1(product)
 write.csv(xval,file = filen)


 tmp <- read.csv(filen)


 # turn plot into a function

 file <- paste0("production_total_SpB_",ctrl$runlabel,".png")
 filename <- filenametopath(resdir,file)
 caption <- "The production curve relative to the depletion level of each population. The vertical lines identify the Depletion level giving rise to the MSY."
 category <- "Production"
 x <- product[,""]

 makexypng <- function(x,y,filen="",legendloc="topright",wid=7,hgt=4,
                         xlab="",ylab="") {
    plotprep(width=wid,height=hgt,filename=filen,verbose=FALSE)

 }



 xval <- findmsy(product)
 numpop <- glb$numpop
 plotprod(product,xname="Deplet",xlab="Population Depletion Level")
 for (pop in 1:numpop) abline(v=xval[pop,"Deplet"],lwd=2,col=pop)
 legend(legendloc,paste0("P",1:numpop),lwd=3,col=c(1:numpop),bty="n",
        cex=1.2)
 if (nchar(filen) > 0) dev.off()

  time <- as.character(Sys.time())
 cat(c(filename,caption,"production",time," \n"),file=plottabfile,sep=",",append=TRUE)

