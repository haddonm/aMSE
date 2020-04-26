# file reading ------------------------------------------------------
starttime <- as.character(Sys.time())
#library(rutilsMH)
library(aMSE)
library(microbenchmark)
rutilsMH::setpalette("R4")

# read data files ----------------------------------------------------

 sourcedir <- "./../../rcode2/aMSE/data-raw"
 source(filenametopath(sourcedir,"sourcer.R"))

 rundir <- "./../../rcode2/aMSEUse/run2"
 outdir <- setupdirs(rundir)
 datadir <- outdir$datadir
 plotdir <- outdir$plotdir

 ctrl <- readctrlfile(datadir,infile="control.csv")
 runname <- ctrl$runlabel
 region1 <- readregionfile(datadir,ctrl$regionfile)
 glb <- region1$globals
 constants <- readdatafile(glb,datadir,ctrl$datafile)

 # Define the Zone without production ---------------------------------
ans <- makeregionC(region1,constants)
regionC <- ans$regionC
popdefs <- ans$popdefs
ans <- makeregion(glb,regionC)
regionC <- ans$regionC  # region constants
regionD <- ans$regionD  # region dynamics
#  str(regionC[[1]])


# estimate production and move regionC to equilibrium-----------------
ans <- findunfished(regionC,regionD,glb)
regionC <- ans$regionC  # region constants
regionD <- ans$regionD  # region dynamics
product <- ans$product


# testing the equilibrium --------------------------------------------
# This runs the dynamics for Nyrs with zero harvest rate, a constant
# larval dispersal rate (from globals), and negligible recruitment
# variation and expects the calculated components to remain constant
# within three decimal places.
 regDe <- testequil(regionC,regionD,glb)
 #  str(regDe)

plottabfile <- setuphtml(plotdir,runname)

# characterize productivity ------------------------------------------
source(filenametopath(sourcedir,"plotproductivity_source.R"))
# end characterize productivity --------------------------------------

#characterize biology ------------------------------------------------
source(filenametopath(sourcedir,"plotbiology_source.R"))
# end characterize biology -------------------------------------------

if (ctrl$initdepl < 1.0) {
  regionD <- dodepletion(regionC, regionD, glb, product)
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

 source(filenametopath(sourcedir,"sourcer.R"))
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
 filename <- filenametopath(plotdir,file)
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




