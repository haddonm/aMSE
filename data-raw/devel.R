# file reading ------------------------------------------------------
starttime <- as.character(Sys.time())
library(rutilsMH)
library(aMSE)
library(microbenchmark)
setpalette("R4")

# read data files ----------------------------------------------------

 sourcedir <- "./../../rcode2/aMSE/data-raw"
 source(filenametopath(sourcedir,"sourcer.R"))

 rundir <- "./../../rcode2/aMSEUse/run2"
 setupdirs(rundir)

 datadir
 plotdir

 ctrl <- readctrlfile(datadir,infile="control.csv")
 runname <- ctrl$runlabel
 region1 <- readregionfile(datadir,ctrl$regionfile)
 glb <- region1$globals
 constants <- readdatafile(datadir,ctrl$datafile,glb)

 setuphtml(plotdir,tabledir,runname)

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

str(regionC[[1]])
str(regionC[[4]])
str(product)

# Summarize MSY and related statistics -------------------------------
approxMSY <- findmsy(product)
approxMSY

# Some summaries -----------------------------------------------------
sapply(regionC,"[[","MSY")           # msy by population
sum(sapply(regionC,"[[","MSY"))      # total msy
sum(sapply(regionC,"[[","effB0"))    # total effective B0
sum(sapply(regionC,"[[","B0"))       # total B0
sum(sapply(regionC,"[[","effExB0"))  # total effective exploitable B0
sum(sapply(regionC,"[[","ExB0"))     # total exploitable B0



# testing the equilibrium --------------------------------------------
# This runs the dynamics for Nyrs with zero harvest rate, a constant
# larval dispersal rate (from globals), and negligible recruitment
# variation and expects the calculated components to remain constant
# within three decimal places.
 regDe <- testequil(regionC,regionD,glb)
 str(regDe)



 source(filenametopath(sourcedir,"sourcer.R"))

pop=1
pick <- findF1(pop=pop,res=product,location=TRUE)
product[pick,,pop]
findmsy(product)[1,]



plotprep(width=7,height=6,newdev=FALSE)
parset(plots=c(2,1))
xval <- findmsy(product)
numpop <- glb$numpop
plotprod(product,xname="MatB",xlab="Spawning Biomass")
for (pop in 1:numpop) abline(v=xval[pop,"MatB"],lwd=2,col=pop)
plotprod(product,xname="AnnH",xlab="Annual Harvest Rate")
for (pop in 1:numpop) abline(v=xval[pop,"AnnH"],lwd=2,col=pop)

# characterize productivity ------------------------------------------

filen <- paste0("product_",runname,".RData")
save(product,file=filenametopath(tabledir,filen))

filen <- paste0("msytable_",runname,".csv")
xval <- findmsy(product)
write.csv(xval,file = filenametopath(tabledir,filen))

filen <- filenametopath(tabledir,paste0("F1table_",runname,".csv"))
xval <- findF1(product)
write.csv(xval,file = filen)


tmp <- read.csv(filen)



# plot of Yield vs Spawning biomass
 file <- paste0("production_SpB_",ctrl$runlabel,".png")
 filename <- filenametopath(plotdir,file)
 plotprep(width=7,height=4,newdev=FALSE,filename=filename)
 xval <- findmsy(product)
 numpop <- glb$numpop
 plotprod(product,xname="MatB",xlab="Spawning Biomass t",
          ylab="Production t")
 for (pop in 1:numpop) abline(v=xval[pop,"MatB"],lwd=2,col=pop)
 legend("topright",paste0("P",1:numpop),lwd=3,col=c(1:numpop),bty="n",
        cex=1.2)
graphics.off()

 caption <- "The production curve relative to each population's spawning biomass. The vertical lines identify the Bmsy values."
 time <- as.character(Sys.time())
 cat(c(filename,caption,"production",time," \n"),file=plottabfile,sep=",",append=TRUE)


 # plot of Yield vs Annual Harvest Rate
 file <- paste0("production_AnnH_",ctrl$runlabel,".png")
 filename <- filenametopath(plotdir,file)
 plotprep(width=7,height=4,newdev=FALSE,filename=filename)
 xval <- findmsy(product)
 numpop <- glb$numpop
 plotprod(product,xname="AnnH",xlab="Annual Harvest Rate")
 for (pop in 1:numpop) abline(v=xval[pop,"AnnH"],lwd=2,col=pop)
 legend("topright",paste0("P",1:numpop),lwd=3,col=c(1:numpop),bty="n",
        cex=1.2)
 dev.off()

 caption <- "The production curve relative to the Annual Harvest Rate applied to each population. The vertical lines identify the Hmsy values."
 time <- as.character(Sys.time())
 cat(c(filename,caption,"production",time," \n"),file=plottabfile,sep=",",append=TRUE)


 # plot of Yield vs population depletion
 file <- paste0("production_Deplet_",ctrl$runlabel,".png")
 filename <- filenametopath(plotdir,file)
 plotprep(width=7,height=4,newdev=FALSE,filename=filename,verbose=FALSE)
 xval <- findmsy(product)
 numpop <- glb$numpop
 plotprod(product,xname="Deplet",xlab="Population Depletion Level")
 for (pop in 1:numpop) abline(v=xval[pop,"Deplet"],lwd=2,col=pop)
 legend("topright",paste0("P",1:numpop),lwd=3,col=c(1:numpop),bty="n",
        cex=1.2)
 dev.off()

 caption <- "The production curve relative to the depletion level of each population. The vertical lines identify the Depletion level giving rise to the MSY."
 time <- as.character(Sys.time())
 cat(c(filename,caption,"production",time," \n"),file=plottabfile,sep=",",append=TRUE)



 # plot of Yield vs population depletion but constrained to within
 # 0.2 and 0.35 levels, to illustrate nearly flat rpoduction curve
 # and more clearly identify the population depletion at MSY
 file <- paste0("production_Deplet_0.2_0.35_",ctrl$runlabel,".png")
 filename <- filenametopath(plotdir,file)
 plotprep(width=7,height=4,newdev=FALSE,filename=filename,verbose=FALSE)
 xval <- findmsy(product)
 numpop <- glb$numpop
 plotprod(product,xname="Deplet",xlab="Population Depletion Level",
          xlimit=c(0.2,0.35))
 for (pop in 1:numpop) abline(v=xval[pop,"Deplet"],lwd=2,col=pop)
 legend("topright",paste0("P",1:numpop),lwd=3,col=c(1:numpop),bty="n",
        cex=1.2)
 dev.off()

 caption <- "The production curve relative to the depletion level of each population. Here the x-axis is shortened to clarify the flatness of the production curve about the MSY points."
 time <- as.character(Sys.time())
 cat(c(filename,caption,"production",time," \n"),file=plottabfile,sep=",",append=TRUE)


# end characterize productivity --------------------------------------

 tablefile <- read.csv(plottabfile,colClasses = "character")
 class(tablefile)

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


 make_html(replist=reportlist,rundir=rundir,width=500,
           openfile=TRUE,runnotes=NULL,verbose=TRUE)














