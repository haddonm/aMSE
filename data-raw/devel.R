# file reading ------------------------------------------------------
starttime <- as.character(Sys.time())
library(rutilsMH)
library(aMSE)
library(microbenchmark)


# read data files ----------------------------------------------------
 resdir <- "./../../rcode2/aMSEUse/out/run1"
 dirExists(resdir,make=TRUE,verbose=TRUE)
 # You now need to ensure that there is a control.csv, reg1smu2pop6.csv
 # and region1.csv file in the data directory
 ctrl <- checkresdir(resdir)
 runname <- ctrl$runlabel
 region1 <- readregionfile(resdir,ctrl$regionfile)
 glb <- region1$globals
 constants <- readdatafile(glb$numpop,resdir,ctrl$datafile)

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
biology_plots(resdir, runname, glb, regionC)
numbersatsize(resdir, runname, glb, regionD)

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
   starttime=starttime,endtime=endtime,
   regionC=regionC, regionD=regionD, product=product,
   glb=glb,constants=constants
 )
 str(reportlist,max.level = 1)

 runnotes <- paste0("The results presented here relate to the included data-sets testregC, ",
                    "testregD, and product. They are for a region made up of 2 SMU and 6 population. ",
                    "These results are currently under development and there are many more needed yet.")

#  source(filenametopath(sourcedir,"sourcer.R"))
 make_html(replist=reportlist,resdir=resdir,width=500,
           openfile=TRUE,runnotes=runnotes,verbose=FALSE)



# end of run ---------------------------------------------------------

# outline a real run--------------------------------------------------

 npop <- glb$numpop
 Nc <- glb$Nclass
 nyrs <- glb$Nyrs
 larvdisp <- glb$larvdisp
 catch <- 300.0

 for (yr in 2:nyrs) {
       regD <- oneyearD(regC=regionC,regD=regionD,Ncl=Nc,
                        inHt=inHarv,year=yr,sigmar=1e-08,npop=npop,
                        deltarec=larvdisp)

 }



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


