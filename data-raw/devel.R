# file reading ------------------------------------------------------
starttime <- as.character(Sys.time())
library(rutilsMH)
library(aMSE)
library(microbenchmark)

Rprof()
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

 out <- setupregion(constants, region1)
 regionC <- out$regionC
 regionD <- out$regionD
 product <- out$product
 glb <- out$glb
Rprof(NULL)
outprof <- summaryRprof()
outprof
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
 B0 <- getvar(regionC,"B0") #sapply(regionC,"[[","B0")
 totB0 <- sum(B0)
 prop <- B0/totB0
 catchpop <- catch * prop

 oneyearC <- function(regC,regD,Ncl,catchp,year,sigmar,npop,deltarec) {
   exB <- regD$exploitB[(year-1),]
   matb <- numeric(npop)
   for (popn in 1:npop) {  # year=2
     out <- oneyear(inpopC=regC[[popn]],inNt=regD$Nt[,year-1,popn],
                    Nclass=Ncl,inH=inHt[popn],yr=year)
     regD$exploitB[year,popn] <- out$ExploitB
     regD$matureB[year,popn] <- out$MatureB
     regD$catch[year,popn] <- out$Catch
     regD$harvestR[year,popn] <- out$Harvest
     regD$cpue[year,popn] <- out$ce
     regD$Nt[,year,popn] <- out$Nt
     regD$catchN[,year,popn] <- out$CatchN
     matb[popn] <- out$MatureB
   }
   steep <- getvect(regC,"steeph") #sapply(regC,"[[","popdef")["steeph",]
   r0 <- getvar(regionC,"R0") #sapply(regC,"[[","R0")
   b0 <- getvar(regionC,"B0") #sapply(regC,"[[","B0")
   recs <- oneyearrec(steep,r0,b0,matb,sigR=sigmar)
   newrec <- driftrec(recs,deltarec)
   regD$recruit[year,] <- newrec
   regD$Nt[1,year,] <- newrec
   regD$deplsB[year,] <- regD$matureB[year,]/b0
   regD$depleB[year,] <- regD$exploitB[year,]/getvar(regionC,"ExB0") #
   return(regD)
 } # end of oneyearC   round(regD$Nt[,year,])






 for (yr in 2:nyrs) {
       regD <- oneyearD(regC=regionC,regD=regionD,Ncl=Nc,
                        inHt=inHarv,year=yr,sigmar=1e-08,npop=npop,
                        movem=glb$move)

 }



 # Some summaries ----------------------------------------------------
 getvar(regionC,"MSY") #sapply(regionC,"[[","MSY")           # msy by population
 sum(getvar(regionC,"MSY"))      # total msy
 sum(getvar(regionC,"B0"))       # total B0
 sum(getvar(regionC,"ExB0"))     # total exploitable B0



 filen <- paste0("product_",runname,".RData")
 save(product,file=filenametopath(tabledir,filen))

 filen <- paste0("msytable_",runname,".csv")
 xval <- findmsy(product)
 write.csv(xval,file = filenametopath(tabledir,filen))

 filen <- filenametopath(tabledir,paste0("F1table_",runname,".csv"))
 xval <- findF1(product)
 write.csv(xval,file = filen)


 tmp <- read.csv(filen)

# barebones ----------------------------------------------------------

 data(constants)
 data(region1)
 ans <- makeregionC(region1,constants) # classical equilibrium
 regionC <- ans$regionC
 glb <- ans$glb
 ans <- makeregion(glb,regionC) # now make regionD
 regionC <- ans$regionC  # region constants
 regionD <- ans$regionD
 numpop <- glb$numpop
 harvest <- rep(0.2,numpop)
 regionD <- oneyearD(regC=regionC,regD=regionD,Ncl=glb$Nclass,
                     inHt=harvest,year=2,sigmar=1e-06,npop=numpop,
                     movem=glb$move)
 str(regionD)
 round(regionD$catchN[60:105,1:5,1],1)

 # end barebones------------------------------------------------------

 # test speed doproduction



 microbenchmark(
   steep1 <- getvect(regionC,"steeph"),
   steep <- sapply(regC,"[[","popdef")["steeph",],
   times=100
 )


















