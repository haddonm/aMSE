# file reading ------------------------------------------------------
starttime <- as.character(Sys.time())
library(rutilsMH)
library(aMSE)
library(microbenchmark)

#Rprof()
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
# Rprof(NULL)
# outprof <- summaryRprof()
# outprof
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

 oneyearcat <- function(inpopC,inNt,Nclass,incat,yr) {  #
   # yr=2; pop=2; inpopC=regC[[pop]]; inNt=regD$Nt[,yr-1,pop];
   # Nclass=glb$Nclass; inH=0.05;
   MatWt <- inpopC$MatWt/1e06
   SelectWt <- inpopC$SelWt[,yr]/1e06
   selyr <- inpopC$Select[,yr]
   Ne <- numeric(Nclass)
   Cat <- numeric(Nclass)
   Os <- exp(-inpopC$Me/2)
   #  MatureB <- sum(MatWt*inNt)
   NumNe <- (Os * (inpopC$G %*% inNt))
   ExploitB <- sum(SelectWt * NumNe) #SelectWt=Select*WtL
   oldExpB <- ExploitB   # ExploitB after growth and 0.5NatM
   estH <- min(incat/ExploitB,0.8) # no more than 0.8 harvest rate
   Fish <- 1-(estH*selyr)
   newNt <- (Os * (Fish * NumNe)) #+ Rec # Nt - catch - 0.5M, and + Rec
   Cat <- (estH*selyr) * NumNe  #numbers at size in the catch
   ExploitB <- (sum(SelectWt * newNt) + oldExpB)/2.0 #av start and end
   MatureB <- sum(MatWt*newNt) #+ MatBC
   Catch <- sum(inpopC$WtL*Cat)/1e06
   Harvest <- Catch/ExploitB  # uses average of the start and end
   ce <- inpopC$popq * ExploitB * 1000.0  #ExploitB
   ans <- list(ExploitB,MatureB,Catch,Harvest,newNt,ce,Cat)
   names(ans) <- c("ExploitB","MatureB","Catch","Harvest","Nt","ce",
                   "CatchN")
   return(ans)
 } # End of oneyear





 storeregC <- regionC
 storeregD <- regionD

 regionC <- storeregC
 regionD <- storeregD
 npop <- glb$numpop
 Nc <- glb$Nclass
 nyrs <- glb$Nyrs
 larvdisp <- glb$larvdisp
 catch <- 360.0
 B0 <- getvar(regionC,"B0") #sapply(regionC,"[[","B0")
 totB0 <- sum(B0)
 prop <- B0/totB0
 catchpop <- catch * prop

 for (yr in 2:nyrs) {
       regionD <- oneyearC(regC=regionC,regD=regionD,Ncl=Nc,
                        catchp=catchpop,year=yr,sigmar=1e-08,npop=npop,
                        movem=glb$move)
  #  cedecline <- regionD$cpue[yr,]/regionD$cpue[1,]
  #  cewt <- cedecline/sum(cedecline)
  #  catchpop <- catch * cewt
 }
str(regionD)
regionD$harvestR[40,]
regionD$exploitB[40,]
regionD$matureB[40,]
regionD$deplsB[40,]
regionD$depleB[40,]
regionD$cpue[1,]
regionD$cpue[40,]

regionD$cpue[40,]/regionD$cpue[1,]

Nt <- regionD$Nt
mids <- glb$midpts
plotprep(width=8, height=9,newdev=FALSE)
parset(plots=c(3,2),margin=c(0.3,0.3,0.05,0.05),outmargin=c(1,1,0,0))
for (pop in 1:numpop) {
  plot(mids[5:105],Nt[5:105,1,pop],type="l",lwd=2,panel.first=grid())
  lines(mids[5:105],Nt[5:105,20,pop],lwd=2,col=3)
  lines(mids[5:105],Nt[5:105,40,pop],lwd=2,col=2)
  abline(v=132,col="grey")
  mtext(paste0("pop",pop),side=3,outer=FALSE,line=-1.1)
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
















