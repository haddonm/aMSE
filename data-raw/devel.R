
library(rutilsMH)
library(aMSE)
library(microbenchmark)
datadir <- "./../../rcode2/aMSE/data-raw/"

  data(condDat)

#makeZone <- function(condDat,uplim=0.4) { #
#   condDat=condDat;  uplim=0.4
  glb <- condDat$globals
  for (i in 1:length(glb)) assign(names(glb)[i],glb[[i]])
  ls()









  doproduction <- function(regC,regD,glob,uplim=0.4) { #
# regC=regionC; regD=regionD; glob=glb; uplim=0.4
    numpop <- glb$numpop
    Nyrs <- glb$Nyrs

    imph <- seq(0.01,uplim,0.01)
    numrow <- length(imph) # hyr <- 1; year <- 2

    columns <- c("ExB","MatB","AnnH","Catch","Deplet","RelCE")
    results <- array(0,dim=c(numrow,length(columns),numpop),
                     dimnames=list(imph,columns,1:numpop))
  #  for (pop in 1:numpop) {  # pop = 1; hyr=1; year=2
      for (hyr in 1:numrow) {    # do the dynamics for numrow diff H values
        for (year in 2:Nyrs) { # always leaving yr1 the same hyr=1;year=2
          ExpB <- dyn$exploitB[year-1,pop]
          catch <- imph[hyr] * ExpB   #inpopC,glob,dyn,incatch,yr,pop
          out <- oneyear(regionC[[pop]],glb,dyn=dyn,catch,year,pop)
          dyn$exploitB[year] <- out$ExploitB
          dyn$matureB[year] <- out$MatureB
          dyn$catch[year] <- out$Catch
          dyn$harvestR[year] <- out$Harvest
          dyn$Nt[,year,pop] <- out$Nt
        }
        results[hyr,"ExB",pop] <- dyn$exploitB[Nyrs,pop]
        results[hyr,"MatB",pop] <- dyn$matureB[Nyrs,pop]
        results[hyr,"AnnH",pop] <- dyn$harvestR[Nyrs,pop]
        results[hyr,"Catch",pop] <- dyn$catch[Nyrs,pop]
        results[hyr,"Deplet",pop] <- results[hyr,"MatB",pop]/inpopC$B0
        results[hyr,"RelCE",pop] <- inpopC$popq * results[hyr,"ExB",pop]
      }
      round(results[,,1],4)

      msy <- max(results[,"Catch"])
      pick <- which(results[,"Catch"] == msy)
      msyDepl <- results[pick,"Deplet"]
      ans <- list(msy,msyDepl,results)
      names(ans) <- c("MSY","MSYDepl","Productivity")
      return(ans)
    } # end of dop

    imph <- seq(0.01,uplim,0.01)   # harvest rate for productivity
    numrow <- length(imph)
    pops <- seq(1,numpop,1)
    columns <- c("ExB","MatB","AnnH","Catch","Deplet","RelCE")
    production <- array(0,dim=c(numpop,numrow,6),
                        dimnames=list(pops,imph,columns))
    out <- doproduction(zone[[pop]],uplim=uplim)
    #  zone[[pop]]$MSY <- out$MSY
    #  zone[[pop]]$MSYDepl <- out$MSYDepl
    #   production[pop,,] <- out$Productivity


# file reading ------------------------------------------------------

library(rutilsMH)
library(aMSE)
library(microbenchmark)
setpalette("R4")

# read data files ----------------------------------------------------
datadir <- "./../../rcode2/aMSE/data-raw/"
ctrlfile <- "control.csv"
source(filenametopath(datadir,"sourcer.R"))

ctrl <- readctrlfile(datadir,ctrlfile)
reg1 <- readregionfile(datadir,ctrl$regionfile)
glb <- reg1$globals
constants <- readdatafile(datadir,ctrl$datafile,glb)

# Define the Zone ----------------------------------------------------
source(filenametopath(datadir,"sourcer.R"))
ans <- makeregionC(reg1,constants)
regionC <- ans$regionC
popdefs <- ans$popdefs
ans <- makeregion(glb,regionC)

regionC <- ans$regionC  # region constants
regionD <- ans$regionD  # region dynamics

#str(regionC[[1]])
# if larvdisp > 0.0 search for equilibrium----------------------------

if (glb$larvdisp > 0.0) ans <- findunfished(regionC,regionD,glb)

inH <- rep(0.0,glb$numpop)
regLC <- ans$regionC  # region constants
regLD <- ans$regionD  # region dynamics

# testing the equilibrium --------------------------------------------
# regDe <- testequil(regLC,regLD,glb,inH)
# regDe$matureB
# regDe$exploitB
# regLD$cpue
#
#
# npop <- glb$numpop
# Nclass <- glb$Nclass
# Nyrs <- glb$Nyrs
# catch <- rep(0.0,npop)
#
# for (yr in 2:Nyrs) {
#   regLD <- oneyearD(regC=regLC,regD=regLD,Ncl=Nclass,incatch=catch,
#                       year=yr,sigmar=1e-08,npop=npop,deltarec=reg1$larvdisp)
# }
# str(regLD)



regC <- regLC;  regD=regLD;  glob=glb; initdepl=0.5; uplim=0.3

numpop <- glob$numpop
Nclass <- glob$Nclass
Nyrs <- glob$Nyrs
larvdisp <- glob$larvdisp
initH <- seq(0.01,uplim,0.005)
nH <- length(initH)
columns <- c("ExB","MatB","AnnH","Catch","Deplet","RelCE")
results <- array(0,dim=c(nH,6,numpop),dimnames=list(initH,columns,1:numpop))
exb <- regLD$exploitB[1,]

for (aH in 1:nH) { # aH=1 ; yr=2
# regD <- restart(oldregD=regD,nyrs=Nyrs,npop=numpop,N=Nclass,zero=TRUE)
 regD <- runthreeH(regC=regC, regD=regD, glob=glb, inHarv=rep(initH[aH],numpop))
  results[aH,"ExB",] <- regD$exploitB[1,]
  results[aH,"MatB",] <- regD$matureB[1,]
  results[aH,"AnnH",] <- regD$harvestR[1,]
  results[aH,"Catch",] <- regD$catch[1,]
  results[aH,"Deplet",] <- regD$deplsB[1,]
  results[aH,"RelCE",] <- regD$cpue[1,]
} # end of yr loop

plotprep(width=7,height=6,newdev=FALSE)
parset(plots=c(2,1))
plot(results[,"AnnH",1],results[,"Catch",1],type="l",lwd=2,col=1,ylim=c(0,160))
for (i in 2:numpop) lines(results[,"AnnH",i],results[,"Catch",i],lwd=2,col=i)
for (i in 1:numpop) {
  catch <- results[,"Catch",i]
  pick <- findF1(i,results)
  abline(v=results[pick,"AnnH",i],lwd=2,col=i)
  print(c(catch[which.max(catch)],catch[pick]))
}
plot(results[,"AnnH",1],results[,"Catch",1],type="l",lwd=2,col=1,ylim=c(0,160))
for (i in 2:numpop) lines(results[,"AnnH",i],results[,"Catch",i],lwd=2,col=i)
for (i in 1:numpop) {
  catch <- results[,"Catch",i]
  pick <- findF1(i,results)
  abline(v=results[which.max(catch),"AnnH",i],lwd=2,col=i)
}

findF1(1,results,location=FALSE)








