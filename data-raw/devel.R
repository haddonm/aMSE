
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
ctrl <- readctrlfile(datadir,ctrlfile)
reg1 <- readregionfile(datadir,ctrl$regionfile)
glb <- reg1$globals
constants <- readdatafile(datadir,ctrl$datafile,glb)

# Define the Zone ----------------------------------------------------
ans <- makeregionC(reg1,constants)
regionC <- ans$regionC
popdefs <- ans$popdefs
ans <- makeregion(reg1$globals,regionC)

regionC <- ans$regionC  # region constants
regionD <- ans$regionD  # region dynamics


# if larvdisp > 0.0 search for equilibrium----------------------------

findunfished <- function(regC,regD,glob) {
#  regC=regionC; regD=regionD; glob=glb
  numpop <- glob$numpop
  catch <- rep(0.0,numpop)
  regD <- runthree(regC,regD,glob,catch)

  for (pop in 1:numpop) {
    regC[[pop]]$B0 <- regD$matureB[1,pop]
    regC[[pop]]$ExB0 <- regD$exploitB[1,pop]
  }

  Nclass <- glob$Nclass
  Nyrs <- glob$Nyrs
  larvdisp <- glob$larvdisp


  regD <- restart(oldregD=regD,nyrs=Nyrs,npop=numpop,N=Nclass)



}



if (glb$larvdisp > 0.0) findunfished(regionC,regionD,glb)





for (yr in 2:Nyrs) {
  regionD <- oneyearD(regC=regionC,regD=regionD,Ncl=Nclass,incatch=catch,
                      year=yr,sigmar=1e-08,npop=npop,deltarec=reg1$larvdisp)
}
str(regionD)

regionD$matureB

plotprep(width=7,height=5,newdev=FALSE)
plot(1:40,regionD$deplsB[,1],type="l",lwd=2,ylim=c(0.950,1.1),panel.first=grid())
for (i in 2:npop) lines(1:40,regionD$deplsB[,i],lwd=2,col=i)




finddepletion <- function(regC,regD,glob,initdepl,uplim=0.45) {
  numpop <- glob$numpop
  Nclass <- glob$Nclass
  Nyrs <- glob$Nyrs
  initH <- seq(0.05,uplim,0.01)
  nH <- length(initH)


} # end of findequil





regionD2 <- restart(oldregD=regionD,nyrs=Nyrs,npop=npop,N=Nclass)

for (yr in 2:Nyrs) {
  regionD2 <- oneyearD(regC=regionC,regD=regionD2,Ncl=Nclass,incatch=catch,
                      year=yr,sigmar=1e-08,npop=npop,deltarec=reg1$larvdisp)
}

plotprep(width=7,height=5,newdev=FALSE)
plot(1:40,regionD2$deplsB[,1],type="l",lwd=2,ylim=c(0.950,1),panel.first=grid())
for (i in 2:npop) lines(1:40,regionD2$deplsB[,i],lwd=2,col=i)

regionD$matureB[1,]
regionD$matureB[Nyrs,]
regionD3$matureB[1,]
regionD4$matureB[1,]
