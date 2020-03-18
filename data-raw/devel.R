
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


  oneyear <- function(inpopC,glob,dyn,incatch,yr,pop) {  #
  #  inpopC=regionC[[1]]; glob=glb; dyn=regionD; incatch=0.0; yr=2; pop=1
    MatWt <- inpopC$MatWt/1e06
    SelectWt <- inpopC$SelWt[,yr]/1e06
    selyr <- inpopC$Select[,yr]
    inNt <- dyn$Nt[,yr-1,pop]
    Nclass <- glob$Nclass
    Ne <- numeric(Nclass)
    Cat <- numeric(Nclass)
    Os <- exp(-inpopC$Me/2)
    MatureB <- sum(MatWt*inNt)
    NumNe <- (Os * (inpop$G %*% inNt))
    ExploitB <- sum(SelectWt * NumNe) #SelectWt=Select*WtL
    oldExpB <- ExploitB   # ExploitB after growth and 0.5NatM
    Ht <- incatch/ExploitB
    if (Ht > 0.9) {
      Ht <- 0.9
      warning(paste0("Harvest rates > 0.9 in year ",yr," in pop ",
                     (inpop$popdef["block"])))
    }
    Fish <- 1-(Ht*selyr)
    newNt <- (Os * (Fish * NumNe)) #+ Rec # Nt - catch - 0.5M, and + Rec
    Cat <- (Ht*selyr) * NumNe  #numbers at size in the catch
    ExploitB <- sum(SelectWt * newNt)
    MatureB <- sum(MatWt*newNt) #+ MatBC
    Catch <- sum(inpopC$WtL*Cat)/1e06
    Harvest <- (2.0 * Catch)/(oldExpB + ExploitB)  # average of the start and end
    ce <- inpopC$popq * ((oldExpB + ExploitB)/2) * 1000.0  #ExploitB
    ans <- list(ExploitB,MatureB,Catch,Harvest,newNt,ce,Cat)
    names(ans) <- c("ExploitB","MatureB","Catch","Harvest","Nt","ce",
                    "CatchN")
    return(ans)
  } # End of oneyear




dooneY <- function(izone,incatch,year,sigmar,npop,deltarec) {
#  izone=zone; incatch=0.0;  year=2; sigmar=1e-08; npop=6 ; deltarec=0.02
  matb <- numeric(npop)
  for (popn in 1:npop) {
    out <- oneyear(inpop=izone[[popn]],incatch=0.0,yr=year)
    izone[[popn]]$ExploitB[year] <- out$ExploitB
    izone[[popn]]$MatureB[year] <- out$MatureB
    izone[[popn]]$Catch[year] <- out$Catch
    izone[[popn]]$HarvestR[year] <- out$Harvest
    izone[[popn]]$Nt[,year] <- out$Nt
    matb[popn] <- out$MatureB
  }
  steep <- sapply(izone,"[[","popdef")["steeph",]
  r0 <- sapply(izone,"[[","R0")
  b0 <- sapply(izone,"[[","B0")
  recs <- oneyearrec(steep,r0,b0,matb,sigR=sigmar)
  newrec <- driftrec(recs,deltarec)
  for (popn in 1:npop) izone[[popn]]$Nt[1,year] <- newrec[popn]
  return(izone)
}





 # doproduction <- function(inpop,uplim=0.4) { #

  out <- makeregionC(condDat)
  regionC <- out$regionC
  popdefs <- out$popdefs

  out2 <- makeregion(glb=glb,regC=regionC)
  regionC <- out2$regionC
  regionD <- out2$regionD

  glb <- glb
  dyn <- regionD
  inpop=regionC[[1]]
  uplim=0.4
  pop=1
  year=2
  hyr=1

    imph <- seq(0.01,uplim,0.01)
    numrow <- length(imph) # hyr <- 1; year <- 2
    numpop <- glb$numpop
    Nyrs <- glb$Nyrs
    columns <- c("ExB","MatB","AnnH","Catch","Deplet","RelCE")
    results <- array(0,dim=c(numrow,6,numpop),dimnames=list(imph,columns,1:numpop))
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
datadir <- "./../../rcode2/aMSE/data-raw/"

ctrlfile <- "control.csv"

ctrl <- readctrlfile(datadir,ctrlfile)

reg1 <- readregionfile(datadir,ctrl$regionfile)





datadir <- "./../../rcode2/aMSE/data-raw/"
ctrlfile <- "control.csv"
ctrl <- readctrlfile(datadir,ctrlfile)
reg1 <- readregionfile(datadir,ctrl$regionfile)
popdefs <- readdatafile(datadir,ctrl$datafile,reg1$globals)
print(popdefs)






