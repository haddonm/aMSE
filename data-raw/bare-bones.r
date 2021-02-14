# Used in readme.Rmd to illustrate the running of aMSE

options("show.signif.stars"=FALSE,
        "stringsAsFactors"=FALSE,
        "max.print"=50000,
        "width"=240)

# use dohistoricC  --------------------------------------------------------
starttime <- (Sys.time())
library(aMSE)
library(rutilsMH)
library(makehtml)
library(knitr)
# Obviously you should modify the resdir to suit your own computer
if (dir.exists("c:/Users/User/DropBox")) {
  ddir <- "c:/Users/User/DropBox/A_code/"
  } else {
  ddir <- "c:/Users/Malcolm/DropBox/A_code/"
}
#resdir <- paste0(ddir,"aMSEUse/conddata/generic2")
resdir <- paste0(ddir,"aMSEUse/conddata/generic")
dirExists(resdir,make=TRUE,verbose=TRUE)

# this should be in resdir, but during development is in data-raw
source(paste0(ddir,"aMSE/data-raw/","TasmanianHS.R"))

#data(zone)
zone <- makeequilzone(resdir,"control2.csv") # normally would read in a file
  equiltime <- (Sys.time())
  # condition on the historic catches
zoneDD <- dohistoricC(zone$zoneD,zone$zoneC,glob=zone$glb,zone$zone1$condC)

  print(equiltime - starttime)

  propD <- getzoneprops(zone$zoneC,zoneDD,zone$glb,year=47)
  round(propD,3)

# Do the replicates ------------------------------------------------------------

midtime <- (Sys.time())
  print(midtime - starttime)

glb <- zone$glb
ctrl <- zone$ctrl
ctrl$sigmaF <- 0.05   #MODIFY INPUT FILE AND CTRLREAD TO INCLUDE sigmaF
zone1 <- zone$zone1
projC <- zone1$projC
condC <- zone1$condC
zoneC <- zone$zoneC

#   source(paste0(ddir,"aMSE/data-raw/","TasmanianHS.R"))


cmcda <- calibrateMCDA(histCE=condC$histCE, saunames=zone1$SAUnames,
                       hsargs=hsargs,sauindex=glb$sauindex,projyrs=projC$projyrs)

str(cmcda)
pms <- cmcda$pms
multTAC <- cmcda$yrmultTAC

out <- prepareprojection(zone1,zoneC,glb,zoneDD)

zoneDP <- out$zoneDP
projC <- out$projC
zoneCP <- out$zoneCP


endtime <- Sys.time()
print(endtime - midtime)


# endcatch the final year's catches from the fishery conditioned zone.
#     These are required to generate the aspiration catches by population for
#     the first year of the projections
#
#
# develop projections-----------------------------------------------------------

zoneCP=zoneCP
zoneDP=zoneDP
glob=glb
ctrl=ctrl
projyrs=projC$projyrs
applyHS=mcdahcr
hsargs=hsargs
projpms=cmcda
# get important constants

npop <- glob$numpop
nsau <- glob$nSAU
Ncl <- glob$Nclass
yrce <- nrow(condC$histCE)
nyrs <- projyrs
endyr <- nrow(zoneDD$catch)
movem <- glob$move
reps <- ctrl$reps
sauindex <- glob$sauindex
matb <- numeric(npop)
origTAC <- sum(zoneDD$catch[yrce,]) # mean sum of catches in last year
if (ctrl$randseedP > 0) set.seed(ctrl$randseedP) # set random seed if desired
endcatch <- tapply(zoneDD$catch[endyr,],sauindex,sum,na.rm=TRUE)
acatch <- endcatch * multTAC[yrce,]

# for (iter in 1:reps) {
iter <- 1
TAC <- sum(acatch)
#for (year in 2:nyrs) { # iter=1; year=11
 year <- 1
 exb <- zoneDD$exploitB[endyr,]
 inN <- zoneDD$Nt[,endyr,]
 zoneCC <- zoneCP
 catchsau=acatch
 sigmar <- ctrl$withsigR # needed to add recruitment variation
 sigmaF <- ctrl$withsigB


 oneyearsauC2 <- function(zoneCC,exb,inN,catchsau,year,sigmar,Ncl,
                         sauindex,movem,sigmaF) {
   popC <- imperr(catchsau,exb,sauindex,sigmaF)
   matb <- numeric(npop)
   ans <- vector("list",npop)
   for (popn in 1:npop) {  # popn=1
     ans[[popn]] <- oneyearcat(inpopC=zoneCC[[popn]],inNt=inN[,popn],
                       Nclass=Ncl,incat=popC[popn],yr=year)
   }
   dyn <- sapply(ans,"[[","vect")
   steep <- getvect(zoneCC,"steeph") #sapply(zoneC,"[[","popdef")["steeph",]
   r0 <- getvar(zoneCC,"R0") #sapply(zoneC,"[[","R0")
   b0 <- getvar(zoneCC,"B0") #sapply(zoneC,"[[","B0")
   recs <- oneyearrec(steep,r0,b0,dyn["matureb",],sigR=sigmar)
   newrecs <- as.numeric(movem %*% recs)
   deplsB <- dyn["matureb",]/b0
   depleB <- dyn["exploitb",]/getvar(zoneCC,"ExB0")
   dyn <- rbind(dyn,newrecs,deplsB,depleB)
   NaL <-  sapply(ans,"[[","NaL")
   NaL[1,] <- newrecs
   catchN <-  sapply(ans,"[[","catchN")
   return(list(dyn=dyn,NaL=NaL,catchN=catchN))
 } # end of oneyearsauC2

out <- oneyearsauC2(zoneCP,exb=zoneDD$exploitB[endyr,],inN=zoneDD$Nt[,endyr,],
                    catchsau=acatch,year=1,sigmar=ctrl$withsigR,Ncl=glb$Nclass,
                    sauindex=sauindex,movem=glb$move,sigmaF=ctrl$withsigR)

zoneDP$exploitB[year,,iter] <- out$dyn["exploitb",]
zoneDP$matureB[year,,iter] <- out$dyn["matureb",]
zoneDP$catch[year,,iter] <- out$dyn["catch",]
zoneDP$catsau[year,,iter] <- acatch
zoneDP$harvestR[year,,iter] <- out$dyn["catch",]/out$dyn["exploitb",]
zoneDP$cpue[year,,iter] <- out$ce
zoneDP$Nt[,year,,iter] <- out$Nt
zoneDP$catchN[,year,,iter] <- out$CatchN




 catbysau <- catchbysau(inexpB=zoneDP$exploitB[(year - 1),,iter],sauindex,
                        TAC,sigmaF=sigmaF) # no error initially
 multh <- apply(zoneDP$cesau[1:(year-1),,iter],2,applyHS,yr=(year-1)) # apply mcdahcr
 TAC <- sum(catbysau * multh)
 divererr <- sauexpB * exp(rnorm(nsau,mean=0,sd=ctrl$withsigB))
 catbysau <- TAC * (divererr/sum(divererr)) # currently no error on TAC
 catbypop <- catbysau[sauindex] * (inexpB/sauexpB[sauindex]) # no error on pops
 for (popn in 1:npop) { # year=11; iter=1; pop=1
   out <- oneyearcat(inpopC=zoneCP[[popn]],inNt=zoneDP$Nt[,year-1,popn,iter],
                     Nclass=Ncl,incat=catbypop[popn],yr=year)
   zoneDP$exploitB[year,popn,iter] <- out$ExploitB
   zoneDP$matureB[year,popn,iter] <- out$MatureB
   zoneDP$catch[year,popn,iter] <- out$Catch
   zoneDP$harvestR[year,popn,iter] <- out$Harvest
   zoneDP$cpue[year,popn,iter] <- out$ce
   zoneDP$Nt[,year,popn,iter] <- out$Nt
   zoneDP$catchN[,year,popn,iter] <- out$CatchN
   matb[popn] <- out$MatureB
 } # pop
 steep <- getvect(zoneCP,"steeph")
 r0 <- sapply(zoneCP,"[[","R0")
 b0 <- sapply(zoneCP,"[[","B0")
 recs <- oneyearrec(steep,r0,b0,matb,sigR=sigmar)
 newrecs <- movem %*% recs
 zoneDP$recruit[year,,iter] <- newrecs
 zoneDP$Nt[1,year,,iter] <- newrecs
 zoneDP$deplsB[year,,iter] <- zoneDP$matureB[year,,iter]/b0
 zoneDP$depleB[year,,iter] <- zoneDP$exploitB[year,,iter]/sapply(zoneCP,"[[","ExB0")
 saucatch[year,,iter] <- tapply(zoneDP$catch[year,,iter],sauindex,sum,na.rm=TRUE)
 wts <- zoneDP$catch[year,,iter]/(saucatch[year,sauindex,iter])
 saucpue[year,,iter] <- tapply((zoneDP$cpue[year,,iter] * wts),sauindex,sum,na.rm=TRUE)















saunames <- zone$zone1$SAUnames
sauindex <- glb$sauindex
pyrs <- projC$projyrs
B0 <- tapply(sapply(zone$zoneC,"[[","B0"),sauindex,sum)
exB0 <- tapply(sapply(zone$zoneC,"[[","ExB0"),sauindex,sum)


# l -----------------------------------------------------------------------------
# depleteSAU and histCE --------------------------------------------------------
options("show.signif.stars"=FALSE,
        "stringsAsFactors"=FALSE,
        "max.print"=50000,
        "width"=240)

starttime <- (Sys.time())
library(aMSE)
library(rutilsMH)
library(makehtml)
library(knitr)
# Obviously you should modify the resdir to suit your own computer
if (dir.exists("c:/Users/User/DropBox")) {
  ddir <- "c:/Users/User/DropBox/A_code/"
} else {
  ddir <- "c:/Users/Malcolm/DropBox/A_code/"
}
#resdir <- paste0(ddir,"aMSEUse/conddata/generic2")
resdir <- paste0(ddir,"aMSEUse/conddata/generic")
dirExists(resdir,make=TRUE,verbose=TRUE)
# this should be in resdir, but during development is in data-raw
source(paste0(ddir,"aMSE/data-raw/","TasmanianHS.R"))
# source(paste0(resdir,"/TasmanianHS.R"))

zone <- makeequilzone(resdir,"control2.csv") # normally would read in a file
equiltime <- (Sys.time())


glb <- zone$glb
zone1 <- zone$zone1
projC <- zone1$projC
condC <- zone1$condC


# histCE <- zone$zone1$condC$histCE; saunames=zone$zone1$SAUnames

# ans <- calibrateMCDA(zoneDD$catch,zoneDD$cpue,projC$inityrs:glb$Nyrs,
#                      nsau=glb$nSAU,sauindex=glb$sauindex,pyrs=projC$projyrs)

deplinit <- c(0.23,0.23,0.23,0.23,0.23,0.23,0.23,0.23)
zoneDD <- depleteSAU(zone$zoneC,zone$zoneD,glob=glb,initdepl=deplinit,zone$product)
propD <- getzoneprops(zone$zoneC,zoneDD,glb,year=1)
round(propD,3)

ctrl <- zone$ctrl
ctrl$withsigR <- 0.1 # all starting points are effectively identical
ans <- prepareprojection(zone1,zone$zoneC,glb,zoneDD,ctrl)
str(ans$zoneDP,max.level = 1)

cmcda <- calibrateMCDA(histCE=condC$histCE,saunames=zone1$SAUnames,
                       hsargs=hsargs,zoneDP=ans$zoneDP,glb$sauindex)

str(cmcda,max.level = 2)
zoneDP <- cmcda$zoneDP
histpms <- cmcda$pms


projpms <- makeprojpm(cmcda$pms,ctrl$reps,projC$projyrs)

ctrl$withsigR <- 0.3

source(paste0(ddir,"aMSE/data-raw/","TasmanianHS.R"))

# Now do Projection -------------------------------------------------------









zoneDP <- cmcda$zoneDP

tac <- colSums(zoneDP$catch[1,,])

plotprep(width=7, height=7,newdev=FALSE)
hist(zoneDP$deplsB[1,,],main="")

#Rprof()
# this should be in resdir, but during development is in data-raw


mseproj <- doprojection(zoneC,zoneDP,glb,ctrl,projC$projyrs,applyHS=mcdahcr,
                        hsargs=hsargs)
sauzoneDP <- asSAU(mseproj,sauindex,saunames,B0,exB0)
zoneproj <- aszone(sauzoneDP,zoneCP)
#Rprof(NULL)
endtime <- (Sys.time())
print(endtime - midtime)

#summaryRprof()

pyrs <- projC$projyrs + projC$inityr
reps <- ctrl$reps
yrs <- 1:pyrs


plotprep(width=9, height=8,newdev=FALSE)
parset(plots=c(3,2))
CIR <- plotzoneproj(zoneproj$zoneR/1000,reps,yrs,"Recruitment",addqnts=TRUE,miny=2000)
CIC <- plotzoneproj(zoneproj$zoneC,reps,yrs,"Catches t",addqnts=TRUE,miny=0)
CIH <- plotzoneproj(zoneproj$zoneH,reps,yrs,"Harvest Rate",addqnts=TRUE)
CIH <- plotzoneproj(zoneproj$zonece,reps,yrs,"CPUE",addqnts=TRUE)
CIsBD <- plotzoneproj(zoneproj$zonedeplsB,reps,yrs,"SpB Depletion",addqnts=TRUE)
CIeB <- plotzoneproj(zoneproj$zoneeB,reps,yrs,"Exploitable Biomass",addqnts=TRUE)



plotC <- function(nsau,saunames,reps,projyrs,plts=c(4,2)) {
  return(list(nsau=nsau,saunames=saunames,reps=reps,projyrs=projyrs,plts=plts))
}

saunames <- zone$zone1$SAUnames
nyrs <- projC$projyrs + projC$inityr
reps <- ctrl$reps
label <- "Catches t"
yrs <- 1:nyrs

plotconst <- plotC(nsau=length(saunames),saunames,reps,pyrs,plts=c(4,2))

plotprep(width=7,height=8,newdev=FALSE)
plotproj(sauzoneDP$cpue,"CPUE",plotconst,vline=projC$inityrs,addqnts=TRUE)

plotproj(sauzoneDP$saudeplsB,"Spawning Biomass Depletion",plotconst,vline=projC$inityrs,addqnts=TRUE)

plotproj(sauzoneDP$recS/1000,"Recruitment '000s",plotconst,vline=projC$inityrs,addqnts=TRUE)

plotproj(sauzoneDP$catS,"Catches t",plotconst,vline=projC$inityrs,addqnts=TRUE)

plotproj(sauzoneDP$harvS,"Harvest Rate",plotconst,vline=projC$inityrs,addqnts=TRUE)



plotprod(zone$product,xname="MatB",xlab="Spawning Biomass t",
         ylab="Production t")


# equilibrium zone characterization---------------------------------------------
# resfile <- setuphtml(resdir)# prepare to save and log results
#
# plotproductivity(resdir,product,glb)
# biology_plots(resdir, glb, zoneC)
# numbersatsize(resdir, glb, zoneD)
#
# endtime <- as.character(Sys.time())
#
# reportlist <- list(starttime=starttime,endtime=endtime,
#                    zoneC=zoneC, zoneD=zoneD, product=product,
#                    glb=glb,constants=constants)
# runnotes <- "This is a bare-bones example."
# # If you unhash this component it will generate a local website inside
# # resdir and open it so you can see the results so far.
# make_html(replist=reportlist,resdir=resdir,width=500,
#           openfile=TRUE,runnotes=runnotes,verbose=FALSE,
#           packagename = "aMSE",htmlname="testrun")

# prepare the HS --------------------------------------------------------------

# if (projC$HS == "constantCatch") {
#   hsFunc <- constCatch
#   inTAC <- projC$HSdetail
# }
# if (projC$HS == "MCDA") {
#   hsFunc <- doMCDA
#   mcdafile <- projC$HSdetail
#   optCE <- MCDAdata(resdir,mcdafile,zone1$SAUnames)
# }
#
#  zoneCP=zoneCP;zoneDP=zoneDR;glob=glb;ctrl=ctrl;projyrs=projC$projyrs;inityrs=projC$inityrs;



condC <- zone$zone1$condC

histce <- condC$histCE

cmcda <- calibrateMCDA(histce,saunames)







