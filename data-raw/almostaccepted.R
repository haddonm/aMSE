



options("show.signif.stars"=FALSE,
        "stringsAsFactors"=FALSE,
        "max.print"=50000,
        "width"=240)
suppressPackageStartupMessages({
  library(aMSE)
  library(TasHS)
  library(hutils)
  library(hplot)
  library(makehtml)
  library(knitr)
})
dropdir <- getDBdir()
prefixdir <- paste0(dropdir,"A_codeUse/aMSEUse/scenarios/")

startime <- Sys.time()
postfixdir <- "HS856"
verbose <- TRUE
rundir <- filenametopath(prefixdir,postfixdir)
controlfile <- paste0("control",postfixdir,".csv")
outdir <- "C:/aMSE_scenarios/M15h7L75/"
confirmdir(rundir)
confirmdir(outdir)


hsargs <- list(mult=0.1,
               wid = 4,
               targqnt = 0.55,
               maxtarg = c(150,150,150,150,150,150,150,150),
               pmwts = c(0.65,0.25,0.1),
               hcr = c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2),
               startCE = 1992,
               endCE = 2012)






# naked do_MSE without producing all the plots etc

postfixdir <- "HS856"
rundir <- rundir
controlfile=controlfile
hsargs=hsargs
hcrfun=mcdahcr
sampleCE=tasCPUE
sampleFIS=tasFIS
sampleNaS=tasNaS
getdata=tasdata
calcpopC=calcexpectpopC
varyrs=7
startyr=38
verbose=TRUE
ndiagprojs=4
savesauout=TRUE
makehcrout=makeouthcr
cutcatchN=56
matureL = c(70,200)
wtatL = c(80,200)
mincount=100


zone <- makeequilzone(rundir=rundir,ctrlfile=controlfile,doproduct=TRUE,verbose=TRUE)

# declare main objects ----------------------------------------------------
glb <- zone$glb
ctrl <- zone$ctrl
zone1 <- zone$zone1
projC <- zone$zone1$projC
condC <- zone$zone1$condC
zoneC <- zone$zoneC
zoneD <- zone$zoneD
production <- zone$product
saudat <- zone$saudat
constants <- zone$constants
# save some equil results -------------------------------------------------
#Condition on Fishery -----------------------------------------------------
if (any(condC$initdepl < 1)) {
  initdepl <- condC$initdepl
  if (verbose) cat("Conducting initial depletions  ",initdepl,"\n")
  zoneD <- depleteSAU(zoneC,zoneD,glb,initdepl=initdepl,production)
}
if (verbose) cat("Conditioning on the Fishery data  \n")
zoneDD <- dohistoricC(zoneD,zoneC,glob=glb,condC,calcpopC=calcpopC,
                      sigR=1e-08,sigB=1e-08)
hyrs <- glb$hyrs
popdefs <- as.data.frame(t(getlistvar(zoneC,"popdef")))
popdefs[,c(1:6,8:18)] <- round(popdefs[,c(1:6,8:18)],3) #so SAU can be txt
popdefs[,"SAU"] <- glb$sauname[glb$sauindex]  #so SAU can be txt
# plot predicted size-comp of catch vs observed size-comps
catchN <- zoneDD$catchN
sauCt <- popNAStosau(catchN,glb)
compdat <- condC$compdat$lfs
if (!is.null(compdat)) {
  for (plotsau in 1:glb$nSAU) {
    lfs <- preparesizecomp(compdat[,,plotsau],mincount=mincount)
    yrsize <- as.numeric(colnames(lfs))
    histyr <- condC$histyr
    pickyr <- match(yrsize,histyr[,"year"])
    LML <- histyr[pickyr,"histLML"]
  }
}
# do projections ------------------------------------------------------------
if (verbose) cat("Doing the projections \n")
outpp <- prepareprojection(projC=projC,condC=condC,zoneC=zoneC,glb=glb,
                           calcpopC=calcpopC,zoneDD=zoneDD,
                           ctrl=ctrl,varyrs=varyrs,lastsigR = ctrl$withsigR)
zoneDP <- outpp$zoneDP
projC <- outpp$projC
zoneCP <- outpp$zoneCP
# naked doprojections -----------------------------------
reps <- ctrl$reps
sigmar <- ctrl$withsigR
sigmab <- ctrl$withsigB
sigce <- ctrl$withsigCE
hyrs <- glb$hyrs
startyr <- hyrs + 1
pyrs <- glb$pyrs
endyr <- hyrs + pyrs
sauindex <- glb$sauindex
saunames <- glb$saunames
yrnames <- c(glb$hyrnames,glb$pyrnames)
Nclass <- glb$Nclass
movem <- glb$move
r0 <- getvar(zoneCP,"R0") #R0 by population
b0 <- getvar(zoneCP,"B0") #sapply(zoneC,"[[","B0")
exb0 <- getvar(zoneCP,"ExB0")
hcrout <- makehcrout(glb,hsargs)
for (iter in 1:reps) {
  if (verbose) if ((iter %% 25) == 0) cat(iter,"   ")
  for (year in startyr:endyr) { # iter=1; year=startyr
    hcrdata <- getdata(sampleCE,sampleFIS,sampleNaS,
                       sauCPUE=zoneDP$cesau[,,iter],
                       sauacatch=zoneDP$acatch[,,iter],
                       sauNAS=list(Nt=zoneDP$Nt[,,,iter],
                                   catchN=zoneDP$catchN[,,,iter],
                                   NumNe=zoneDP$NumNe[,,,iter]),
                       year=year)
    hcrout <- hcrfun(hcrdata,hsargs,saunames=glb$saunames)
    popC <- calcpopC(hcrout,exb=zoneDP$exploitB[year-1,,iter],
                     sauindex,sigmab=sigmab)
    outy <- oneyearsauC(zoneCC=zoneCP,inN=zoneDP$Nt[,year-1,,iter],
                        popC=popC,year=year,Ncl=Nclass,sauindex=sauindex,
                        movem=movem,sigmar=sigmar,sigce,
                        r0=r0,b0=b0,exb0=exb0)
    dyn <- outy$dyn
    saudyn <- poptosauCE(dyn["catch",],dyn["cpue",],sauindex)
    zoneDP$exploitB[year,,iter] <- dyn["exploitb",]
    zoneDP$midyexpB[year,,iter] <- dyn["midyexpB",]
    zoneDP$matureB[year,,iter] <- dyn["matureb",]
    zoneDP$catch[year,,iter] <- dyn["catch",]
    zoneDP$acatch[year,,iter] <- hcrout$acatch
    zoneDP$catsau[year,,iter] <- saudyn$saucatch
    zoneDP$harvestR[year,,iter] <- dyn["catch",]/dyn["midyexpB",]
    zoneDP$cpue[year,,iter] <- dyn["cpue",]
    zoneDP$cesau[year,,iter] <- saudyn$saucpue
    zoneDP$recruit[year,,iter] <- dyn["recruits",]
    zoneDP$deplsB[year,,iter] <- dyn["deplsB",]
    zoneDP$depleB[year,,iter] <- dyn["depleB",]
    zoneDP$Nt[,year,,iter] <- outy$NaL
    zoneDP$catchN[,year,,iter] <- outy$catchN
    zoneDP$NumNe[,year,,iter] <- outy$NumNe
    zoneDP$TAC[year,iter] <- hcrout$TAC
  } # year loop
}   # iter loop
outproj <- list(zoneDP=zoneDP,hcrout=hcrout)

# outproj <- doprojections(ctrl=ctrl,zoneDP=zoneDP,zoneCP=zoneCP,glb=glb,
#                          hcrfun=hcrfun,hsargs=hsargs,sampleCE=sampleCE,
#                          sampleFIS=sampleFIS,sampleNaS=sampleNaS,
#                          getdata=getdata,calcpopC=calcpopC,
#                          makehcrout=makeouthcr,verbose=TRUE)
#  Rprof(NULL)
if (verbose) cat("Now generating final plots and tables \n")
zoneDP=outproj$zoneDP
hcrout <- outproj$hcrout; #str(hcrout)

NAS <- list(Nt=zoneDP$Nt,catchN=zoneDP$catchN)
# NumNe=zoneDP$NumNe, mid-year numbers-at-size removed to save space
zoneDP <- zoneDP[-c(17,16,15)]
# histCE <- condC$histCE
B0 <- getvar(zoneC,"B0")
ExB0 <- getvar(zoneC,"ExB0")
sauout <- sauplots(zoneDP,NAS,glb,rundir,B0,ExB0,
                   startyr=startyr,addCI=TRUE,histCE=condC$histCE)
outzone <- poptozone(zoneDP,NAS,glb,
                     B0=sum(getvar(zoneC,"B0")),
                     ExB0=sum(getvar(zoneC,"ExB0")))
NAS$catchN <- NAS$catchN[(cutcatchN:glb$Nclass),,,]

# calculate HS performance statistics
sum5 <- getprojyrC(catsau=zoneDP$catsau,glb=glb,period=5)
sum10 <- getprojyrC(catsau=zoneDP$catsau,glb=glb)
HSstats <- list(sum10=sum10,sum5=sum5)

out <- list(tottime=tottime,projtime=projtime,starttime=starttime,glb=glb,
            ctrl=ctrl,zoneCP=zoneCP,zoneD=zoneD,zoneDD=zoneDD,zoneDP=zoneDP,
            NAS=NAS,projC=projC,condC=condC,sauout=sauout,outzone=outzone,
            hcrout=hcrout,production=production,condout=condout,
            HSstats=HSstats,saudat=saudat,constants=constants,hsargs=hsargs)




getDBdir2 <- function() {
  if (dir.exists("./../../../Dropbox")) {
    prefixdir <- "./../../../Dropbox/"
  } else {
    prefixdir <- NULL
    stop("No dropbox directory found \n")
  }
  return(prefixdir)
}

x <- getDBdir2()
