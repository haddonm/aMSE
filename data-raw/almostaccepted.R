



library(aMSE)
library(codeutils)
library(hplot)


outdir <- "C:/aMSE_scenarios/EG/"

rundir=rundir;postfixdir=postfixdir;outdir=outdir;files=files;pickfiles=c(2,4,5,6)
 verbose=TRUE; intensity=100; zero=FALSE; altscenes=NULL
 juris="";ribbonleg="topleft"; Q90=TRUE; scencol=c(1,2,3,4)
# get files -------------------
files2 <- files[pickfiles]
nfile <- length(pickfiles)
label <- vector(mode="character",length=nfile)
for (i in 1:nfile) label[i] <- unlist(strsplit(files2[i],".",fixed=TRUE))[1]
ans <- makelist(label) # vector(mode="list",length=nfile)
dyn <- makelist(label)
glbc <- makelist(label)
ctrlc <- makelist(label)
condCc <- makelist(label)
prods <- makelist(label)
scenes <- vector(mode="character",length=nfile)
scores <- makelist(label)
zone <- makelist(label)
for (i in 1:nfile) { # i = 1
  if (nchar(outdir) == 0) {
    filename <- files2[i]
  } else {
    filename <- pathtopath(outdir,files2[i])
  }
  if (verbose)
    cat("Loading ",files2[i]," which may take time, be patient  \n")
  out <- NULL # so the function knows an 'out' exists
  load(filename)
  ans[[i]] <- out
  dyn[[i]] <- out$sauout
  glbc[[i]] <- out$glb
  ctrlc[[i]] <- out$ctrl
  condCc[[i]] <- out$condC
  prods[[i]] <- t(out$sauprod)
  scenes[i] <- out$ctrl$runlabel
  scores[[i]] <- out$outhcr
  zone[[i]] <- out$outzone
}





for (i in 1:8) print(getabdevs(catchS[1:58,i,1],years=yrs[1:58],n=7))

sau <- 8  # = 11
sauC <- catchS[59:88,sau,]
tmp <- apply(sauC,2,getabdevs,years=yrs[59:88],n=7)
tmp2 <- apply(sauC,2,getabdevs,years=yrs[59:88],n=7,smooth=FALSE)
plotprep(width=8, height=5)
parset(plots=c(2,1),margin=c(0.45,0.45,0.1,0.05))
hist(tmp,breaks=seq(0,25,0.5),main=round(mean(tmp),3),xlab="Abs Deviations")
hist(tmp2,breaks=seq(0,25,0.5),main=round(mean(tmp2),3),xlab="Abs Deviations")

# test how well the movav works on replicates




devout <- plotdevs(x,scenes,saunames,filen="")




glb <- glbc[[1]]
saunames <- glb$saunames
nsau <- glb$nSAU
nscen <- length(scenes)

x <- makelist(scenes)
for (scen in 1:length(scenes)) x[[scen]] <- dyn[[scen]]$catch[59:88,,]

devout <- plotdevs(x,scenes,saunames,filen="")



meandevs
sddevs


plotprep(width=8, height=9)
parset(plots=c(4,2))
for (i in 9:16) {
  yval <- sauC[59:88,i]
  years <- yrs[59:88]
  plot1(x=years,y=yval,lwd=2,defpar=FALSE)
  ma <- loess(yval ~ years)
 # ma <- movav(yval,n=7)
  lines(x=ma$x,y=ma$fitted,lwd=2,col=2)
  abline(v=2020.5,lwd=1,col=1)

}

sauC[,1:8]


tmp <- movav(rep1,n=7)

plot1(yrs[1:58],rep1[1:58])
lines(yrs[1:58],tmp[1:58],lwd=2,col=2)


pick <- which(!is.na(tmp[1:58]))

sum(abs(rep1[pick] - tmp[pick]))/length(pick)











# naked do_MSE without producing all the plots etc
#
postfixdir <- "EG"
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


options("show.signif.stars"=FALSE,
        "stringsAsFactors"=FALSE,
        "max.print"=50000,
        "width"=240)
suppressPackageStartupMessages({
  library(aMSE)
  library(TasHS)
  library(codeutils)
  library(hplot)
  library(makehtml)
  library(knitr)
})
dropdir <- getDBdir()
prefixdir <- paste0(dropdir,"A_codeUse/aMSEUse/hsargs/")

startime <- Sys.time()
postfixdir <- "EG"
verbose <- TRUE
rundir <- filenametopath(prefixdir,postfixdir)
controlfile <- paste0("control",postfixdir,".csv")
outdir <- "C:/aMSE_scenarios/hsargs/"
confirmdir(rundir)
confirmdir(outdir)


hsargs <- list(mult=0.1, # expansion factor for cpue range when calc the targqnt
               wid = 4, # number of years in the grad4 PM
               targqnt = 0.55, # quantile defining the cpue target
               maxtarg = c(150,150,150,150,150,150,150,150), # max cpue Target
               pmwts = c(0.65,0.25,0.1),  # relative weights of PMs
               hcr = c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2), # hcr multipliers
               startCE = 1992, # used in constant reference period HS
               endCE = 2020,   # used in constant reference period HS
               metRunder = 2,  # should the metarules be used. o =
               metRover = 2,   # use metarules
               decrement=1, # use fishery data up to the end of the time series
               pmwtSwitch = 0,
               stablewts = c(0.4, 0.5, 0.1),
               hcrname="mcdahcr")




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
zoneDD <- dohistoricC(zoneD,zoneC,glob=glb,condC,calcpopC=calcexpectpopC,
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
    lfs <- preparesizecomp(compdat[,,plotsau],mincount=120)
    yrsize <- as.numeric(colnames(lfs))
    histyr <- condC$histyr
    pickyr <- match(yrsize,histyr[,"year"])
    LML <- histyr[pickyr,"histLML"]
  }
}
# do projections ------------------------------------------------------------
if (verbose) cat("Doing the projections \n")
outpp <- prepareprojection(projC=projC,condC=condC,zoneC=zoneC,glb=glb,
                           calcpopC=calcexpectpopC,zoneDD=zoneDD,
                           ctrl=ctrl,varyrs=7,lastsigR = ctrl$withsigR)
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
hcrout <- makeouthcr(glb,hsargs)
impactyrs <- NULL
impactsurv <- NULL
if (!is.null(glb$envimpact)) {
  envimpact <- glb$envimpact
  impactyrs <- sapply(envimpact,"[[","eyr")
  impactsurv <- cbind(sapply(envimpact,"[[","proprec"),
                      sapply(envimpact,"[[","propNt"))
  colnames(impactsurv) <- c("proprec","propNt")
}
for (year in startyr:endyr) { # iter=1; year=startyr
  if (verbose) cat(year,"   ")
  for (iter in 1:reps) {
    hcrdata <- getdata(sampleCE,sampleFIS,sampleNaS,
                       sauCPUE=zoneDP$cesau[,,iter],
                       sauacatch=zoneDP$acatch[,,iter],
                       sauNAS=list(Nt=zoneDP$Nt[,,,iter],
                                   catchN=zoneDP$catchN[,,,iter],
                                   NumNe=zoneDP$NumNe[,,,iter]),
                       year=year,decrement=hsargs$decrement)
    hcrout <- hcrfun(hcrdata,hsargs,saunames=glb$saunames)
    popC <- calcpopC(hcrout,exb=zoneDP$exploitB[year-1,,iter],
                     sauindex,sigmab=sigmab,year=year)
    outy <- oneyearsauC(zoneCC=zoneCP,inN=zoneDP$Nt[,year-1,,iter],
                        popC=popC,year=year,Ncl=Nclass,sauindex=sauindex,
                        movem=movem,sigmar=sigmar,sigce,r0=r0,b0=b0,
                        exb0=exb0,envyr=impactyrs,envsurv=impactsurv)
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



plotconditioning(zoneDD=zoneDD,glb=glb,zoneC=zoneCP,
                 histCE=condC$histCE,histCatch=condC$histCatch,rundir=rundir,
                 recdevs=condC$recdevs,console=TRUE)


addccf <- function(condC,sau) {
  label <- "Correlation"
  cedat <- condC$histCE[,sau]
  cpuedat <- cedat[which(cedat > 0)]
  pickC <- which(condC$histyr %in% as.numeric(names(cpuedat)))
  catdat <- condC$histCatch[pickC,sau]
  ccfout <- ccf(x=catdat,y=cpuedat,type="correlation",
                ylab=label,plot=TRUE,xlab="Lag Years")
} # end of plotccf


addccf(condC,sau=1)




