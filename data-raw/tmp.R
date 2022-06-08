



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
postfixdir <- "M15h7L75_lowvar"
verbose <- TRUE
#hsfile <- "TasHS1_Tas.R"
rundir <- filenametopath(prefixdir,postfixdir)
controlfile <- paste0("control",postfixdir,".csv")
outdir <- "C:/aMSE_scenarios/M15h7L75/"
confirmdir(outdir)

#ctrlfiletemplate(rundir,filename="testcontrolfile.csv")
#source(paste0(dropdir,"A_Code/aMSE/data-raw/almostaccepted.R"))

#source(paste0(rundir,"/",hsfile))
hsargs <- list(mult=0.1,
               wid = 4,
               targqnt = 0.55,
               maxtarg = c(150,150,150,150,150,150,150,150),
               pmwts = c(0.65,0.25,0.1),
               hcr = c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2),
               startCE = 1992)
doproduct=FALSE

#zone <- makeequilzone(rundir,controlfile,verbose=verbose)

zone1 <- readctrlfile(rundir,infile=ctrlfile,verbose=verbose)
ctrl <- zone1$ctrl
glb <- zone1$globals     # glb without the movement matrix
bysau <- ctrl$bysau
if (is.null(bysau)) bysau <- 0
if (bysau) {
  saudata <- readsaudatafile(rundir,ctrl$datafile)
  constants <- saudata$constants
  saudat <- saudata$saudat
} else {
  constants <- readdatafile(glb$numpop,rundir,ctrl$datafile)
  saudat <- constants
}
if (verbose) cat("Files read, now making zone \n")
#out <- setupzone(constants,zone1,doproduct,verbose=verbose) # make operating model
constants=constants; zone1=zone1; doproduct=FALSE; uplim=0.4; inc=0.001; verbose=TRUE
ans <- makezoneC(zone1,constants) # initiates zoneC
zoneC <- ans$zoneC
glb <- ans$glb
ans <- makezone(glob=glb,zoneC=zoneC) # make zoneD, add cpue, qest to zoneC
zoneC <- ans$zoneC  # zone constants
zoneD <- ans$zoneD  # zone dynamics
product <- NULL
if (doproduct) {
  if (verbose) cat("Now estimating population productivity \n")
  # adds productivity, and MSY, MSYdepl to zoneC if doproduct=TRUE
  ans <- modzoneC(zoneC=zoneC,zoneD=zoneD,glob=glb,uplim=uplim,inc=inc)
  zoneC <- ans$zoneC  # zone constants
  product <- ans$product  # productivity by population
}
out <- list(zoneC=zoneC, zoneD=zoneD, product=product,glb=glb)

zoneC <- out$zoneC
zoneD <- out$zoneD
glb <- out$glb             # glb now has the movement matrix
product <- out$product     # important bits usually saved in rundir
zone1$globals <- glb
zone <- list(zoneC=zoneC,zoneD=zoneD,glb=glb,constants=constants,
             saudat=saudat,product=product,ctrl=ctrl,zone1=zone1)
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

str1(zoneD)

pop=1
SurvE <- exp(-zoneC[[pop]]$Me)
hSurv <- exp(-zoneC[[pop]]$Me/2.0)
recr <- rep(0,N)
recr[1] <- newrecs[pop]  # change this once I have imported move
UnitM <- matrix(0,nrow=N,ncol=N)
diag(UnitM) <- 1.0
G <- zoneC[[pop]]$G
Minv <- solve(UnitM - (SurvE * G))
Nt[,1,pop] <- Minv %*% recr # initial unfished numbers-at-size
MatB[1,pop] <- sum(zoneC[[pop]]$MatWt*Nt[,1,pop])/1e06
zoneC[[pop]]$B0 <- MatB[1,pop] # mature biomass at start of year
newNt1 <- (hSurv * (G %*% Nt[,1,pop]))
newNt2 <- (hSurv * newNt1)
popSel <- zoneC[[pop]]$SelWt[,1]
preexB <- sum(popSel*newNt1)/1e06
postexB <- sum(popSel*newNt2)/1e06
ExplB[1,pop] <- (preexB + postexB)/2 # mean exploitB before/after halfyear
zoneC[[pop]]$ExB0 <- ExplB[1,pop]
deplExB[1,pop] <- 1.0  # no depletion when first generating zones
deplSpB[1,pop] <- 1.0
Recruit[1,pop] <- recr[1]
#qcalc <- as.numeric(zoneC[[pop]]$popdef["MaxCE"])
zoneC[[pop]]$scalece<- 0 # used to scale all pops to maxce



dyn <- out$condout$sauZone

rownames(dyn$expB) <- glb$hyrnames

head(round(dyn$expB,3),20)

head(round(dyn$harvestR,3),20)


x <- 1:200

sm <- maturity(-16,0.154,x)

am <- maturity(-22.2536,0.19119,x)


sm <- WtatLen(5.62e-05,3.161963,x)
am <- WtatLen(6.868252e-05,3.147,x)


plotprep(width=8,height=4.5,newdev=FALSE)
parset()
plot1(x,sm,lwd=2)
lines(x,am,lwd=2,col=2)


round(cbind(sm,am),4)







# qest,3.15,1.50,3.3,0.975,0.78,0.44,0.445,0.995
# "LnR0",327321,475662,298599,1188198,975141,1977949,1682092,609722
# "MaxDL",19.6196945721412,26.5190213803115,24.3452875729363,29.1005068921024,29.2957741405674,28.6797788699979,26.7580375179069,21.3828278617031
# "L95",183.789931015386,172.739696769617,183.450201400732,179.267909313624,183.701518577507,182.377795633067,178.929810469837,184.034627956569
# "qest",4.35479652151046,2.20203182483022,5.10468337057917,1.63447806684747,1.27888832322373,0.799525848626615,0.613938218148408,1.43244227081343
# "seldelta",5.48519006329542,5.54217168224989,3.86051291751843,3.26973770119845,2.86997981991041,4.08887960529036,4.62259595059192,3.90026307554492
# L50C ,126.4222,126.4222,126.4222,126.4222,126.4222,126.4222,126.4222,126.4222, length at 50% emergent
# sL50C ,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,




getparam <- function(rundir,datafile,nsau,varname) {
  filen <- filenametopath(rundir,datafile)
  dat <- readLines(filen)
  pickA <- grep(varname,dat)[1] # ignore sAvRec
  param <- getConst(dat[pickA],nsau)
  return(param)
} # end of getavrec


saussq <- function(param,rundir,controlfile,datafile,varname,
                   calcpopC,fillq,picksau,nsau,outplot=FALSE) {
# param=2.381915; rundir=rundir;controlfile=controlfile;datafile=datafile;varname="qest"
# calcpopC=calcexpectpopC;fillq=fillq;picksau=sau;nsau=8;outplot=FALSE
  nsaum1 <- nsau - 1
  if (picksau == 1)
    replacetxt <- paste0(varname,",",as.character(param[1]),",",
                         paste0(as.character(fillq[1:nsaum1]),collapse=","),
                         collapse=",")
  if ((picksau > 1) & (picksau < nsau))
    replacetxt <- paste0(varname,",",
                         paste0(as.character(fillq[1:(picksau-1)]),collapse=","),
                         ",",as.character(param[1]),",",
                         paste0(as.character(fillq[picksau:nsaum1]),collapse=","))
  if (picksau == nsau)
    replacetxt <- paste0(varname,",",
                         paste0(as.character(fillq[1:nsaum1]),collapse=","),
                         ",",as.character(param[1]),collapse=",")
  changeline(rundir,datafile,varname,replacetxt)
  zone <- makeequilzone(rundir,controlfile,doproduct=FALSE,
                        verbose=FALSE)
  glb <- zone$glb
  condC <- zone$zone1$condC
  zoneC <- zone$zoneC
  zoneD <- zone$zoneD
  zoneDD <- dohistoricC(zoneD,zoneC,glob=glb,condC,calcpopC=calcpopC,
                        sigR=1e-08,sigB=1e-08)
  hyrs <- glb$hyrs
  sauindex <- glb$sauindex
  popB0 <- getlistvar(zoneC,"B0")
  B0 <- tapply(popB0,sauindex,sum)
  popExB0 <- getlistvar(zoneC,"ExB0")
  ExB0 <- tapply(popExB0,sauindex,sum)
  sauZone <- getsauzone(zoneDD,glb,B0=B0,ExB0=ExB0)
  if (outplot) {
    filename=paste0("compareCPUE_",picksau,"_AvRec_cond.png")
    ssq <- plotcondCPUE(condC$histCE,sauZone$cpue,glb,rundir,filen=filename)
  } else {
    ssq <- getCPUEssq(condC$histCE,sauZone$cpue,glb)
  }
  return(ssq[picksau])
} # end of sauavrecssq


startime <- Sys.time()
postfixdir <- "M15h7L75_lowvar"
verbose <- TRUE
#hsfile <- "TasHS1_Tas.R"
rundir <- filenametopath(prefixdir,postfixdir)
controlfile <- paste0("control",postfixdir,".csv")
#source(paste0(datadir,"/",hsfile))
#source(paste0(DBdir,"A_Code/aMSE/data-raw/almostaccepted.R"))
datafile <- "saudataM15h7L75_lowvar.csv"
nsau <- 8
final <- numeric(nsau)


for (sau in 1:2) { #  sau=1
  qest <- getparam(rundir,datafile,nsau=nsau,varname="qest")
  param <- qest[sau]
  fillq <- qest[-sau]
  low <- param * 0.6
  high <- param * 1.4
  cat("Running for sau ",sau,"\n")
  origssq <- saussq(param,rundir,controlfile,
                      datafile=datafile,varname="qest",
                      calcpopC=calcexpectpopC,fillq=fillq,
                      picksau=sau,nsau=nsau)
  cat("oldssq ",origssq,"   orig param ",param,"\n")
  ans <- optim(param,saussq,method="Brent",lower=low,upper=high,
               rundir=rundir,
               controlfile=controlfile,datafile=datafile,varname="qest",
               calcpopC=calcexpectpopC,fillq=fillq,picksau=sau,nsau=nsau,
               control=list(maxit=30))
  cat("ssq = ",ans$value,"  qest value = ",ans$par,"\n")
  if (((ans$par - low) < 0) | ((high - ans$par) < 0))
    warning(cat("Boundary reached for param ",sau,low,high,ans$par,"\n"))
  final[sau] <- round(ans$par,5)

  cat(sau,"   ",low,"    ",param,"   ",round(ans$par,4),"   ",high,"\n")
  cat("old ",round(origssq,1),"     new ",round(ans$value,1),"\n\n")
  cat("old ",round(param,1),"     new ",round(ans$par,1),"\n\n")
}
endtime <- Sys.time()
print(initial); print(final)
print(endtime - startime)



out <- do_condition(rundir,controlfile,calcpopC=calcexpectpopC,
                    verbose = TRUE,
                    doproduct = FALSE)

makeoutput(out,rundir,postfixdir,controlfile,hsfile="TasHS Package",
           openfile=TRUE,verbose=FALSE)



# test SA recdevs----------------------------

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
postfixdir <- "SA"
verbose <- TRUE
#hsfile <- "TasHS1_Tas.R"
rundir <- filenametopath(prefixdir,postfixdir)
controlfile <- paste0("control",postfixdir,".csv")
outdir <- "C:/aMSE_scenarios/M15h7L75/"
confirmdir(outdir)

#ctrlfiletemplate(rundir,filename="testcontrolfile.csv")
#source(paste0(dropdir,"A_Code/aMSE/data-raw/almostaccepted.R"))

#source(paste0(rundir,"/",hsfile))
hsargs <- list(mult=0.1,
               wid = 4,
               targqnt = 0.55,
               maxtarg = c(150,150,150,150,150,150,150,150),
               pmwts = c(0.65,0.25,0.1),
               hcr = c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2),
               startCE = 1992)
doproduct=FALSE
wtsc=c(7e-07,2e-06,4e-05,7e-07,5e-07,3e-07,1e-07,1e-07)

zone <- makeequilzone(rundir,controlfile,doproduct=FALSE,
                      verbose=FALSE)
# declare main objects
glb <- zone$glb
condC <- zone$zone1$condC
zoneC <- zone$zoneC
zoneD <- zone$zoneD
condC <- zone$zone1$condC
zoneDD <- dohistoricC(zoneD,zoneC,glob=glb,condC,calcpopC=calcpopC,
                      sigR=1e-08,sigB=1e-08)
hyrs <- glb$hyrs
sauindex <- glb$sauindex
popB0 <- getlistvar(zoneC,"B0")
B0 <- tapply(popB0,sauindex,sum)
popExB0 <- getlistvar(zoneC,"ExB0")
ExB0 <- tapply(popExB0,sauindex,sum)
sauZone <- getsauzone(zoneDD,glb,B0=B0,ExB0=ExB0)
ssq <- getCPUEssq(condC$histCE,sauZone$cpue,glb)
names(ssq) <- glb$saunames
nsau <- glb$nSAU
LFlog <- numeric(nsau); names(LFlog) <- glb$saunames
minsccount <- mincount
scwt <- wtsc
for (sau in 1:nsau) {
  if (length(mincount) > 1) minsccount <- mincount[sau]
  if (length(wtsc) > 1) scwt <- wtsc[sau]
  obsLFs <- preparesizecomp(condC$compdat$lfs[,,sau],mincount=minsccount)
  LFlog[sau] <- getLFlogL(zoneDD$catchN,obsLFs,out$glb,wtsc=scwt,sau=sau)
}
totssq <- ssq + LFlog
ans <- rbind(totssq,ssq,LFlog)


