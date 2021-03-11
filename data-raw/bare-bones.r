










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
# Obviously you should modify the rundir to suit your own computer
if (dir.exists("c:/Users/User/DropBox")) {
  ddir <- "c:/Users/User/DropBox/A_code/"
} else {
  ddir <- "c:/Users/Malcolm/DropBox/A_code/"
}
#rundir <- paste0(ddir,"aMSEUse/conddata/generic2")
rundir <- paste0(ddir,"aMSEUse/conddata/generic")
dirExists(rundir,make=TRUE,verbose=TRUE)
# this should be in rundir, but during development is in data-raw
source(paste0(ddir,"aMSE/data-raw/","TasmanianHS.R"))
# source(paste0(rundir,"/TasmanianHS.R"))

zone <- makeequilzone(rundir,"control2.csv") # normally would read in a file
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
# this should be in rundir, but during development is in data-raw


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
# resfile <- setuphtml(rundir)# prepare to save and log results
#
# plotproductivity(rundir,product,glb)
# biology_plots(rundir, glb, zoneC)
# numbersatsize(rundir, glb, zoneD)
#
# endtime <- as.character(Sys.time())
#
# reportlist <- list(starttime=starttime,endtime=endtime,
#                    zoneC=zoneC, zoneD=zoneD, product=product,
#                    glb=glb,constants=constants)
# runnotes <- "This is a bare-bones example."
# # If you unhash this component it will generate a local website inside
# # rundir and open it so you can see the results so far.
# make_html(replist=reportlist,rundir=rundir,width=500,
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
#   optCE <- MCDAdata(rundir,mcdafile,zone1$SAUnames)
# }
#
#  zoneCP=zoneCP;zoneDP=zoneDR;glob=glb;ctrl=ctrl;projyrs=projC$projyrs;inityrs=projC$inityrs;



condC <- zone$zone1$condC

histce <- condC$histCE

cmcda <- calibrateMCDA(histce,saunames)







