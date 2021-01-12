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
#data(zone)
zone <- makeequilzone(resdir,"control2.csv") # normally would read in a file
  equiltime <- (Sys.time())
zoneDD <- dohistoricC(zone$zoneD,zone$zoneC,glob=zone$glb,zone$zone1$condC)
  midtime <- (Sys.time())

  print(equiltime - starttime)
  print(midtime - equiltime)

  propD <- getzoneprops(zone$zoneC,zoneDD,zone$glb,year=47)
  round(propD,3)

# Do the replicates ------------------------------------------------------------

midtime <- (Sys.time())

glb <- zone$glb
zone1 <- zone$zone1
projC <- zone1$projC
condC <- zone1$condC
# histCE <- zone$zone1$condC$histCE; saunames=zone$zone1$SAUnames

# ans <- calibrateMCDA(zoneDD$catch,zoneDD$cpue,projC$inityrs:glb$Nyrs,
#                      nsau=glb$nSAU,sauindex=glb$sauindex,pyrs=projC$projyrs)

cmcda <- calibrateMCDA(histCE=condC$histCE, saunames=zone1$SAUnames)

str(ans)

saucpue <- ans$saucpue
saucatch <- ans$saucatch
lastyr <- glb$Nyrs
histyr <- projC$inityrs:lastyr



ctrl <- zone$ctrl
zoneC=zone$zoneC

hsargs <- list(wid = 4,targqnt = 0.55, pmwts = c(0.65, 0.25,0.1),
               hcr = c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2),only=TRUE)
applyHS=mcdahcr;
HSargs=hsargs;


















saunames <- zone$zone1$SAUnames
sauindex <- glb$sauindex
pyrs <- projC$projyrs
B0 <- tapply(sapply(zone$zoneC,"[[","B0"),sauindex,sum)
exB0 <- tapply(sapply(zone$zoneC,"[[","ExB0"),sauindex,sum)



# use depleteSAU --------------------------------------------------------
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
#data(zone)
zone <- makeequilzone(resdir,"control2.csv") # normally would read in a file
equiltime <- (Sys.time())


glb <- zone$glb
zone1 <- zone$zone1
projC <- zone1$projC
condC <- zone1$condC
HSargs <- list(wid = 4,targqnt = 0.55, pmwts = c(0.65, 0.25,0.1),
               hcr = c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2),only=TRUE)

# histCE <- zone$zone1$condC$histCE; saunames=zone$zone1$SAUnames

# ans <- calibrateMCDA(zoneDD$catch,zoneDD$cpue,projC$inityrs:glb$Nyrs,
#                      nsau=glb$nSAU,sauindex=glb$sauindex,pyrs=projC$projyrs)

cmcda <- calibrateMCDA(histCE=condC$histCE, saunames=zone1$SAUnames,
                       hcrargs=HSargs)
str(cmcda)




#Rprof()
source(paste0(resdir,"/alternative_HS.R"))

HSargs <- list(wid = 4,targqnt = 0.55, pmwts = c(0.65, 0.25,0.1),
               hcr = c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2),only=TRUE)

mseproj <- doprojection(zoneC,zoneDP,glb,ctrl,projC$projyrs,applyHS=mcdahcr,
                        HSargs)
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







