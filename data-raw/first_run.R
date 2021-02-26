
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
starttime <- (Sys.time())
zone <- makeequilzone(resdir,"control2.csv") # normally would read in a file
equiltime <- (Sys.time()); print(equiltime - starttime)

glb <- zone$glb
ctrl <- zone$ctrl
zone1 <- zone$zone1
projC <- zone1$projC
condC <- zone1$condC
zoneC <- zone$zoneC
zoneD <- zone$zoneD
product <- zone$product
# save objects to resdir
save(ctrl,file=filenametopath(resdir,"ctrl.RData"))
save(glb,file=filenametopath(resdir,"glb.RData"))
save(product,file=filenametopath(resdir,"product.RData"))
save(zoneC,file=filenametopath(resdir,"zoneC.RData"))
save(condC,file=filenametopath(resdir,"condC.RData"))

# condition on the historic catches
zoneDD <- dohistoricC(zoneD,zoneC,glob=glb,condC,sigR=1e-08,sigB=1e-08)
midtime <- (Sys.time())
print(equiltime - starttime)
save(zoneDD,file=filenametopath(resdir,"zoneDD.RData"))
# Illustrate productivity
propD <- getzoneprops(zoneC,zoneDD,glb,year=47)
addtable(round(propD,4),"propertyDD.csv",resdir,category="zoneDD",caption=
           "Properties of zoneD after conditioning on historical catches.")
addtable(round(t(zoneDD$harvestR[45:47,]),4),"final_harvestR.csv",resdir,
         category="zoneDD",caption="Last three years of harvest rate.")
popdefs <- getlistvar(zone$zoneC,"popdef")
addtable(round(t(popdefs),3),"popdefs.csv",resdir,category="zoneDD",caption=
           "Population specific definitions")

# Do the replicates ------------------------------------------------------------
midtime <- (Sys.time()); print(midtime - starttime)

cmcda <- mcdahcr(arrce=condC$histCE,hsargs=hsargs,
                 yearnames=rownames(condC$histCE),saunames=glb$saunames)
#  str(cmcda)
pms <- cmcda$pms
multTAC <- cmcda$multTAC

out <- prepareprojection(projC,zoneC,glb,zoneDD,ctrl,multTAC)

zoneDP <- out$zoneDP
projC <- out$projC
zoneCP <- out$zoneCP
save(zoneCP,file=filenametopath(resdir,"zoneCP.RData"))
save(projC,file=filenametopath(resdir,"projC.RData"))

endtime <- Sys.time(); print(endtime - midtime)

# doprojections
# str(cmcda)
begintime <- Sys.time()

zoneDP <- doTASprojections(ctrl,zoneDP,zoneCP,condC$histCE,glb,mcdahcr,hsargs)

projtime <- Sys.time()
print(projtime - begintime); print(projtime - starttime)
save(zoneDP,file=filenametopath(resdir,"zoneDP.RData"))


result <- alltosau(zoneDP,glb)


plotprep(width=8, height=8,newdev=FALSE)
plotsau(invar=result$cesau,glb=glb,plots=c(4,2),ylab="CPUE",medcol=1,addCI=TRUE)


plotprep(width=8, height=8,newdev=FALSE)
plotsau(invar=zoneDP$catsau,glb=glb,plots=c(4,2),ylab="sau catch",medcol=1,addCI=TRUE)






plotprep(width=7, height=4,newdev=FALSE)
parset()
hist(zoneDP$catsau[1,1,])



















str(zoneDP,max.level = 1)



plotprep(width=8, height=8,newdev=FALSE)
parset(plots=c(4,2),byrow=FALSE)
label <- glb$saunames
arrce <- zoneDP$cesau[,,1]
rownames(arrce) <- 2020:2049
yrs <- 2020:2049
reps <- 100
pm <- mcdahcr(arrce,hsargs,yearnames=yrs,saunames=glb$saunames)
first <- pm$multTAC
for (sau in 1:8) {
  plot(yrs,first[,sau],type="l",lwd=1,col="grey",xlab="",ylab=label[sau],
       ylim=c(0,1.2),panel.first=grid())
  for (i in 2:reps) {
    arrce <- zoneDP$cesau[,,i]
    pm <- mcdahcr(arrce,hsargs,yearnames=yrs,saunames=glb$saunames)
    lines(yrs,pm$multTAC[,sau])
  }
}

invar=zoneDP$catch;

lt5 <- function(invect) {  invect <- sumc
  n <- length(invect)
  divvect <- abs(1- invect[1:(n-1)]/invect[2:n])
  pick <- which(divvect < 0.05)
  return(pick)
}


fivep <- function(invar,glb) { #  invar <- zoneDP$catsau; glb=glb
  numpop <- glb$numpop
  nsau <- glb$nSAU
  nyrs <- dim(invar)[1]
  reps <- dim(invar)[3]
  meta5 <- matrix(0,nrow=nyrs,ncol=reps,dimnames=list(1:nyrs,1:reps))
  TAC <- meta5
  for (iter in 1:reps) {  # iter=1
    sumc <- rowSums(invar[,,iter])
    TAC[,iter] <- sumc
    meta5[lt5(sumc),iter] <- 1
  }
  return(list(meta5=meta5,TAC=TAC))
}

out <- fivep(zoneDP$catch,glb)










