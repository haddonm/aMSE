
# Used in readme.Rmd to illustrate the running of aMSE
# BEGIN ------------------------------------------------------------------------
options("show.signif.stars"=FALSE,
        "stringsAsFactors"=FALSE,
        "max.print"=50000,
        "width"=240)
# declare libraries ------------------------------------------------------------
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
resdir <- paste0(ddir,"aMSEUse/conddata/generic")
dirExists(resdir,make=TRUE,verbose=TRUE)
# this should be in resdir, but during development is in data-raw
source(paste0(ddir,"aMSE/data-raw/","TasmanianHS.R"))
# generate equilibrium zone ----------------------------------------------------
starttime2 <- (Sys.time())
zone <- makeequilzone(resdir,"control2.csv",cleanslate = TRUE) # normally would read in a file
equiltime <- (Sys.time()); print(equiltime - starttime2)
# declare main objects ---------------------------------------------------------
glb <- zone$glb
ctrl <- zone$ctrl
zone1 <- zone$zone1
projC <- zone1$projC
condC <- zone1$condC
zoneC <- zone$zoneC
zoneD <- zone$zoneD
product <- zone$product
# save equil results -----------------------------------------------------------
# save objects to resdir
# save(ctrl,file=filenametopath(resdir,"ctrl.RData"))
# save(glb,file=filenametopath(resdir,"glb.RData"))
# save(product,file=filenametopath(resdir,"product.RData"))
# save(zoneC,file=filenametopath(resdir,"zoneC.RData"))
# save(condC,file=filenametopath(resdir,"condC.RData"))
biology_plots(resdir, glb, zoneC)
plotproductivity(resdir,product,glb)
numbersatsize(resdir, glb, zoneD)
# Condition on Fishery ---------------------------------------------------------
zoneDD <- dohistoricC(zoneD,zoneC,glob=glb,condC,sigR=1e-08,sigB=1e-08)
# save conditioned results -----------------------------------------------------
save(zoneDD,file=filenametopath(resdir,"zoneDD.RData"))
# Illustrate productivity
propD <- getzoneprops(zoneC,zoneDD,glb,year=47)
addtable(round(t(propD),4),"propertyDD.csv",resdir,category="zoneDD",caption=
           "Properties of zoneD after conditioning on historical catches.")
addtable(round(t(zoneDD$harvestR[45:47,]),4),"final_harvestR.csv",resdir,
         category="zoneDD",caption="Last three years of harvest rate.")
popdefs <- getlistvar(zone$zoneC,"popdef")
addtable(round(t(popdefs),3),"popdefs.csv",resdir,category="zoneDD",caption=
           "Population specific definitions")
# Prepare projections ----------------------------------------------------------
cmcda <- mcdahcr(arrce=condC$histCE,hsargs=hsargs,
                 yearnames=rownames(condC$histCE),saunames=glb$saunames)
pms <- cmcda$pms
multTAC <- cmcda$multTAC
out <- prepareprojection(projC,zoneC,glb,zoneDD,ctrl,multTAC)
# out <- prepareprojectionnew(projC=projC,condC=condC,zoneC=zoneC,glb=glb,
#                             zoneDep=zoneDD,ctrl=ctrl,varyrs=6,lastsigR = 0.075)

zoneDP <- out$zoneDP
projC <- out$projC
zoneCP <- out$zoneCP
#zoneDDR <- out$zoneDDR

save(zoneCP,file=filenametopath(resdir,"zoneCP.RData"))
save(projC,file=filenametopath(resdir,"projC.RData"))

# do projections ---------------------------------------------------------------

zoneDP <- doTASprojections(ctrl,zoneDP,zoneCP,condC$histCE,glb,mcdahcr,hsargs)

#save(zoneDP,file=filenametopath(resdir,"zoneDP.RData"))

# save results to resdir -------------------------------------------------------
out <- plotbysau(zoneDP,glb)

projtime <- Sys.time()
print(projtime - starttime)

# make results webpage ---------------------------------------------------------
replist <- list(starttime=as.character(starttime),endtime=as.character(projtime))
make_html(
  replist = replist,
  resdir = resdir,
  width = 500,
  openfile = TRUE,
  runnotes = NULL,
  verbose = FALSE,
  packagename = "aMSE",
  htmlname = "aMSE"
)

# END --------------------------------------------------------------------------


cats <- poptosau(zoneDDR[["catch"]],glb)
exB <- poptosau(zoneDDR[["exploitB"]],glb)


str(zoneDP,max.level = 1)


str(zoneDDR,max.level = 1)

invar <- "cpue"
startyr <- 35
prerep <- poptosau(zoneDDR[[invar]],glb)
plotprep(width=8, height=8,newdev=FALSE)
parset(plots=c(4,2))
for (sau in 1:nsau) {
  ymax <- getmax(prerep[startyr:preyrs,sau,])
  plot(startyr:47,prerep[startyr:preyrs,sau,1],type="l",lwd=1,col="grey",
       panel.first=grid(),ylab=paste0(invar,"  ",label[sau]),ylim=c(0,ymax))
  for (i in 2:reps)
    lines(startyr:47,prerep[startyr:preyrs,sau,i],lwd=1,col="grey")
  if (invar == "cpue") abline(h=150,col=2)
}


plot(startyr:47,zoneDDR$harvestR[startyr:47,,1],type="l")

cats[startyr:47,,1]/exB[startyr:47,,1]


allrep <- poptosau(zoneDP[[invar]],glb)
preyrs <- dim(prerep)[1]
postyrs <- ctrl$projection
allyrs <- preyrs + postyrs
nsau <- glb$nSAU
reps <- dim(prerep)[3]
label <- glb$saunames

plotprep(width=8, height=8,newdev=FALSE)
parset(plots=c(4,2))
for (sau in 1:nsau) {
  ymax <- getmax(rbind(prerep[startyr:preyrs,sau,],allrep[,sau,]))
  traj <- c(prerep[startyr:preyrs,sau,1],allrep[,sau,1])
  plot(startyr:allyrs,traj,type="l",lwd=1,col="grey",panel.first=grid(),
       ylab=paste0(invar,"  ",label[sau]),ylim=c(0,ymax))
  for (i in 2:reps)
    lines(startyr:allyrs,c(prerep[startyr:preyrs,sau,i],allrep[,sau,i]),lwd=1,col="grey")
  abline(v=preyrs,col=2)
 if (invar == "cpue") abline(h=150,col=2)
}


range(allrep[,3,])



plotprep(width=8, height=8,newdev=FALSE)
plotsau(zoneDP$cpue,glb,c(4,2),ylab="CPUE")


zoneDD$cpue[,2]







str(zoneDP,max.level = 1)

# l ----------------------------------------------------------------------------
# plot total score -------------------------------------------------------------
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
  abline(h=1.0,col=2,lwd=2)
}

# plot cesau -------------------------------------------------------------------
plotprep(width=8, height=8,newdev=FALSE)
parset(plots=c(4,2),byrow=FALSE)
label <- glb$saunames
yrs <- 2020:2049
reps <- 100
for (sau in 1:8) {
  ymax <- getmax(zoneDP$cesau[,sau,])
  plot(yrs,zoneDP$cesau[,sau,1],type="l",lwd=1,col="grey",xlab="",ylab=label[sau],
       ylim=c(0,ymax),panel.first=grid())
  for (i in 2:reps) {
    arrce <-
    lines(yrs,zoneDP$cesau[,sau,i],lwd=1,col="grey")
  }
  abline(h=150,lwd=2,col=2)
}











invar=zoneDP$catch;
#
# lt5 <- function(invect) { # invect <- sumc
#   n <- length(invect)
#   divvect <- abs(1- invect[1:(n-1)]/invect[2:n])
#   pick <- which(divvect < 0.05)
#   return(pick)
# }
#
#
# fivep <- function(invar,glb) { #  invar <- zoneDP$catsau; glb=glb
#   numpop <- glb$numpop
#   nsau <- glb$nSAU
#   nyrs <- dim(invar)[1]
#   reps <- dim(invar)[3]
#   meta5 <- matrix(0,nrow=nyrs,ncol=reps,dimnames=list(1:nyrs,1:reps))
#   TAC <- meta5
#   for (iter in 1:reps) {  # iter=1
#     sumc <- rowSums(invar[,,iter])
#     TAC[,iter] <- sumc
#     meta5[lt5(sumc),iter] <- 1
#   }
#   return(list(meta5=meta5,TAC=TAC))
# }
#
# out <- fivep(zoneDP$catch,glb)










