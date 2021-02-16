





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
library(MQMF)
# Obviously you should modify the resdir to suit your own computer
if (dir.exists("c:/Users/User/DropBox")) {
  ddir <- "c:/Users/User/DropBox/A_code/"
} else {
  ddir <- "c:/Users/Malcolm/DropBox/A_code/"
}
#resdir <- paste0(ddir,"aMSEUse/conddata/generic2")
resdir <- paste0(ddir,"aMSEUse/conddata/generic")
dirExists(resdir,make=TRUE,verbose=TRUE)


zone1 <- readctrlfile(resdir,infile="control2.csv")

histC <- zone1$condC$histCatch
histyr <- zone1$condC$histyr
histCE <- zone1$condC$histCE


sau <- 8

catch <- histC[,sau]
cpue <- histCE[,sau]
nce <- length(cpue)
nca <- length(catch)
yrs <- as.numeric(names(catch))
pick <- match(as.numeric(names(cpue)),yrs)
nce <- length(pick)
if (nce < nca) {
  cpue <- c(rep(NA,(nca-nce)),cpue)
}

saudat <- cbind(year=yrs,catch=catch,cpue=cpue)
saudat
saudat <- as.matrix(saudat[20:47,])
yrs <- as.numeric(rownames(saudat))
plotprep(width=7, height=8)
parset(plots=c(3,1))
ymax <- getmax(saudat[,"catch"])
plot(yrs,saudat[,"catch"],type="l",lwd=2,ylim=c(0,ymax),yaxs="i",ylab="catches",
     panel.first = grid())
pick <- which(saudat[,"cpue"] > 0)
ccf(x=saudat[pick,"catch"],y=saudat[pick,"cpue"],type="correlation",
    ylab="Correlation",plot=TRUE)
ymax2 <- getmax(saudat[,"cpue"])
plot(yrs,saudat[,"cpue"],type="l",lwd=2,ylim=c(0,ymax2),yaxs="i",ylab="cpue",
     panel.first=grid())




# sau 7 param <- log(c(r=0.3,K=1000,Binit=500,sigma=0.5))  # SAU7
# sau10 param <- log(c(r=0.5,K=4000,Binit=1500,sigma=0.5)) # SAU10
# sau11 param <- log(c(r=0.35,K=5000,Binit=2000,sigma=0.5)) # SAU11
# SAU12 param <- log(c(r=0.35,K=5000,Binit=2000,sigma=0.5)) # SAU12
param <- log(c(r=0.35,K=3000,Binit=1000,sigma=0.5)) # SAU13
best <- optim(param,fn=negLL,funk=simpspm,indat=saudat,logobs=log(saudat[,"cpue"]),
              method="BFGS")
outfit(best,digits=4,title="SAU 9")
cat("\n")
best2 <- nlm(negLL,best$par,funk=simpspm,indat=saudat,logobs=log(saudat[,"cpue"]),
             steptol=1e-5)
outfit(best2,digits=4,title="SAU 9")
best3 <- nlm(negLL,best2$estimate,funk=simpspm,indat=saudat,logobs=log(saudat[,"cpue"]),
             steptol=1e-5)
outfit(best3,digits=4,title="SAU 9")
ans <- spm(best2$estimate,saudat)
ans$msy

plotprep(width=7, height=8)
parset(plots=c(3,1))
ymax <- getmax(saudat[,"catch"])
plot(yrs,saudat[,"catch"],type="l",lwd=2,ylim=c(0,ymax),yaxs="i",ylab="catches",
     panel.first = grid())
abline(h=ans$msy,lwd=2,col=2)
pick <- which(saudat[,"cpue"] > 0)
ccf(x=saudat[pick,"catch"],y=saudat[pick,"cpue"],type="correlation",
    ylab="Correlation",plot=TRUE)
ymax2 <- getmax(c(saudat[,"cpue"],ans$outmat[1:nce,"predCE"]))
plot(yrs,saudat[,"cpue"],type="l",lwd=2,ylim=c(0,ymax2),yaxs="i",ylab="cpue",
     panel.first=grid())
lines(yrs,ans$outmat[1:nce,"predCE"],lwd=2,col=2)

ans$msy


