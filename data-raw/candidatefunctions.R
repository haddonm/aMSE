


doprojection <- function(zoneCP,zoneDP,glob,ctrl,projyrs,applyHS,hsargs,
                         projpms) {
  zoneCP=zoneCP
  zoneDP=zoneDP
  glob=glb
  ctrl=ctrl
  projyrs=projC$projyrs
  applyHS=mcdahcr
  hsargs=hsargs
  projpms=cmcda
  # get important constants
  sigmaR <- ctrl$withsigR # needed to add recruitment variation
  sigmaF <- ctrl$withsigB
  npop <- glob$numpop
  nsau <- glob$nSAU
  Ncl <- glob$Nclass
  nyrs <- projyrs
  movem <- glob$move
  reps <- ctrl$reps
  sauindex <- glob$sauindex
  matb <- numeric(npop)
  origTAC <- colSums(zoneDP$catch[1,,]) # mean sum of catches in last year
  # now do replicates, updating saucatch and saucpue each year
  if (ctrl$randseedP > 0) set.seed(ctrl$randseedP) # set random seed if desired
  for (iter in 1:reps) {
    TAC <- origTAC[iter]
    for (year in 2:nyrs) { # iter=1; year=11
      catbysau <- catchbysau(inexpB=zoneDP$exploitB[(year - 1),,iter],sauindex,TAC) # no error initially
      multh <- apply(saucpue[1:(year-1),,1],2,applyHS,yr=(year-1)) # apply mcdahcr
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
    }   # year loop        zoneDR$matureB[,,1]
  }     # rep loop
  zoneDP$cesau <- saucpue
  return(zoneDP=zoneDP)
} # end of doprojection




doproj <- function(zoneC,zoneDD,glb,ctrl,projC,applyHS=mcdahcr,HSargs,
                   histpms) {
  ans <- calcprojsel(zoneC,projC,glb) # calculate selectivity for projections
  sel <- ans$projSel
  selwt <- ans$projSelWt
  reps <- ctrl$reps
  pyrs <- projC$projyrs
  zoneDP <- makezoneDP(pyrs,reps,glb,zoneDD) # make object to hold projections
  sigmar <- ctrl$withsigR # needed to add recruitment variation
  npop <- glb$numpop
  nsau <- glb$nSAU
  Ncl <- glb$Nclass

  lastyr <- glb$Nyrs
  histyrs <- projC$inityrs:lastyr
  movem <- glb$move
  reps <- ctrl$reps
  matb <- numeric(npop)
  origTAC <- sum(zoneDD$catch[lastyr,]) # Assumes the TAC is always caught!
  sauindex <- glb$sauindex
  saucatch <- array(0,dim=c(pyrs,nsau,reps))
  saucpue <- saucatch
  grad4 <- array(0,dim=c(pyrs,nsau,reps)) # declare arrays to store PMs
  grad1 <- grad4
  targsc <- grad4
  targetce <- numeric(nsau)
  for (iter in 1:reps) {
    for (year in 1:pyrs) {  # iter=1; year = 1
      ans <- getmcdadatap(zoneDD$catch,zoneDD$cpue,
                          zoneDP$catch[,,iter],zoneDP$cpue[,,iter],
                          histyrs,year)


    }
  }
  #calibrate the HCR on zoneDD before the projections
  for (yr in 1:inityrs) { # iter=1; yr=4
    saucatch[yr,,iter] <- tapply(zoneDP$catch[yr,,iter],sauindex,sum,na.rm=TRUE)
    wts <- zoneDP$catch[yr,,iter]/(saucatch[yr,sauindex,iter])
    saucpue[yr,,iter] <- tapply((zoneDP$cpue[yr,,iter] * wts),sauindex,sum,na.rm=TRUE)
    if (yr >= wid) {
      grad1[yr,,iter] <- apply(saucpue[1:yr,,iter],2,getgradone,yr=yr)
      grad4[yr,,iter] <- apply(saucpue[1:yr,,iter],2,getgradwid,yr=yr)
    }
  }
  for (sau in 1:nsau)
    targetce[sau] <- targscore(saucpue[1:inityrs,sau,iter],qnt=targqnt)$result
  targsc[1:inityrs,,iter] <- targetce

}









#' @title plotCPUE plots the scaled historic cpue against the predicted CPUe
#'
#' @description plotCPUE rescales the predicted CPUE from the OM to match the
#'     duration of each time-series of observed historical CPUE, currently it
#'     assumes the same number of years of CPUE for each observed series, which
#'     is wrong for blocks 6 and 13W!#'
#'
#' @param predCE the predicted cpue for each SAU from zoneD object
#' @param histCE the historical cpue for each SAU from the fishery
#' @param pickyr the indices of the years of overlap between the historical and
#'     predicted time-series of cpue
#' @param sau the SAU to be plotted. default=1:8 to suit the west coast of TAS
#'
#' @return invisibly the rescaled predCE matrix.
#' @export
#'
#' @examples
#' print("wait on some data sets and this function being adopted")
plotCPUE <- function(predCE,histCE,pickyr,sau=1:8) {
  plotprep(width=6,height=7,newdev=FALSE)
  parset(plots=c(4,2),margin=c(0.25,0.25,0.05,0.01),outmargin = c(1,1,0,0.4))
  for (i in sau) { # i = 1
    pCE <- predCE[,i]/mean(predCE[pickyr,i])
    ymax <- getmax(pCE)
    plot(yrs,pCE,type="l",lwd=2,ylim=c(0,ymax),xlab="",ylab="",col=4,
         panel.first=grid())
    lines(yrs[pickyr],scaletoOne(histCE[,i]),lwd=2,col="darkorange")
    mtext(text=zone1$SAUnames[i],side=3,line=-1.1,cex=1.0)
  }
  mtext("Predicted vs Observed CPUE",side=2,cex=1.0,outer=TRUE,line=-0.35)
  return(invisible(pCE))
} # end of plotCPUE

#' @title plotSAUdepl plots the depletion time-series for all SAU
#'
#' @description plotSAUdepl is a convenient plotting routine designed for the
#'     west coast of Tasmania. It is currently not a general SAU depletion
#'     plotting routine, but can potentially be generalized.
#'
#' @param outB the matrix of depletion values from the zoneD object
#'
#' @return nothing but it does plot a set of graphs
#' @export
#'
#' @examples
#'  print("wait on some data sets")
plotSAUdepl <- function(outB) {
  plotprep(width=7,height=6,newdev = FALSE)
  parset(plots=c(3,1))
  ymax <- getmax(outB[,c(1,2,3,8)])
  plot(yrs,outB[,1],type="l",lwd=2,ylim=c(0,ymax),panel.first=grid(),ylab="6,7,8,13")
  for (i in c(2,3,8)) lines(yrs,outB[,i],lwd=2,col=i)
  legend("topright",legend=c(6,7,8,13),col=c(1,2,3,8),lwd=3,bty="n",cex=1.0)
  ymax <- getmax(outB[,c(4,5,7)])
  plot(yrs,outB[,4],type="l",lwd=2,ylim=c(0,ymax),col=4,panel.first=grid(),
       ylab=c("9,10,12"))
  for (i in c(5,7)) lines(yrs,outB[,i],lwd=2,col=i)
  legend("topright",legend=c(9,10,12),col=c(4,5,7),lwd=3,bty="n",cex=1.0)
  ymax <- getmax(outB[,6])
  plot(yrs,outB[,6],type="l",lwd=2,ylim=c(0,ymax),col=6,panel.first=grid(),
       ylab="11")
  legend("topright",legend=c(11),col=c(6),lwd=3,bty="n",cex=1.0)
} # end of plotSAUdepl

# cpue <- zoneDRp$cpue[1:10,1,1]
# years <- 1:10
# getgrad1(cpue)
# getgrad4(cpue)
#
# plotprep()
# parset()
# plot1(years,cpue,lwd=2)



#condcpue <- MCDAdata(resdir,projC$HS,zone1$SAUnames)


vectce <- c(0.867862817,0.885256733,0.877541942,1.008936892,1.169236848,
            1.231699353,1.226890881,1.237515477,1.238744573,1.156193697,
            1.069769077,1.086021344,0.954323058,0.956561516,0.982457905,
            1.001020146,0.951562922,1.040046115,0.970610956,1.016168543,
            0.978065368,0.862886265,0.876245926,0.832028769,0.901777782,
            0.914997995,0.925753845,0.779823254)
names(vectce) <- 1992:2019

quantile(ce,probs=c(0.55))

one <- getgrad4(vectce)
two <- getgrad4(scaleto1(vectce))
cbind(one,two,two-one)

one <- targscore(vectce)
two <- targscore(scaleto1(vectce))
cbind(one$ans[,"targsc"],two$ans[,"targsc"],two$ans[,"targsc"]-one$ans[,"targsc"])


# hcr <- c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2)
# names(hcr) <- c(1:10)
# wid=4
# targqnt=0.55
# pmwts=c(0.65,0.25,0.1)
#
#
# mcda <- mcdahcr(vectce,only=FALSE)
#
# plot(vectce,type="l",ylim=c(0,1.4))
# abline(h=c(limrp,targ,uprp),col=c(2,3,2))
#
  # requirements for doprojection
  zoneCP=zoneCP;zoneDP=zoneDR;glob=glb;ctrl=ctrl;projyrs=projC$projyrs;inityrs=projC$inityrs;

  yr=28
  hcr <- c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2)
  names(hcr) <- c(1:10)
  wid=4
  targqnt=0.55
  pmwts=c(0.65,0.25,0.1)
  only=TRUE


library(microbenchmark)

microbenchmark(
mcdahcr(vectce,28),
oldmcdahcr(vectce,28)
)


zoneDP=zoneDRp
year=11; iter=1

ce <- zoneDP$cpue[1:(year - 1),,iter]
ce
out <- apply(ce,2,mcdahcr)

ces <- apply(ce,2,scaleto1)
outs <- apply(ces,2,mcdahcr)


library(aMSE)
library(rutilsMH)
library(makehtml)

resdir <- "C:/Users/User/Dropbox/A_Code/aMSEUse/conddata/generic"
ctrlzonetemplate(resdir,"control2.csv")


zone <- readctrlfile(resdir,infile="control2.csv")






x <- c(1,2,3,4,5)
wts <- c(1,1,1,1,1)

wtedmean(x,wts)


library(microbenchmark)

microbenchmark(
wtedmean(saucpue,saucatch),
wtmean(saucpue,saucatch)
)



reps=10
hr <- zoneDD$harvestR[1,]
ans <- matrix(0,nrow=reps,ncol=16)
for (pop in 1:16) ans[,pop] <- rnorm(reps,mean=hr[pop],sd=0.003)

ans





