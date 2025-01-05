


options("show.signif.stars"=FALSE,"stringsAsFactors"=FALSE,
        "max.print"=50000,"width"=240)
suppressPackageStartupMessages({
  library(aMSE)
  library(TasHS)
  library(codeutils)
  library(hplot)
  library(makehtml)
  library(knitr)
})
dropdir <- getDBdir()
prefixdir <- pathtopath(dropdir,"A_codeR/aMSEGuide/runs/")
postfixdir <- "EGcompare"
rundir <- pathtopath(prefixdir,postfixdir)
outdir <- "C:/aMSE_scenarios/EG/"
files <- dir(outdir)
#printV(files)
pickfiles=c(2,3,5)
scencol=c(1,2,4)
rundir=rundir;postfixdir=postfixdir;outdir=outdir;files=files;
verbose=TRUE; intensity=100; zero=FALSE; altscenes=c("NoMR","MR12","MR123")
juris="";ribbonleg="topleft"; Q90=TRUE;
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

nscen <- length(scenes)
glb <- glbc[[1]]
saunames <- glb$saunames
nsau <- glb$nSAU
projyrs <- glb$pyrnames
pyrs <- glb$pyrs

catmult <- makelist(scenes)
medcatmult <- array(0,dim=c(pyrs,nscen,nsau),
                    dimnames=list(projyrs,scenes,saunames))
for (scen in 1:nscen) { # scen = 1
  catm <- scores[[scen]]$catchmult
  catmult[[scen]] <- catm
  for (sau in 1:nsau)
    medcatmult[,scen,sau] <- apply(catm[,sau,],1,median)
}


# When catch multiplier first gets >=1.0, where TACs should stabilize or increase

mult1 <- matrix(0,nrow=nsau,ncol=nscen,dimnames=list(saunames,scenes))
for (scen in 1:nscen) {
  meds <- medcatmult[,scen,]
  for (sau in 1:nsau) {
    pick <- which(meds[,sau] >= 1.0)
    mult1[sau,scen] <- pick[1]
  }
}
mult1


poss <- c(0.25,0.75,0.8,0.85,0.9,1.0,1.05,1.1,1.15,1.2,1.25,1.3)
nposs <- length(poss)
dummy <- 1:12
countmult <- array(0,dim=c(nscen,nsau,nposs),
                   dimnames=list(scenes,saunames,poss))
invar <- catmult
filen <- ""
plotprep(width=12,height=6,newdev=FALSE,filename=filen,verbose=FALSE)
parset(plots=c(nscen,nsau),margin=c(0.25,0.25,0.05,0.05),
       outmargin = c(2.2,1.2,0.1,0.1),byrow=TRUE)
for (scen in 1:nscen) {
  for (sau in 1:nsau) {
    # scen = 1; sau=1
    sauC <- invar[[scen]][mult1[sau,scen],sau,]
    intnum <- numeric(nposs)
    counts <- table(sauC)
    index <- as.numeric(names(counts))
    pickpos <- match(index,poss)
    intnum[pickpos] <- counts
    intHout <- inthist(cbind(dummy,intnum),col=2,border=3,
                                    xaxis=FALSE,xmin=1,xmax=12)
    countmult[scen,sau,] <- intHout[,"counts"]
    axis(1,at=1:nposs,labels=poss)
    if (sau == 1) mtext(scenes[scen],side=2,line=1.25,cex=1.0)
    if (scen == nscen) mtext(saunames[sau],side=1,line=1.25,cex=1)
  }
}
label <- "Catch Multiplier Across Replicates when Median first => 1)"
mtext(label,side=1,outer=TRUE,line=0.75,cex=1.25)



startyr <- matrix(0,nrow=nsau,ncol=nscen,dimnames=list(saunames,scenes))
for (scen in 1:nscen) { # scen = 2
  for (sau in 1:nsau) {
    vectmult <- medcatmult[,scen,sau]
    pick1 <- which(vectmult == 1)
    np <- length(pick1)
    diffmult <-  pick1[2:np]- pick1[1:(np-1)]
    if (any(diffmult > 1)) {
      gt1 <- which(diffmult > 1)
      firstyr <- as.numeric(names(gt1[length(gt1)])) + 1
    } else {
      firstyr <- projyrs[pick1[1]]
    }
    startyr[sau,scen] <- firstyr
  }
}
startyr


# use of aav--------------------------------------------------

str(dyn[[1]])


nscen <- length(scenes)
glb <- glbc[[1]]
saunames <- glb$saunames
nsau <- glb$nSAU
projyrs <- glb$pyrnames
reps <- glb$reps


#aavbyrep <- matrix(0,nrow=reps,ncol=nsau,dimnames=list(1:reps,saunames))
aavbyrep <- array(0,dim=c(nscen,nsau,reps),
                  dimnames=list(scenes,saunames,1:reps))
for (scen in 1:nscen) { # scen = 1
  scenrep <- dyn[[scen]]$catch
  for (sau in 1:nsau) {
    aavbyrep[scen,sau,] <- apply(scenrep[59:70,sau,],2,getaav)
  }
}
str(aavbyrep)


getlim <- function(invar,inc=2) { #  invar = aavbyrep[1,6,]; inc=2
  rge <- range(invar)
  tmp <- trunc(rge[1]/inc)
  low <- tmp * inc
  tmp <- ceiling(rge[2]/inc)
  upp <- tmp * inc
  return(c(low,upp))
}



filen <- ""
plotprep(width=12,height=7,newdev=FALSE,filename=filen,verbose=FALSE)
parset(plots=c(nscen,nsau),margin=c(0.25,0.25,0.05,0.05),
       outmargin = c(2.2,2.2,0.1,0.1),byrow=TRUE)
for (scen in 1:nscen) {
  for (sau in 1:nsau) {
    # scen = 1; sau=1
    limx <- getlim(aavbyrep[,sau,])
    aavs <- aavbyrep[scen,sau,]
    hist(aavs,breaks=seq(limx[1],limx[2],2),main="",xlab="",ylab="",xlim=limx)
    # dmed <- density(sauC[mult1[sau,scen],])
    # plot(dmed,lwd=2,main="",xlab="",ylab="",xlim=c(0.2,1.75))
    # abline(v=c(5,6),lwd=1,col=2)
    mtext(round(mean(aavs),2),side=3,line=-1.1,cex=0.9,adj=1)
    mtext(round(sd(aavs),2),side=3,line=-2.1,cex=0.9,adj=1)
    if (sau == 1) mtext(scenes[scen],side=2,line=1.25,cex=1.0)
    if (scen == nscen) mtext(saunames[sau],side=1,line=1.25,cex=1)
  }
}
mtext("AAV by SAU",side=1,outer=TRUE,line=0.5,cex=1.2)
mtext("Tasmanian HS Scenario",side=2,outer=TRUE,line=0.75,cex=1.2)






# AAV explorations---------------------------------

devmovav <- function(invar,n=5,centred=TRUE) { # invar = catch; n= 5;centred=TRUE
  mov <- movav(invar,n=n,centred=centred)
  pick <- which(!is.na(mov))
  dev <- catch[pick] - mov[pick]
  return(list(dev=dev,mov=mov,pick=pick))
}


inc <- seq(0,30,1)
osc <- sin(inc)
posc <- osc + 1.0
main <- c(1:15,15,15:1)
mult <- 5
catch <- (mult * main) + posc
#trend <- loess(catch ~ inc,span=0.5)
outdev <- devmovav(catch,n=3)   #(catch - trend$fitted)
trend <- outdev$mov
dev <- outdev$dev
pick <- outdev$pick
detrend <- dev - min(dev)
aosc <- osc[pick] - min(dev)
#tmp <- loess(detrend ~ inc,span=0.75)
incp <- inc[pick]

plotprep(width=9, height=8)
parset(plots=c(2,1))
plot1(inc,catch,defpar=FALSE)
lines(inc,trend,lwd=2,col=3)
plot1(incp,detrend,lwd=2,col=2,defpar=FALSE)
lines(incp,posc[pick],lwd=2,col=4)
#lines(inc,tmp$fitted,lwd=3,col="purple")

getaav(catch)
getaav(posc)
getaav(detrend)
getaav(trend[pick])



getaav(tmp$fitted)



getaav(tmp$fitted)















