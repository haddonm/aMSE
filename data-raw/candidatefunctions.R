

# Compare recovery rates vs catch --------------------------------------------

scenes <- result$scenes
nscenes <- length(scenes)
nsau <- glb$nSAU
saunames <- glb$saunames

matB <- array(0,c(glb$pyrs,nsau,nscenes),
               dimnames=list(glb$pyrnames,glb$saunames,scenes))
acatch <- array(0,c(glb$pyrs,nsau,nscenes),
              dimnames=list(glb$pyrnames,glb$saunames,scenes))

for (scen in 1:nscenes) {
  matB[,,scen] <- getmedbysau(result$dyn[[scen]]$deplsB,glb)
  acatch[,,scen] <- getmedbysau(result$dyn[[scen]]$acatch,glb)
}

yrs <- glb$pyrnames

legloc="bottomright"


#single
sau=1

plotprep(width=10,height=5)
parset()
maxx <- getmax(acatch[,sau,])
minx <- getmin(acatch[,sau,])
maxy <- getmax(matB[,sau,])
miny <- getmin(matB[,sau,])
plot(acatch[,sau,1],matB[,sau,1],type="l",lwd=2,xlim=c(minx,maxx),
     ylim=c(miny,maxy),panel.fist=grid(),xlab="",ylab="")
for (scen in 2:nscenes) lines(acatch[,sau,scen],matB[,sau,scen],lwd=2,col=scen)
mtext(saunames[sau],side=3,outer=FALSE,cex=1.1,line=-1.1,adj=0.05)
mtext("Aspirational Catch (t)",side=1,outer=TRUE,cex=1.1,line=-1.2)
mtext("Depletion of Mature Biomass (t)",side=2,outer=TRUE,cex=1.1,line=-1.2)
legend(legloc,legend=scenes,col=1:nscenes,lwd=3,bty="n",cex=1.1)



plotprep(width=9, height=10)
parset(plots=c(4,2),outmargin=c(1,1,0,0),margin=c(0.3,0.3,0.05,0.1))
for (sau in 1:nsau) {
  maxx <- getmax(acatch[,sau,])
  minx <- getmin(acatch[,sau,])
  maxy <- getmax(matB[,sau,])
  miny <- getmin(matB[,sau,])
  plot(acatch[,sau,1],matB[,sau,1],type="l",lwd=2,xlim=c(minx,maxx),
       ylim=c(miny,maxy),panel.fist=grid(),xlab="",ylab="")
  for (scen in 2:nscenes) lines(acatch[,sau,scen],matB[,sau,scen],lwd=2,col=scen)
  mtext(saunames[sau],side=3,outer=FALSE,cex=1.1,line=-1.1,adj=0.05)
}
mtext("Aspirational Catch (t)",side=1,outer=TRUE,cex=1.1,line=-0.2)
mtext("Depletion of Mature Biomass (t)",side=2,outer=TRUE,cex=1.1,line=-0.2)
legend(legloc,legend=scenes,col=1:nscenes,lwd=3,bty="n",cex=1.1)



c1 <- dyn[[1]]$deplsB[59:88,1,]
c2 <- dyn[[2]]$acatch[59:88,1,]

ymax=getmax(c(c1,c2))
plotprep(width=10,height=5)
parset()
plot(yrs,c1[,1],type="l",lwd=1,col=RGB("black",alpha=100),ylim=c(0,ymax),
     xlab="year",ylab="Mature Biomass Depletion")
for (i in 2:50) lines(yrs,c1[,i],lwd=1,col=RGB("black",alpha=100))
for (i in 1:50) lines(yrs,c2[,i],lwd=1,col=RGB("red",alpha=100))
legend("topleft",scenes,col=c(1,2),lwd=4,bty="n",cex=1.5)




depl <- projectiononly(sauout$deplsB,glb)

dep6 <- apply(depl[,1,],1,median)


cdepl <- projectiononly(result$dyn[[1]]$deplsB,glb)

cdep6 <- apply(cdepl[,1,],1,median)

cbind(dep6,cdep6)









cols <- c("#FF0000","#FF8000","#FFFF00","#80FF00","#00FF00","#00FF80","#00FFFF",
          "#0080FF","#0000FF","#8000FF","#FF00FF","#FF0080")

popcol <- c("red","#FF7790","orange","yellow","#80FE90","green","cyan",
            "#0080FF","blue","#8000FF","magenta","#FF0080")


plotprep(width=10,height=5)
parset()
plot(1:100,rep(1,100),type="l",lwd=2,col=popcol[1],ylim=c(0,13))
for (i in 2:12) lines(1:100, rep(i,100),lwd=2,col=popcol[i])
for (i in 1:12)  lines(1:100, rep((i+0.2),100),lwd=2,col=cols[i])





# headtail


headtail <- function(x) { # x <- sauout$matureB[,,1]
  if ("tbl" %in% class(x)) {
    x <- as.data.frame(unclass(x), stringsAsFactors = FALSE,
                           check.names=FALSE)
  }
  if (length(dim(x)) > 2) x <- x[,,1]
  numcol <- ncol(x)
  outans <- rbind(head(x),rep(NA,numcol),tail(x))
  return(outans)
}


# scenarioproperties


filename <- fullfiles[4]
load(filename)


outhcr <- out$outhcr
g4s <- outhcr$g4s
yrs <- as.numeric(dimnames(g4s)[[1]])

plotprep(width=8, height=4.5)
parset()
plot(yrs,g4s[,7,1],type="l",lwd=2)
lines(yrs,g4s[,7,2],lwd=2,col=2)
lines(yrs,g4s[,7,3],lwd=2,col=4)

expB <- out$sauout$exploitB
yrsindex <- 59:88
yrs <- as.numeric(dimnames(expB)[[1]])[yrsindex]
catch <- out$sauout$catch
cpue <- out$sauout$cpue

iny1 <- expB
iny2 <- catch
iny3 <- cpue
plotprep(width=8, height=8)
parset(plots=c(3,1),margin=c(0.3,0.45,0.05,0.05))
plot(yrs,iny1[yrsindex,7,1],type="l",lwd=2,xlab="")
lines(yrs,iny1[yrsindex,7,2],lwd=2,col=2)
lines(yrs,iny1[yrsindex,7,3],lwd=2,col=4)
plot(yrs,iny2[yrsindex,7,1],type="l",lwd=2,xlab="")
lines(yrs,iny2[yrsindex,7,2],lwd=2,col=2)
lines(yrs,iny2[yrsindex,7,3],lwd=2,col=4)
plot(yrs,iny3[yrsindex,7,1],type="l",lwd=2,xlab="")
lines(yrs,iny3[yrsindex,7,2],lwd=2,col=2)
lines(yrs,iny3[yrsindex,7,3],lwd=2,col=4)



   result <- do_comparison(rundir=rundir,postfixdir=postfixdir,outdir=outdir,
                           files=files,pickfiles=c(1,2,3),verbose=TRUE)
   out <- do_comp_outputs(result,projonly=TRUE)
   catch <- scenebyvar(dyn=out$dyn,byvar="catch",glb=out$glbc[[1]])
   catqnts <- sauquantbyscene(catch,out$glbc[[1]])
   sauribbon("",scenes=out$scenes,sau=8,varqnts=catqnts,glb=out$glbc[[1]],
   varname="Catch",console=TRUE,q90=TRUE,intens=100,addleg="bottomright")




 #' @title plotscenproj plots a 3D array of data by sau in a single plot
 #'
 #' @description plotsceneproj uses a 3D array of years x sau x replicates to
 #'     produce a plot by sau of those projections all in a single plot. Options
 #'     exist to include a horizontal dashed line and to include a median line
 #'     and either 90 or 95 quantiles across the plot.
 #'
 #' @param rundir the directory in which all results are held for a scenario or
 #'     comparison of scenarios
 #' @param inarr the 3D array of projections for a given scenario derived from
 #'     catch or cpue or whatever
 #' @param glb the global object relating to the particular acenario
 #' @param scene the scenario name
 #' @param filen the filename in which to store the plot, default="" which draws
 #'     the plot to the console
 #' @param label what is the name of the variable being plotted, default=""
 #' @param maxy what constant maximum y-axis to use? default=0 which uss the
 #'     maximum for each sau
 #' @param Q which quantile to use, default = 90. If set to 0 no median or
 #'     quantiles are plotted, 95 will lead to the 95 percent quantiles being
 #'     plotted
 #' @param hline should a horizontal dashed line be included. default=NA, which
 #'     means that no line is added. Otherwise whichever value is given this
 #'     argument will lead to a dashed horizontal black line of width 1.
 #'
 #' @return invisibly a list of the quantiles for each sau
 #' @export
 #'
 #' @examples
 #' print("wait on data sets")
 #' #  rundir=rundir; inarr=out$sauout$deplsB[58:88,,];glb=out$glb;scene="BCmeta"
 #' #  filen="";label="Mature Biomass Depletion";maxy=0;Q=0;hline=0.2; eg=3
 #' #
plotsceneproj <- function(rundir,inarr,glb,scene,filen="",label="",
                         maxy=0,Q=90,hline=NA,eg=0) {
  nsau <- glb$nSAU
  saunames <- glb$saunames
  outmed <- makelist(saunames)
  egtraj <- NULL
  reps <- glb$reps
  if (nchar(filen) > 0) {
   filen <- filenametopath(rundir,filen)
   caption <- paste0("Projections of ",label," for ",scene," for each SAU.")
  }
  yrs <- as.numeric(names(inarr[,1,1]))
  plotprep(width=8, height=8,newdev=FALSE,filename = filen,verbose=FALSE)
  parset(plots=pickbound(nsau),margin=c(0.3,0.4,0.05,0.1),outmargin=c(0,1,0.25,0),
        byrow=FALSE)
  for (sau in 1:nsau) { # sau=1; i = 1
     dat <- inarr[,sau,]
     meds <- apply(dat,1,quants)
     outmed[[sau]] <- meds
     ymax <- maxy
     if (ymax == 0) ymax <- getmax(dat)
     plot(yrs,dat[,1],type="l",lwd=1,col="grey",panel.first=grid(),
          ylim=c(0,ymax),ylab=saunames[sau],xlab="",yaxs="i")
   for (i in 1:reps) lines(yrs,dat[,i],lwd=1,col="grey")
   if (eg > 0) {
     trajeg <- sort(trunc(runif(eg,min=1,max=reps)))
     for (i in 1:eg) lines(yrs,dat[,trajeg[i]],lwd=2,col=2)
     egtraj <- rbind(egtraj,trajeg)
   }
   if (Q > 0) {
       lines(yrs,meds["50%",],lwd=2,col=4)
       if (Q == 90) {
         lines(yrs,meds["5%",],lwd=1,col=2)
         lines(yrs,meds["95%",],lwd=1,col=2)
       } else {
         lines(yrs,meds["2.5%",],lwd=1,col=2)
         lines(yrs,meds["97.5%",],lwd=1,col=2)
       }
   }
   if (!is.na(hline)) abline(h=hline,lwd=1,col="black",lty=2)
  }
  if (eg > 0) rownames(egtraj) <- saunames; colnames(egtraj) <- 1:eg
  mtext(paste0(scene,"  ",label),side=2,line=-0.2,outer=TRUE,cex=1.1)
  if (nchar(filen) > 0)
   addplot(filen=filen,rundir=rundir,category="C_vs_MSY",caption)
  return(invisible(list(outmed=outmed, egtraj=egtraj)))
} # end of plotsceneproj


plotsceneproj(rundir,out$deplsB,glb=out$glb,scene="BCmeta",filen="",
              label="Mature Biomass Depletion",maxy=0,Q=90,hline=NA)




#' @title onecattraj is a one category trajectory plot of yrs x variable
#'
#' @description onecattraj expects a 2D matrix of yrs x reps of some variable.
#'     For example, if one has 50 years x 250 replicates of the depletion level
#'     of spawning biomass then dat would be dat[1:50,1:250]. Options exist
#'     for including quantiles of the median trajectory plus either the inner
#'     90th or 95th quantiles, designated by setting Q = 90 or 95. If Q is left
#'     = 0 no quantiles are plotted. If eg is set > 0 then eg random replicates
#'     will be selected and plotted as example trajectories. defpar defaults to
#'     TRUE which will define the plot characteristics. To override this with
#'     expernal definitions set defpar = FALSE
#'
#' @param dat a 2d matrix of yrs x replicates of some dynamic variable
#' @param label The y-axis label for the plot, default = 'ylabe
#' @param maxy  allows one to define the maximum y value for the plot. The
#'     default = 0, which means the maximum of he data will be used.
#' @param miny default = 0, which means the plot with have a y-axis starting at
#'     zero. This argument provides the option of changing this, which would be
#'     important if there are negative values to your dynamic variable. If miny
#'     is set then maxy must also be set
#' @param rundir the directory in which to store the plot file if a filename is
#'     given. This is only required if filen is defined.
#' @param filen the full name of the file
#' @param Q valid values can be 0, 90, or 95, meaning no CI, add the 90th CI,
#'     or add the 05th CI in red to the blue to the median line.
#' @param eg Should example single trajectories be added to the plot? More than
#'     5 tends to be very busy, 3 or less are good. More than 2 appears to many
#'     when the CI are also plotted. If plotted long they are red, if plotted
#'     with the CI these lines are yellow.
#' @param defpar default = TRUE, which defines the plot properties. If you wish
#'     to add multiple plots to a single graph then define your own pars
#'     outside the function and set defpar = FALSE
#' @param xlab xlabel omitted if the x-axis is years, otherwise name it here
#'     and the plot will have room added for the xlabel.
#' @param hline should a horizontal dashed line be included. default=NA, which
#'     means that no line is added. Otherwise whichever value is given this
#'     argument will lead to a dashed horizontal black line of width 1.
#'
#' @return a list of the median values and a vector of the indices of the
#'     random trajectories selected if added.
#' @export
#'
#' @examples
#'  delta <- 2+rnorm(100,mean=0.02,sd=0.05)
#'  y <- matrix(0,nrow=101,ncol=100)
#'  for (i in 1:100) y[,i] <-sin(seq(0,delta[i]*pi,length=101))
#'  rownames(y) <- seq(0,2*pi,length=101)
#'  catout <- onecattraj(dat=y,label="Randomized sine waves",maxy=1.1,miny=-1.1,
#'                       rundir=NA,filen="",Q=95,eg=0,xlab="pi-value")
#'  catout <- onecattraj(dat=y,label="Randomized sine waves",maxy=1.1,
#'                       miny=-1.1,rundir=NA,filen="",Q=0,eg=3,xlab="pi-value")
onecattraj <- function(dat,label="ylabel",maxy=0,miny=0,rundir=NA,filen="",
                       Q=0,eg=0,defpar=TRUE,xlab="",hline=NA) {
  yrs <- as.numeric(rownames(dat))
  reps <- ncol(dat)
  meds <- apply(dat,1,quants)
  if (miny != 0) {
    yax <- c(miny,maxy)
   } else {
     if (maxy > 0) {
       yax <- c(0,maxy)
       } else {
       yax <- c(0,getmax(dat))
     }
  }
  filelen <- nchar(filen)
  if (filelen > 0) {
    if (substr(filen,filelen-2,filelen) != "png") filen <- paste0(filen,".png")
    filen <- pathtopath(rundir,filen)
    caption <- paste0("Projections of ",label)
  }
  if (defpar) {
    plotprep(width=8, height=4.5,newdev=FALSE,filename = filen,verbose=FALSE)
    xmarg <- ifelse(nchar(xlab > 0),0.45,0.3)
    parset(plots=c(1,1),margin=c(xmarg,0.4,0.05,0.1))
  }
  plot(yrs,dat[,1],type="l",lwd=1,col="grey",panel.first=grid(),
       ylim=yax,ylab=label,xlab=xlab,yaxs="i")
  for (i in 1:reps) lines(yrs,dat[,i],lwd=1,col="grey")
  if (Q > 0) {
    lines(yrs,meds["50%",],lwd=2,col=4)
    if (Q == 90) {
      lines(yrs,meds["5%",],lwd=2,col=2)
      lines(yrs,meds["95%",],lwd=2,col=2)
    } else {
      lines(yrs,meds["2.5%",],lwd=2,col=2)
      lines(yrs,meds["97.5%",],lwd=2,col=2)
    }
  }
  if (!is.na(hline)) abline(h=hline,lwd=1,col="black",lty=2)
  egtraj <- NULL
  if (eg > 0) {
    lincol <- ifelse(Q > 0,"yellow","red")
    egtraj <- as.matrix(sort(trunc(runif(eg,min=1,max=reps))))
    colnames(egtraj) <- label; rownames(egtraj) <- 1:eg
    for (i in 1:eg) lines(yrs,dat[,egtraj[i]],lwd=2,col=lincol)
  }
  return(invisible(list(meds=meds,egtraj=egtraj)))
} # end of onecatproj

catout <- onecattraj(dat=out$sauout$deplsB[58:88,6,],label="sau12 - BCmeta",maxy=0,
           rundir=NA,filen="",Q=95,eg=2,hline=0.2)

catout <- onecattraj(dat=y,label="Randomized sine waves",maxy=1.1,miny=-1.1,
                     rundir=NA,filen="",Q=95,eg=0,xlab="pi")

catout <- onecattraj(dat=y,label="Randomized sine waves",maxy=1.1,miny=-1.1,
                     rundir=NA,filen="",Q=0,eg=3,xlab="pi")

catout



dat3 <- out$sauout$deplsB[59:88,,];label="Mature Biomass Depletion"
catnames=out$glb$saunames;maxy=0;miny=0:rundir=NA;filen="";Q=0;eg=0;xlab="";wide=8
hline=NA


#' @title multicattraj plots a 3D array of data by category in a single graphic
#'
#' @description multicattraj uses a 3D array of years x sau x replicates
#'     'steps x category x replciates' to produce a plot by category of those
#'     replicates all in a single plot. Options exist to include a horizontal
#'     dashed line and to include a median line with either 90 or 95 quantiles,
#'     and or a highlighted set of randomly chosen replicatges to illustrate
#'     individual variation.
#'
#' @param dat3 a 3D array of steps 'years", x category 'sau' x replicates
#' @param label The y-axis label for the plot, default = 'ylabe
#' @param catnames names for each category, eg the sau names in aMSE
#' @param maxy  allows one to define the maximum y value for the plot. The
#'     default = 0, which means the maximum of he data will be used.
#' @param miny default = 0, which means the plot with have a y-axis starting at
#'     zero. This argument provides the option of changing this, which would be
#'     important if there are negative values to your dynamic variable. If miny
#'     is set then maxy must also be set
#' @param rundir the directory in which to store the plot file if a filename is
#'     given. This is only required if filen is defined.
#' @param filen the full name of the file
#' @param Q valid values can be 0, 90, or 95, meaning no CI, add the 90th CI,
#'     or add the 05th CI in red to the blue to the median line.
#' @param eg Should example single trajectories be added to the plot? More than
#'     5 tends to be very busy, 3 or less are good. More than 2 appears to many
#'     when the CI are also plotted. If plotted long they are red, if plotted
#'     with the CI these lines are yellow.
#' @param defpar default = TRUE, which defines the plot properties. If you wish
#'     to add multiple plots to a single graph then define your own pars
#'     outside the function and set defpar = FALSE
#' @param xlab xlabel omitted if the x-axis is years, otherwise name it here
#'     and the plot will have room added for the xlabel.
#' @param hline should a horizontal dashed line be included. default=NA, which
#'     means that no line is added. Otherwise whichever value is given this
#'     argument will lead to a dashed horizontal black line of width 1.
#' @param wide how wide, in inches, should the graphic be, default = 8
#' @param high how tall, in inches, should the graphic be, default = 8
#'
#' @return a list of the quantiles and random trajectories, if selected, for
#'     each category
#' @export
#'
#' @examples
#'  y <- array(0,dim=c(101,4,100),dimnames=list(seq(0,2*pi,length=101),1:4,1:100))
#'  for (j in 1:4) {  # make some 3D data
#'    delta <- 2+rnorm(100,mean=0.02,sd=0.05)
#'    for (i in 1:100) y[,j,i] <-sin(seq(0,delta[i]*pi,length=101))
#'  }
#'  multicattraj(dat3=y,label="Randomized Sines",catnames=c("A","B","C","D"),
#'               maxy=1.1,miny=-1.1,rundir=NA,filen="",
#'               Q=95,eg=0,hline=NA,xlab="pi-value",wide=8,high=6)
multicattraj <- function(dat3,label="ylabel",catnames,maxy=0,miny=0,rundir=NA,
                         filen="",Q=0,eg=0,xlab="",hline=NA,wide=8,high=8) {
  datdim <- dim(dat3)
  ninst <- datdim[1]
  ncat <- datdim[2]
  inst <- as.numeric(dimnames(dat3)[[1]])
  outmulti <- makelist(catnames)
  filelen <- nchar(filen)
  if (filelen > 0) {
    if (substr(filen,filelen-2,filelen) != "png") filen <- paste0(filen,".png")
    filen <- pathtopath(rundir,filen)
    caption <- paste0("Projections of ",label)
  }
  plotprep(width=wide,height=high,newdev=FALSE,filename = filen,verbose=FALSE)
  xmarg <- ifelse(nchar(xlab) > 0,0.45,0.3)
  parset(plots=pickbound(ncat),margin=c(xmarg,0.4,0.05,0.1),
         outmargin=c(0,1,0.25,0),byrow=FALSE)
  for (cat in 1:ncat) { # cat=1; i = 1
    dat <- dat3[,cat,]
    outmulti[[cat]] <- onecattraj(dat=dat,label=paste0(catnames[cat],"_",label),
                         maxy=maxy,miny=miny,rundir=rundir,filen=filen,
                         Q=Q,eg=eg,defpar=FALSE,xlab=xlab,hline=hline)
  }
  return(invisible(outmulti))
} # en do fmulticattraj


multicattraj(dat3=out$sauout$deplsB[55:88,,],label="MatureB Depletion",
              catnames=out$glb$saunames,maxy=0,miny=0,rundir=NA,filen="",
              Q=95,eg=0,hline=NA,xlab="",wide=8,high=9)




# compare hsargs-----------------------------------



comptashsargs <- function(ans) {
  label <- names(ans)
  hsargs <- makelist(label)
  nscen <- length(label)
  for (i in 1:nscen) hsargs[[i]] <- ans[[i]]$hsargs
  nargs <- length(hsargs[[1]])
  names(ans[[1]])
}



hlines=list(catch=outprod[,"MSY"],spawnB=outprod[,"Bmsy"],harvestR=0,
            cpue=outprod[,"CEmsy"])




# find sauMSY-----------------------------


findsaumsy <- function(product,glb) {  # product=prody
  harv <- as.numeric(rownames(product[,"Catch",]))
  nh <- length(harv)
  nsau <- glb$nSAU
  saunames <- glb$saunames
  sauindex <- glb$sauindex
  label <- dimnames(product)
  sauyield <- matrix(0,nrow=nh,ncol=nsau,dimnames=list(label[[1]],saunames))
  saumatB <- matrix(0,nrow=nh,ncol=nsau,dimnames=list(label[[1]],saunames))
  for (i in 1:nh) {
    saumatB[i,] <- tapply(product[i,"MatB",],sauindex,sum,na.rm=TRUE)
    sauyield[i,] <- tapply(product[i,"Catch",],sauindex,sum,na.rm=TRUE)
  }

  catch <- product[,"Catch",]
  numpop <- ncol(catch)
  label <- c(colnames(product),"index")
  xval <- matrix(0,nrow=numpop,ncol=length(label),
                 dimnames=list(1:numpop,label))
  for (pop in 1:numpop) { # pop=1
    pick <- which.max(catch[,pop])
    xval[pop,] <- c(product[pick,,pop],pick)
  }
  return(xval)
}





comp <- getfilestocompare(outdir=outdir,filenames=files[c(12,2,16)],
                          altlabel=c("TT_BC","TV_BC","TS_BC"),verbose=TRUE,
                          listtoenv=TRUE)



# population properties-------------------------

# run aMSE ---------------------------
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
prefixdir <- pathtopath(dropdir,"/A_codeUse/aMSETAS/NW/")
#prefixdir <- pathtopath(dropdir,"/A_codeUse/aMSEUse/scenarios/")


postfixdir <- "sau56"
verbose <- TRUE
rundir <- path.expand(filenametopath(prefixdir,postfixdir))



biol <- read.csv(pathtopath(rundir,"zonebiology.csv"),header=TRUE)

pickP <- which(biol[,"MSY"] < 50)
biol1 <- biol[pickP,]

plotprep(width=10, height=5)
parset()
plot(biol1[,"MaxDL"],biol1[,"MSY"],type="p",pch=16,cex=1.0)

model <- lm(biol1[,"MSY"] ~ biol1[,"MaxDL"])
abline(model,lwd=3,col=2)
summary(model)


var1 <- "MaxDL"
var2 <- "MSYDepl"

plotprep(width=10, height=5)
parset()
plot(biol1[,var1],biol1[,var2],type="p",pch=16,cex=1.0,xlab=var1,ylab=var2)
model <- lm(biol1[,var2] ~ biol1[,var1])
summary(model)
abline(model,lwd=3,col=2)

pairs(biol1[,c(1:3,5:10,12)],pch=16,gap=0.2)


MSYDepl ~ MaxDL


# makePOPmatrix -----------------------------------------------

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
prefixdir <- pathtopath(dropdir,"/A_codeUse/aMSETAS/NW/")
#prefixdir <- pathtopath(dropdir,"/A_codeUse/aMSEUse/scenarios/")


postfixdir <- "sau56"
verbose <- TRUE
rundir <- path.expand(filenametopath(prefixdir,postfixdir))
ctrlfile <- "controlsau56.csv"


draft <- getdraftpops(rundir=rundir,ctrlfile=ctrlfile)


round(draft$pops[,1:8],4)




filename <- pathtopath(rundir,"popconstants.csv")
write.csv(draft$pops,filename)


msy <- sapply(out$zoneC,"[[","MSY")
bLML <- sapply(out$zoneC,"[[","bLML")
B0 <- sapply(out$zoneC,"[[","B0")




displaypopprops <- function(rundir,x,verbose) {
  pops <- draft$pops
  out <- x$out
  condC <- x$condC
  glb <- out$glb
  ctrl <- x$ctrl
  zoneD <- x$out$zoneD
  zoneC <- x$out$zoneC
  setuphtml(rundir)

  notes <- c(paste0("Populations = ",glb$numpop),paste0("SAU = ",glb$nSAU),
             paste0("Randomseed for conditioning = ",ctrl$randseed),
             paste0("Recruitment variability sigR         = ",ctrl$withsigR),
             paste0("Exploitable Biomass variability sigB = ",ctrl$withsigB),
             paste0("Catch-Rate variability sigCE         = ",ctrl$withsigCE))

  label <- paste0(ctrl$runlabel," population properties. This is a bigtable",
                  " so, if you scroll across then use the topleft arrow to",
                  " return to the home page.")
  addtable(round(pops,4),filen="popprops.csv",rundir=rundir,category="poptable",
           caption=label,big=TRUE)
  # OrigComp tab--------------------------------------------
  saucompdata(allcomp=condC$compdat$lfs,glb=glb,horizline=140,console=FALSE,
              rundir=rundir,ylabel="Size-Composition of Catches",
              tabname="OrigComp")
  # What size comp data is there
  palfs <- condC$compdat$palfs
  label <- paste0(ctrl$runlabel," Number of observations of numbers-at-size in ",
                  "the catch for each SAU.")
  addtable(palfs,filen="sizecompnumbers.csv",rundir=rundir,category="OrigComp",
           caption=label)
  # plot initial equilibrium size-comp by population
  Nt <- zoneD$Nt[,1,]
  rownames(Nt) <- glb$midpts
  colnames(Nt) <- paste0(glb$saunames[glb$sauindex],"_",1:glb$numpop)
  draftnumbersatsize(rundir, glb, Nt, ssc=5)


  make_html(replist = NULL,  rundir = rundir,
            controlfile=ctrl$controlfile, datafile=ctrl$datafile, hsfile=NULL,
            width = 500, openfile = TRUE,  runnotes = notes,
            verbose = verbose, packagename = "aMSE",  htmlname = ctrl$runlabel)
} # end of displaypopprops




displaypopprops(rundir,draft,verbose=FALSE)





library(codeutils)

lmls <- c(2000,145,210,2001,145,210,2002,145,210,2003,145,210,2004,145,185,
          2005,150,185,2006,150,185,2007,150,185,2008,145,210,2009,145,210,
          2010,150,185,2011,150,185,2012,145,210)
numrow <- length(lmls)/3
projL <- matrix(lmls,nrow=numrow,ncol=3,
                dimnames=list(1:numrow,c("year","LML","Max")),byrow=TRUE)
projL
outp <- uniquepairs2(x=projL,col1=2,col2=3,yrs="year")
outp

yearcomb <- outyrs


ngrp <- length(yearcomb)
nobs <- sapply(yearcomb,length)
numrow <- sum(nobs) - 2*ngrp
columns <- c("year","LML","maxLML")
LMLtable <- matrix(0,nrow=numrow,ncol=length(columns))
colnames(LMLtable) <- columns
count <- 0
begin <- 1
for (i in 1:ngrp) { # i = 2
  num <- nobs[i]
  val <- yearcomb[[i]]
  lml <- val[1]
  maxlml <- val[2]
  yrs <- val[3:num]
  nyr <- length(yrs)
  count <- count + nyr
  LMLtable[begin:count,1] <- yrs
  LMLtable[begin:count,2] <- rep(val[1],length(begin:count))
  LMLtable[begin:count,3] <- rep(val[2],length(begin:count))
  begin <- begin+nyr
}
LMLs <- LMLtable[order(LMLtable[,"year"]),]
rownames(LMLs) <- LMLs[,"year"]


for (i in 1:ngrp) {
  wrk <- yearcomb[[i]]
  n <- length(wrk)
   yearvect <- c(yearvect,)

}

combnames <- names(outyrs)





popgrowth(rundir="",zoneC=out$zoneC,glb=out$glb,verbose=TRUE,console=TRUE,
          maxage=30,startsize= 2.0)






