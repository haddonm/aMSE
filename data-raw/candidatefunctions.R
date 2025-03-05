

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


postfixdir <- "sau22"
verbose <- TRUE
rundir <- path.expand(filenametopath(prefixdir,postfixdir))
ctrlfile <- "controlsau22.csv"
controlfile=ctrlfile

draft <- getdraftpops(rundir=rundir,ctrlfile=ctrlfile)

pops <- as.matrix(draft$pops)
round(pops[,1:12],4)

tapply(pops[,"msy"],pops[,"SAU"],sum)

str2(draft)

# filename <- pathtopath(rundir,"popconstants.csv")
# write.csv(draft$pops,filename)



displaypopprops(rundir,draft,verbose=FALSE)

#pops

pops <- as.matrix(draft$pops)
glb <- draft$out$glb

colpts <-  c(1,2   ,1,1,1,1,1,1,1,1,1,1,1,1,1,1,4,1,1,1,1,3)
cexpts <-  c(1,1.75,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1.5,1,1,1,1,1.75)

plotprep(width=10,height=9)
parset()
pairs(pops[1:glb$SAUpop[1],c(2,3,4,6,7,10,13,19,24,25)],gap=0.2,pch=16,
      cex=cexpts,col=colpts)


pickcol <- c(1,2,3,4,6,7,19,9,10,13,17,24,25,26)

pops5 <- pops[1:glb$SAUpop[1],]
pops6 <- pops[(glb$SAUpop[1]+1):(glb$SAUpop[1]+glb$SAUpop[2]),]

pops5[order(pops5[,"msy"],decreasing=TRUE),pickcol]


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



# Fix sMaxDL to sDLMax----------------------




findinfluence <- function(ctrl,pops,reps,pickline,pickcol,calcpopC,
                          deleteyrs) {
  columns <- c("DLMax","L50","L95","SaMa","L50mat","Wta","Wtb","AvRec","SAU","msy",
               "B0")
  result <- matrix(0,nrow=reps,ncol=length(columns),dimnames=list(1:reps,columns))
  pickR <- match(columns,rownames(pops))
  whchline <- pickline
  whchcol <- pickcol
  dfile <- pathtopath(ctrl$rundir,ctrl$datafile)
  chgeloc <- c(3,4,6,8,9,3)
  wchline <- c(45,46,45,45,45,45)
  delvar <- matrix(0,nrow=6,ncol=reps)
  delvar[1,] <- runif(reps,min=0.1,max=0.5) # loc = 3 line 45 AvRec pop1
  delvar[2,] <- runif(reps,min=95,max=104) # loc = 4 L50mat
  delvar[3,] <- runif(reps,min=19.5,max=25) # loc = 6  DLMax
  delvar[4,] <- runif(reps,min=22,max=42)  # loc = 7   L95
  delvar[5,] <- runif(reps,min=3.1,max=3.3)   # loc = 8  Wtb
  delvar[6,] <- 1 - delvar[1,]              # loc = 3 line 46 AvRec pop2
  for (i in 1:reps) {  # i = 1
    indatr <- readLines(dfile)
    vals <- as.numeric(removeEmpty(unlist(strsplit(indatr[45],split=","))))
    for (varloc in 1:5) vals[chgeloc[varloc]] <- delvar[varloc,i]
    indatr[45] <- paste0(vals,collapse=" ,")
    vals <- as.numeric(removeEmpty(unlist(strsplit(indatr[46],split=","))))
    vals[chgeloc[6]] <- delvar[6,i]
    indatr[46] <- paste0(vals,collapse=" ,")
    writeLines(indatr,con=dfile)
    out <- do_condition(rundir=ctrl$rundir,controlfile=ctrl$controlfile,
                        calcpopC=calcpopC,
                        verbose=FALSE,doproduct=TRUE,dohistoric=TRUE,
                        matureL=c(70,200),wtatL=c(80,200),mincount=100,
                        uplimH=0.4,incH=0.005,deleteyrs=deleteyrs)
    pops <- t(out$pops)
    result[i,] <- pops[pickR,whchcol]
    print(round(pops[pickR,whchcol],3))
    print(paste0("replicate ",i),quote=FALSE)
  } # end of reps for loop
  return(invisible(result))
} # end of findinfluence



deleteyrs <- matrix(c(0,0,0,0,0,0,0,1993,1994,1995,2014,2015,2016,2020),
                    nrow=7,ncol=2,byrow=FALSE)
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
  library(qmdutils)
})
dropdir <- getDBdir()
prefixdir <- pathtopath(dropdir,"/A_codeUse/aMSETAS/NW/")
#prefixdir <- pathtopath(dropdir,"/A_codeUse/aMSETAS/scenarios/")
#prefixdir <- pathtopath(dropdir,"/A_codeUse/aMSEUse/SAHS/")

#postfixdir <- "saBC"
postfixdir <- "sau22"
verbose <- TRUE
rundir <- path.expand(filenametopath(prefixdir,postfixdir))
controlfile <- paste0("control",postfixdir,".csv")
outdir <- "C:/aMSE_scenarios/BC/"
confirmdir(rundir,ask=FALSE)
confirmdir(outdir,ask=FALSE)

load(file=paste0(outdir,"sau22",".RData"))


start <- Sys.time()
result <- findinfluence(ctrl=out$ctrl,pops=t(out$pops),reps=200,pickline=45,
                        pickcol=1,calcpopC=calcexpectpopC,deleteyrs=deleteyrs)

finish <- Sys.time()
cat("Time Taken = ",finish - start, "\n\n")

require(codeutils)
savedir <- "C:/Users/Malcolm/Dropbox/projects/AIRF/NWLML/"
filename <- pathtopath(savedir,"productivitydrivers.csv")
write.csv(result,file=filename)






chgevar <- "AvRec"

plotprep(width=10,height=5)
parset(plots=c(1,1))
model1 <- lm(result[,"msy"] ~ result[,chgevar])
plot(result[,chgevar],result[,"msy"],type="p",pch=16,cex=1.2,panel.first=grid)
abline(model1,lwd=2,col=2)
summary(model1)


plotprep(width=10,height=8)
parset(plots=c(1,1))
pairs(result[,c(1,3,5,7,8,10,11)])
#pairs(result[,c(1:4,6,7)])

head(result)  #+

model <- lm(result[,"msy"] ~ result[,"AvRec"]  + result[,"Wtb"] + result[,"DLMax"]
                             +  result[,"L95"] + result[,"L50mat"])
modout <- summary(model)
modout

model <- lm(result[,"msy"] ~ result[,"L50mat"])
modout <- summary(model)
modout



model <- lm(result[,"msy"] ~ result[,"AvRec"] + result[,"Wtb"]+ result[,"DLMax"]
             +  result[,"L95"]+ result[,"L50mat"])

modout <- summary(model)
modout
round(modout$coefficients,6)


)

cor(result[,"msy"],result[,"AvRec"])
cor(result[,"msy"],result[,"Wtb"])
cor(result[,"msy"],result[,"DLMax"])
cor(result[,"msy"],result[,"L95"])
cor(result[,"msy"],result[,"L50mat"])

require(codeutils)
rundir <- "C:/Users/Malcolm/Dropbox/projects/AIRF/NWLML/"
filename <- pathtopath(rundir,"productivitydrivers.csv")
#write.csv(result,file=filename)

result <- read.csv(filename,header=TRUE)




# population production

mod1 <- formula(result[,"msy"] ~ result[,"AvRec"])




summod3



printout(summod3)



apply(tabmat,1,mean)

apply(tabmat,1,mean) * 1.05
columns <- c("AvRec","Wtb","DLMax","L95","L50mat","msy")
rows <- c("BC","Wtb","DLMax","L95","L50mat")
result5 <- matrix(0,nrow=5,ncol=6,dimnames=list(rows,columns))
oldvals <- c(0,3.2,21.75,31.75,99.5)
newvals <- c(0,3.34879,22.8375,38.9375,104.475) # Wtb, DLMax,L95,L50mat
chgeloc <- c(0,9,6,8,4)
load(file=paste0("C:/aMSE_scenarios/BC/","sau22.RData"))
ctrl <- out$ctrl
dfile <- pathtopath(ctrl$rundir,ctrl$datafile)
deleteyrs <- matrix(c(0,0,0,0,0,0,0,1993,1994,1995,2014,2015,2016,2020),
                    nrow=7,ncol=2,byrow=FALSE)
require(TasHS)
out <- do_condition(rundir=ctrl$rundir,controlfile=ctrl$controlfile,
                    calcpopC=calcexpectpopC,
                    verbose=TRUE,doproduct=TRUE,dohistoric=TRUE,
                    matureL=c(70,200),wtatL=c(80,200),mincount=100,
                    uplimH=0.4,incH=0.005,deleteyrs=deleteyrs)
pops <- t(out$pops)
pickout <- match(c("AvRec","Wtb","DLMax","L95","L50mat","msy"),rownames(pops))
result5[1,] <- pops[pickout,1]
for (i in 2:5) { # i = 2
  indatr <- readLines(dfile)
  vals <- as.numeric(removeEmpty(unlist(strsplit(indatr[45],split=","))))
  vals[chgeloc[i]] <- newvals[i]
  indatr[45] <- paste0(vals,collapse=" ,")
  writeLines(indatr,con=dfile)
  out <- do_condition(rundir=ctrl$rundir,controlfile=ctrl$controlfile,
                      calcpopC=calcexpectpopC,
                      verbose=TRUE,doproduct=TRUE,dohistoric=TRUE,
                      matureL=c(70,200),wtatL=c(80,200),mincount=100,
                      uplimH=0.4,incH=0.005,deleteyrs=deleteyrs)
  pops <- t(out$pops)
  pickout <- match(c("AvRec","Wtb","DLMax","L95","L50mat","msy"),rownames(pops))
  result5[i,] <- pops[pickout,1]
  indatr <- readLines(dfile)
  vals <- as.numeric(removeEmpty(unlist(strsplit(indatr[45],split=","))))
  vals[chgeloc[i]] <- oldvals[i]
  indatr[45] <- paste0(vals,collapse=" ,")
  writeLines(indatr,con=dfile)
}
result5


result5[,6]/result5[1,6]


library(codeutils)

extractRcode(indir="C:/Users/Malcolm/Dropbox/A_Code/aMSEGuide",
             rmdfile="04_Using_aMSE.qmd",
             filename="usingaMSE.R")



# Zone phaseplot------------------------

zonemeds <- out$zonesummary$zonemeds
harv <- zonemeds[,"HarvestR"]
depl <- zonemeds[,"depl"]
console=FALSE
filen=""
targdepl=0.4
limdepl=0.2
limH=0.15   # F = M
yrs <- as.numeric(rownames(zonemeds))
nyrs <- length(yrs)
startyr <- yrs[1]
condy <- out$glb$hyrs


plotprep(width=7,height=6,filename=filen,cex=0.9,verbose=FALSE)
parset(cex.lab=1.25)

maxH <- getmax(harv)
maxdepl <- getmax(depl,mult=1.2)
plot(depl,harv,type="p",pch=16,xlim=c(0,maxdepl),ylim=c(0,maxH),xaxs="i",yaxs="i",
     xlab="Zone Scale Mature Biomass Depletion ",ylab="Harvest Rate")
phaseplotpolygons(maxy=maxH,maxx=maxdepl,targx=targdepl,limx=limdepl,limy=limH)
suppressWarnings(arrows(x0=depl[1:(nyrs-1)],y0=harv[1:(nyrs-1)],
                        x1=depl[2:nyrs],y1=harv[2:nyrs],lwd=2,
                        length=0.075,code=2))
# when H values repeatedly zero one gets zero length arrows and a warning
# this is not a problem so it is suppressed.
points(depl[c(1,condy,nyrs)],harv[c(1,condy,nyrs)],pch=16,cex=3,col=c(1,6,4))
legend("topright",legend=c(startyr,yrs[condy],yrs[nyrs]),
       col=c(1,6,4),lwd=5,cex=1.5,bty="n")
if ((!console) & (setpar)) {
  addplot(filen,rundir=rundir,category="zonescale",caption)
}




# plotproduction------------------------------

plotproductivity2 <- function(rundir,product,glb,hsargs) {
  # rundir=rundir; product=production;glb=glb; hsargs=hsargs
  maxtarg <- hsargs$maxtarg
  xval <- findmsy(product)
  numpop <- glb$numpop
  if (numpop <= 16) { # otherwise too much detail to see. Ok for SA not TAS
    # by pop  popn<= 16-----
    # Yield vs Spawning biomass-----
    filen <- filenametopath(rundir,"production_SpB.png")
    plotprod(product,xname="MatB",xlab="Spawning Biomass t",
             ylab="Production t",filename = filen,devoff=FALSE)
    caption <- paste0("The production curve relative to each population's ",
                      "spawning biomass. The vertical lines identify the ",
                      "Bmsy values.")
    addplot(filen,rundir=rundir,category="Production",caption)
    # Yield vs Annual Harvest Rate-----
    filen <- filenametopath(rundir,"production_AnnH.png")
    plotprod(product,xname="AnnH",xlab="Annual Harvest Rate",filename = filen,
             devoff=FALSE)
    caption <- paste0("The production curve relative to the Annual ",
                      "Harvest Rate applied to each population. The ",
                      "vertical lines identify the Hmsy values.")
    addplot(filen,rundir=rundir,category="Production",caption)
    # Yield vs depletion-----
    filen <- filenametopath(rundir,"production_Deplet.png")
    plotprod(product,xname="Deplet",xlab="Population Depletion Level",
             filename = filen,devoff=FALSE)
    for (pop in 1:numpop) abline(v=xval[pop,"Deplet"],lwd=2,col=pop)
    caption <- paste0("The production curve relative to the depletion ",
                      "level of each population. The vertical lines identify ",
                      "the Depletion level giving rise to the MSY.")
    addplot(filen,rundir=rundir,category="Production",caption)
    # plot of Yield vs population depletion but constrained to within
    # 0.2 and 0.35 levels, to illustrate nearly flat rpoduction curve
    # and more clearly identify the population depletion at MSY
    filen <- filenametopath(rundir,"production_Deplet_0.2_0.35.png")
    plotprod(product,xname="Deplet",xlab="Population Depletion Level",
             xlimit=c(0.2,0.35),filename = filen,devoff=FALSE)
    for (pop in 1:numpop) abline(v=xval[pop,"Deplet"],lwd=2,col=pop)
    caption <- paste0("The production curve relative to the depletion ",
                      "level of each population. Here the x-axis is ",
                      "shortened to clarify the flatness of the production ",
                      "curve about the MSY points.")
    addplot(filen,rundir=rundir,category="Production",caption)
  } # end of numpop <= 16 if statement
  # what CPUE at Bmsy would be exhibited at MSY
  # by SAU-----------------------
  filen <- filenametopath(rundir,"production_CPUE_at_Bmsy.png")
  nsau <- glb$nSAU
  npop <- glb$numpop
  nh <- dim(product)[1]
  label <- dimnames(product)
  harv <- as.numeric(label[[1]])
  sauindex <- glb$sauindex
  saunames <- glb$saunames
  wts <- matrix(0,nrow=nh,ncol=npop,dimnames=list(label[[1]],label[[3]]))
  sauyield <- matrix(0,nrow=nh,ncol=nsau,dimnames=list(label[[1]],saunames))
  saumatB <- matrix(0,nrow=nh,ncol=nsau,dimnames=list(label[[1]],saunames))
  sauexB <- saumatB
  saucpue <- saumatB
  # Now do sau production
  for (i in 1:nh) {
    saumatB[i,] <- tapply(product[i,"MatB",],sauindex,sum,na.rm=TRUE)
    sauexB[i,] <- tapply(product[i,"ExB",],sauindex,sum,na.rm=TRUE)
    sauyield[i,] <- tapply(product[i,"Catch",],sauindex,sum,na.rm=TRUE)
    wts[i,] <- product[i,"Catch",]/sauyield[i,sauindex]
    saucpue[i,] <- tapply((product[i,"RelCE",] * wts[i,]),sauindex,sum,na.rm=TRUE)
  }
  label <- glb$saunames
  # yield x cpue---------------
  rows <- c("B0","Bmsy","MSY","Dmsy","CEmsy","Hmsy","Bexmsy")
  sauprod <- matrix(0,nrow=length(rows),ncol=nsau,dimnames=list(rows,saunames))
  plotprep(width=8,height=7,newdev=FALSE,filename=filen,verbose=FALSE)
  parset(plots=pickbound(nsau),margin=c(0.3,0.3,0.05,0.05),outmargin=c(1,1,0,0))
  for (i in 1:nsau) { # i=1
    ymax <- getmax(sauyield[,i])
    if (!is.null(maxtarg)) {
      xmax <- getmax(c(saucpue[2:nh,i],maxtarg[i]))
      xmin <- getmin(c(saucpue[2:nh,i],maxtarg[i]))
    } else {
      xmax <- getmax(saucpue[2:nh,i])
      xmin <- getmin(saucpue[2:nh,i])
    }
    pick <- which.max(sauyield[,i])
    msyce <- saucpue[pick,i]
    plot(saucpue[2:nh,i],sauyield[2:nh,i],type="l",lwd=2,xlab="",ylab="",
         panel.first=grid(),ylim=c(0,ymax),yaxs="i",xlim=c(xmin,xmax))
    abline(v=msyce,col=2,lwd=2)
    if (!is.null(maxtarg)) abline(v=maxtarg[i],lwd=2,col=4)
    msylab <- paste0("MSY = ",round(sauyield[pick,i],1))
    text(0.7*max(saucpue[,i]),0.92*ymax,msylab,cex=1.2,pos=4)
    text(0.8*max(saucpue[,i]),0.75*ymax,label[i],cex=1.5,pos=4)
    text(1.1*msyce,0.15*ymax,round(msyce,2),cex=1.25,pos=4)
    sauprod[,i] <- c(saumatB[1,i],saumatB[pick,i],sauyield[pick,i],
                     (saumatB[pick,i]/saumatB[1,i]),saucpue[pick,i],
                     harv[pick],sauexB[pick,i])
  }
  mtext("CPUE at Bmsy (note different scales)",side=1,outer=TRUE,cex=1.0,
        line = -0.1)
  mtext("Equilibrium Yield  (note different scales)",side=2,outer=TRUE,cex=1.0,
        line=-0.2)
  caption <- paste0("The equilibrium production vs CPUE for each SAU. The red ",
                    "vertical lines and related number represent the expected",
                    "CPUE when biomass is at Bmsy. ")
  if (!is.null(maxtarg)) caption <- paste0(caption,"The blue vertical line is ",
                                           "the maxtarg for that SAU.")
  addplot(filen,rundir=rundir,category="Production",caption)
  # SAU prod plots-------------------
  if (nsau > 1) {
    yield <- rowSums(product[,"Catch",])
    spb <- rowSums(product[,"MatB",])
    expB <- rowSums(product[,"ExB",])
  } else {
    yield <- product[,"Catch",]
    spb <- product[,"MatB",]
    expB <- product[,"ExB",]
  }
  Ht <- harv
  depletMSY <- spb/spb[1]
  pickmsy <- which.max(yield)
  maxy <- getmax(yield)
  filen <- pathtopath(rundir,"production_SpB_Total.png")
  plotprep(width=7,height=8,newdev=FALSE,filename=filen,verbose=FALSE)
  parset(plots=c(4,2),cex=0.9)
  plot(spb,yield,type="l",lwd=2,col=1,xlab="Spawning Biomass t",
       ylab="Production t",panel.first = grid(),
       ylim=c(0,maxy),yaxs="i")
  abline(v=spb[pickmsy],col=2,lwd=2)

  plot(Ht,spb,type="l",lwd=2,xlab="Annual Harvest Rate",
       ylab="Spawning Biomass t",panel.first = grid(),
       ylim=c(0,getmax(spb)),yaxs="i")
  abline(h=spb[pickmsy],col=2,lwd=2)
  abline(v=Ht[pickmsy],col=2,lwd=2)

  plot(expB,yield,type="l",lwd=2,col=1,xlab="Exploitable Biomass t",
       ylab="Production t",panel.first = grid(),
       ylim=c(0,maxy),yaxs="i")
  abline(v=expB[pickmsy],col=2,lwd=2)

  plot(Ht,expB,type="l",lwd=2,xlab="Annual Harvest Rate",
       ylab="Exploitable Biomass t",panel.first = grid(),
       ylim=c(0,getmax(expB)),yaxs="i")
  abline(h=expB[pickmsy],col=2,lwd=2)
  abline(v=Ht[pickmsy],col=2,lwd=2)

  plot(Ht,yield,type="l",lwd=2,col=1,xlab="Annual Harvest Rate",
       ylab="Production t",panel.first = grid(),
       ylim=c(0,maxy),yaxs="i")
  abline(v=Ht[pickmsy],col=2,lwd=2)

  plot(spb,depletMSY,type="l",lwd=2,ylab="Total Depletion Level",
       xlab="Spawning Biomass t",panel.first = grid(),
       ylim=c(0,1.05),yaxs="i")
  abline(h=depletMSY[pickmsy],col=2,lwd=2)
  abline(v=spb[pickmsy],col=2,lwd=2)

  plot(depletMSY,yield,type="l",lwd=2,col=1,xlab="Total Depletion Level",
       ylab="Production t",panel.first = grid())
  abline(v=depletMSY[pickmsy],col=2,lwd=2)

  plot(Ht,depletMSY,type="l",lwd=2,col=1,xlab="Annual Harvest Rate",
       ylab="Total Depletion Level",panel.first = grid(),
       ylim=c(0,1.05),yaxs="i")
  abline(h=depletMSY[pickmsy],col=2,lwd=2)
  abline(v=Ht[pickmsy],col=2,lwd=2)


  caption <- paste0("The production curves for the zone. Also the ",
                    "relationships between spawning biomass depletion and ",
                    "harvest rate.")
  addplot(filen,rundir=rundir,category="Production",caption)
  zoneprod <- cbind(catch=yield,matB=spb,harv=Ht,deplsB=depletMSY,
                    expB=expB)
  # add zone production to sauprod
  zonemsy <- zoneprod[which.max(zoneprod[,"catch"]),]
  Pzone <- numeric(nrow(sauprod)); names(Pzone) <- rownames(sauprod)
  Pzone[c(1:4,6:7)] <- c(B0=zoneprod[1,"matB"],Bmsy=zonemsy["matB"],
                         MSY=zonemsy["catch"],Dmsy=zonemsy["expB"],
                         Hmsy=zonemsy["harv"],Bexmsy=zonemsy["expB"])
  wts <- sauprod["MSY",]/Pzone["MSY"]; wts <- wts/sum(wts)
  Pzone[5] <- sum(sauprod["CEmsy",]*wts)
  # Pzone[c(1:3,6:7)] <- rowSums(sauprod[c(1:3,7),],na.rm=TRUE)
  #
  #  Pzone[4:6] <- rowSums(sauprod[4:6,]*wts)
  sauprodzone <- cbind(sauprod,zone=Pzone)
  # add sauprod table to production tab
  filen <- "Production_by_sau.csv"
  caption <- paste0("Productivity properties: B0, Bmsy, MSY, Dmsy, Cemsy,",
                    " Hmsy, and Bexmsy for each sau and zone. CEmsy is the",
                    " predicted CPUE by SAU at MSY weighted by the sauMSY. ",
                    "They must,therefore be treated as approximate.")
  addtable(sauprodzone,filen,rundir=rundir,category="Production",caption)
  # add zone production table
  filen <- "production_across_zone.csv"
  caption <- paste0("Zone production: Catch, matureB, HarvestR, and Mature ",
                    "depletion level. MSY at H = ",Ht[pickmsy])
  addtable(zoneprod,filen,rundir=rundir,category="Production",caption)
  return(invisible(list(sauprod=sauprod,zoneprod=zoneprod,
                        sauprodzone=sauprodzone)))
} # end of plotproductivity



outdir <- "C:/aMSE_scenarios/EG/"
out <- NULL
load(file=pathtopath(outdir,"EGMRall.RData"))
glb <- out$glb
hsargs <- out$hsargs
rundir <- out$ctrl$rundir
product <- out$production

plotproductivity2(rundir,product,glb,hsargs)

