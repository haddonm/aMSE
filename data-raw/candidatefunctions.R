

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
 #' #  rundir=rundir; inarr=cdivmsy[[1]];glb=glbc[[1]];scenes=scenes[1];
 #' #  filen="";label="Catch / MSY";maxy=0
plotsceneproj <- function(rundir,inarr,glb,scene,filen="",label="",
                         maxy=0,Q=90,hline=NA) {
  nsau <- glb$nSAU
  saunames <- glb$saunames
  outmed <- makelist(saunames)
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
        ylim=c(0,ymax),ylab=saunames[sau],xlab="")
   for (i in 1:reps) lines(yrs,dat[,i],lwd=1,col="grey")
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
  mtext(paste0(scene,"  ",label),side=2,line=-0.2,outer=TRUE,cex=1.1)
  if (nchar(filen) > 0)
   addplot(filen=filen,rundir=rundir,category="C_vs_MSY",caption)
  return(invisible(outmed))
} # end of plotsceneproj




























