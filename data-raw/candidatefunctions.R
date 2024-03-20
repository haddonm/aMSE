

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






