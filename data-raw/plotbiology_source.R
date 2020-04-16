
# characterize productivity ------------------------------------------
# plot of Yield vs Spawning biomass
mids <- glb$midpts
numpop <- glb$numpop

# maturation ---------------------------------------------------------
maturity <- getlistvar(regionC,"Maturity")
rownames(maturity) <- mids

file <- paste0("maturity_v_Length_",ctrl$runlabel,".png")
filename <- filenametopath(plotdir,file)  #  filename=""
plotprep(width=7,height=4,newdev=FALSE,filename=filename,cex=0.9,verbose=FALSE)
plot(mids,maturity[,1],type="l",lwd=2,xlab="Shell Length mm",
     ylab="Proportion Mature",panel.first=grid(),xlim=c(50,210))
for (pop in 2:numpop) lines(mids,maturity[,pop],lwd=2,col=pop)
legend("topright",paste0("P",1:numpop),lwd=3,col=c(1:numpop),bty="n",
       cex=1.2)
if (nchar(filename) > 0) dev.off()

caption <- "The maturity vs length for each population."
addfilename(filename,caption,"Biology")

# weight-at-length ---------------------------------------------------
WtL <- getlistvar(regionC,"WtL")
rownames(WtL) <- mids

file <- paste0("Weight_at_Length_",ctrl$runlabel,".png")
filename <- filenametopath(plotdir,file)  #  filename=""
plotprep(width=7,height=4,newdev=FALSE,filename=filename,cex=0.9,verbose=FALSE)
plot(mids,WtL[,1],type="l",lwd=2,xlab="Shell Length mm",
     ylab="Weight Kg",panel.first=grid(),xlim=c(110,210))
for (pop in 2:numpop) lines(mids,WtL[,pop],lwd=2,col=pop)
legend("topleft",paste0("P",1:numpop),lwd=3,col=c(1:numpop),bty="n",
       cex=1.2)
if (nchar(filename) > 0) dev.off()

caption <- "The weight-at-length for each population. The x-axis is constrained to encompass legal sizes."
addfilename(filename,caption,"Biology")

# emergence  ---------------------------------------------------------
emerg <- getlistvar(regionC,"Emergent")
rownames(emerg) <- mids

file <- paste0("Emergence_at_Length_",ctrl$runlabel,".png")
filename <- filenametopath(plotdir,file)  #  filename=""
plotprep(width=7,height=4,newdev=FALSE,filename=filename,cex=0.9,verbose=FALSE)
plot(mids,emerg[,1],type="l",lwd=2,xlab="Shell Length mm",
     ylab="Weight Kg",panel.first=grid(),xlim=c(105,150))
for (pop in 2:numpop) lines(mids,emerg[,pop],lwd=2,col=pop)
legend("topleft",paste0("P",1:numpop),lwd=3,col=c(1:numpop),bty="n",
       cex=1.2)
if (nchar(filename) > 0) dev.off()

caption <- "The emergence-at-length for each population. The x-axis is constrained to emphasize differences."
addfilename(filename,caption,"Biology")

# initial numbers-at-size --------------------------------------------

mids <- glb$midpts
numpop <- glb$numpop
Nt <- regionD$Nt[,1,]/1000.0
Ntt <- rowSums(regionD$Nt[,1,])/1000.0

file <- paste0("Total_Initial_Numbers-at-Size_",ctrl$runlabel,".png")
filename <- filenametopath(plotdir,file)  #  filename=""
plotprep(width=7,height=6,newdev=FALSE,filename=filename,cex=0.9,verbose=FALSE)
parset(plots=c(2,1))
plot(mids[5:105],Ntt[5:105],type="l",lwd=2,xlab="Shell Length mm (5 - 210mm)",
     ylab="Numbers-at_size '000s",panel.first=grid())
maxy <- getmax(Nt[5:105,])
plot(mids[5:105],Nt[5:105,1],type="l",lwd=2,xlab="Shell Length mm (5 - 210mm)",
     ylab="Numbers-at_size '000s",panel.first=grid(),ylim=c(0,maxy))
for (pop in 2:numpop) lines(mids[5:105],Nt[5:105,pop],lwd=2,col=pop)
abline(h=0.0,col="darkgrey")
legend("topright",paste0("P",1:numpop),lwd=3,col=c(1:numpop),bty="n",
       cex=1.2)

if (nchar(filename) > 0) dev.off()

caption <- "The numbers-at-size for the whole region and for each population separately. The recruitment numbers are omitted for clarity."
addfilename(filename,caption,"Biology")











