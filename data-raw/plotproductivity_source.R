# All these plots only use the product array
# characterize productivity ------------------------------------------
# plot of Yield vs Spawning biomass
xval <- findmsy(product)
numpop <- glb$numpop

file <- paste0("production_SpB_",ctrl$runlabel,".png")
filename <- filenametopath(plotdir,file)  #  filename=""
plotprep(width=7,height=4,newdev=FALSE,filename=filename,cex=0.9,
         verbose=FALSE)
plotprod(product,xname="MatB",xlab="Spawning Biomass t",
         ylab="Production t")
for (pop in 1:numpop) abline(v=xval[pop,"MatB"],lwd=2,col=pop)
legend("topright",paste0("P",1:numpop),lwd=3,col=c(1:numpop),bty="n",
       cex=1.2)
if (nchar(filename) > 0) dev.off()

caption <- "The production curve relative to each population's spawning biomass. The vertical lines identify the Bmsy values."
addfilename(filename,tabfile=plottabfile,"Production","plot",caption)

# time <- as.character(Sys.time())
#  cat(c(filename,caption,"Production",time," \n"),file=plottabfile,sep=",",append=TRUE)


# plot of Yield vs Annual Harvest Rate
file <- paste0("production_AnnH_",ctrl$runlabel,".png")
filename <- filenametopath(plotdir,file)  #  filename=""

plotprep(width=7,height=4,newdev=FALSE,filename=filename,cex=0.9,
         verbose=FALSE)
plotprod(product,xname="AnnH",xlab="Annual Harvest Rate")
for (pop in 1:numpop) abline(v=xval[pop,"AnnH"],lwd=2,col=pop)
legend("topright",paste0("P",1:numpop),lwd=3,col=c(1:numpop),bty="n",
       cex=1.2)
if (nchar(filename) > 0) dev.off()

caption <- "The production curve relative to the Annual Harvest Rate applied to each population. The vertical lines identify the Hmsy values."
addfilename(filename,tabfile=plottabfile,"Production","plot",caption)



# plot of Yield vs population depletion
file <- paste0("production_Deplet_",ctrl$runlabel,".png")
filename <- filenametopath(plotdir,file)
plotprep(width=7,height=4,newdev=FALSE,filename=filename,cex=0.9,
         verbose=FALSE)
xval <- findmsy(product)
numpop <- glb$numpop
plotprod(product,xname="Deplet",xlab="Population Depletion Level")
for (pop in 1:numpop) abline(v=xval[pop,"Deplet"],lwd=2,col=pop)
legend("topright",paste0("P",1:numpop),lwd=3,col=c(1:numpop),bty="n",
       cex=1.2)
if (nchar(filename) > 0) dev.off()

caption <- "The production curve relative to the depletion level of each population. The vertical lines identify the Depletion level giving rise to the MSY."
addfilename(filename,tabfile=plottabfile,"Production","plot",caption)


# plot of Yield vs population depletion but constrained to within
# 0.2 and 0.35 levels, to illustrate nearly flat rpoduction curve
# and more clearly identify the population depletion at MSY
file <- paste0("production_Deplet_0.2_0.35_",ctrl$runlabel,".png")
filename <- filenametopath(plotdir,file)
plotprep(width=7,height=4,newdev=FALSE,filename=filename,cex=0.9,
         verbose=FALSE)
xval <- findmsy(product)
numpop <- glb$numpop
plotprod(product,xname="Deplet",xlab="Population Depletion Level",
         xlimit=c(0.2,0.35))
for (pop in 1:numpop) abline(v=xval[pop,"Deplet"],lwd=2,col=pop)
legend("topright",paste0("P",1:numpop),lwd=3,col=c(1:numpop),bty="n",
       cex=1.2)
if (nchar(filename) > 0) dev.off()

caption <- "The production curve relative to the depletion level of each population. Here the x-axis is shortened to clarify the flatness of the production curve about the MSY points."
addfilename(filename,tabfile=plottabfile,"Production","plot",caption)

# Now do total production --------------------------------------------
yield <- rowSums(product[,"Catch",])
spb <- rowSums(product[,"MatB",])
Ht <- yield/spb
depletMSY <- spb/spb[1]
pickmsy <- which.max(yield)
maxy <- getmax(yield)

file <- paste0("production_SpB_Total_",ctrl$runlabel,".png")
filename <- filenametopath(plotdir,file)  #  filename=""
plotprep(width=7,height=6,newdev=FALSE,filename=filename,verbose=FALSE)
parset(plots=c(3,2),cex=0.9)
plot(spb,yield,type="l",lwd=2,col=1,xlab="Spawning Biomass t",
     ylab="Production t",panel.first = grid(),
     ylim=c(0,maxy),yaxs="i")
abline(v=spb[pickmsy],col=2,lwd=2)

plot(Ht,spb,type="l",lwd=2,xlab="Annual Harvest Rate",
     ylab="Spawning Biomass t",panel.first = grid(),
     ylim=c(0,getmax(spb)),yaxs="i")
abline(h=spb[pickmsy],col=2,lwd=2)
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

if (nchar(filename) > 0) dev.off()

caption <- "The production curves for the region. Also the relationships between spawning biomass depletion and harvest rate."
addfilename(filename,tabfile=plottabfile,"Production","plot",caption)


















