# These plots use both regionC and regionD
# characterize productivity ------------------------------------------
# plot of Yield vs Spawning biomass
mids <- glb$midpts
numpop <- glb$numpop

# maturation uses regionC---------------------------------------------
maturity <- getlistvar(regionC,"Maturity")
rownames(maturity) <- mids

file <- paste0("maturity_v_Length_",ctrl$runlabel,".png")
filename <- filenametopath(resdir,file)  #  filename=""
plotprep(width=7,height=4,newdev=FALSE,filename=filename,cex=0.9,
         verbose=FALSE)
plot(mids,maturity[,1],type="l",lwd=2,xlab="Shell Length mm",
     ylab="Proportion Mature",panel.first=grid(),xlim=c(50,210))
for (pop in 2:numpop) lines(mids,maturity[,pop],lwd=2,col=pop)
legend("topright",paste0("P",1:numpop),lwd=3,col=c(1:numpop),bty="n",
       cex=1.2)
if (nchar(filename) > 0) dev.off()

caption <- "The maturity vs length for each population."
addfilename(filename,resfile=resfile,"Biology","plot",caption)

# weight-at-length uses regionC---------------------------------------
WtL <- getlistvar(regionC,"WtL")
rownames(WtL) <- mids

file <- paste0("Weight_at_Length_",ctrl$runlabel,".png")
filename <- filenametopath(resdir,file)  #  filename=""
plotprep(width=7,height=4,newdev=FALSE,filename=filename,cex=0.9,
         verbose=FALSE)
plot(mids,WtL[,1],type="l",lwd=2,xlab="Shell Length mm",
     ylab="Weight Kg",panel.first=grid(),xlim=c(110,210))
for (pop in 2:numpop) lines(mids,WtL[,pop],lwd=2,col=pop)
legend("topleft",paste0("P",1:numpop),lwd=3,col=c(1:numpop),bty="n",
       cex=1.2)
if (nchar(filename) > 0) dev.off()

caption <- "The weight-at-length for each population. The x-axis is constrained to encompass legal sizes."
addfilename(filename,resfile=resfile,"Biology","plot",caption)

# emergence uses regionC----------------------------------------------
emerg <- getlistvar(regionC,"Emergent")
rownames(emerg) <- mids

file <- paste0("Emergence_at_Length_",ctrl$runlabel,".png")
filename <- filenametopath(resdir,file)  #  filename=""
plotprep(width=7,height=4,newdev=FALSE,filename=filename,cex=0.9,
         verbose=FALSE)
plot(mids,emerg[,1],type="l",lwd=2,xlab="Shell Length mm",
     ylab="Weight Kg",panel.first=grid(),xlim=c(105,150))
for (pop in 2:numpop) lines(mids,emerg[,pop],lwd=2,col=pop)
legend("topleft",paste0("P",1:numpop),lwd=3,col=c(1:numpop),bty="n",
       cex=1.2)
if (nchar(filename) > 0) dev.off()

caption <- "The emergence-at-length for each population. The x-axis is constrained to emphasize differences."
addfilename(filename,resfile=resfile,"Biology","plot",caption)

# initial numbers-at-size uses regionD--------------------------------

mids <- glb$midpts
numpop <- glb$numpop
Nt <- regionD$Nt[,1,]/1000.0
Ntt <- rowSums(regionD$Nt[,1,])/1000.0

file <- paste0("Total_Initial_Numbers-at-Size_",ctrl$runlabel,".png")
filename <- filenametopath(resdir,file)  #  filename=""
plotprep(width=7,height=6,newdev=FALSE,filename=filename,cex=0.9,
         verbose=FALSE)
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
addfilename(filename,resfile=resfile,"Biology","plot",caption)


# Tabulate biological properties uses regionC-------------------------

getvar <- function(indexvar,regC=regionC) {
        invar <- getlistvar(regC,indexvar)
        return(c(invar,sum(invar)))
}
numpop <- glb$numpop
rows <- c("M","R0","B0","effB0","ExB0","effExB0","MSY","MSYDepl","bLML")
columns <- c(paste0("p",1:numpop),"region")
numrow <- length(rows)
numcol <- length(columns)
results <- matrix(0,nrow=numrow,ncol=numcol,dimnames=list(rows,columns))
results["B0",] <- getvar("B0")
M <- getlistvar(regionC,"Me")
wtr <- (results["B0",1:numpop]/results["B0",(numpop+1)])
results["M",] <- c(M,sum(M*wtr))
results["R0",] <- getvar("R0")
results["effB0",] <- getvar("effB0")
results["ExB0",] <- getvar("ExB0")
results["effExB0",] <- getvar("effExB0")
results["MSY",] <- getvar("MSY")
MSYD <- getlistvar(regionC,"MSYDepl")
results["MSYDepl",] <- c(MSYD,sum(MSYD * wtr))
bLML <- getlistvar(regionC,"bLML")
results["bLML",] <- c(bLML,sum(bLML * wtr))
res <- round(results,3)
filename <- filenametopath(resdir,"regionbiology.csv")
write.table(res,file = filename,sep=",")
#  use tmp <- read.csv(file=filename,header=TRUE,row.names=1)


caption <- paste("Population amd Regional Biological Properties.",
                         "Where the regional total is an average it is weighted",
                         "relative to the proportion of total B0.",collapse=" ")
addfilename(filename,resfile=resfile,"Tables","table",caption)


# total regional productivity ----------------------------------------

# numrow <- length(product[,1,1])
# numvar <- length(product[1,,1])
# columns <- names(product[1,,1])
# rows <- names(product[,1,1])
# result <- matrix(0,nrow=numrow,ncol=numvar,dimnames=list(rows,columns))
# cols <- c(1,2,4)
# for (i in cols) result[,i] <- rowSums(product[,i,])
# result[,"AnnH"] <- result
# result




























