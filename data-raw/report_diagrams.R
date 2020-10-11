library(diagram)
library(rutilsMH)


# Initiation alternatives-------------------------------------------------------

plotprep(width=7, height=7,newdev=FALSE)
par(mar=c(1,1,1,1))
openplotmat()

locdat <- c(0.5,0.96,
            0.2,0.8,
            0.8,0.8,
            0.5,0.65,
            0.2,0.5,
            0.5,0.5,
            0.8,0.5,
            0.5,0.325,
            0.5,0.175,
            0.5,0.04)
locs <- matrix(locdat,nrow=10,ncol=2,byrow=TRUE)

elpos <- coordinates(locs,N=10)

textrect(elpos[1,], 0.4,0.04, lab = "Unfished Equilibrium",
         shadow.size = 0.0, cex = 1.5)
textrect(elpos[2,], 0.1, 0.04,lab = "No Fishing", box.col = "lightblue",
         shadow.size = 0.0, cex = 1.5)
textrect(elpos[3,], 0.2, 0.06,lab = c("Constant Harvest Rate","unknown catches"),
         box.col = "lightblue",shadow.size = 0.0, cex = 1.5)
textrect(elpos[4,], 0.4, 0.04, lab = c("Depletion before Conditioning on Catch"),
         shadow.size = 0.0, cex = 1.5)
textrect(elpos[5,], 0.1, 0.04, lab = c("No Fishing"), box.col = "lightblue",
         shadow.size = 0.0, cex = 1.5)
textrect(elpos[6,], 0.15, 0.04, lab = c("Set Harvest Rate"),box.col = "lightblue",
         shadow.size = 0.0, cex = 1.5)
textrect(elpos[7,], 0.125, 0.04, lab = c("Known Catches"),box.col = "lightblue",
         shadow.size = 0.0, cex = 1.5)
textrect(elpos[8,], 0.25, 0.04, lab = c("Operating Model Initiated"),
         shadow.size = 0.0, cex = 1.5)
textrect(elpos[9,], 0.15, 0.04, lab = c("Application of HCR"),
         shadow.size = 0.0, cex = 1.5)
textrect(elpos[10,], 0.125, 0.04, lab = c("Final State"),
         shadow.size = 0.0, cex = 1.5)

straightarrow(from=c(0.2,0.92),to=c(0.2,0.84),lwd=2,arr.pos=0.6)
straightarrow(from=c(0.2,0.76),to=c(0.2,0.69),lwd=2,arr.pos=0.6)
straightarrow(from=c(0.8,0.92),to=c(0.8,0.86),lwd=2,arr.pos=0.6)
straightarrow(from=c(0.8,0.74),to=c(0.8,0.69),lwd=2,arr.pos=0.6)
straightarrow(from=c(0.2,0.61),to=c(0.2,0.54),lwd=2,arr.pos=0.6)
straightarrow(from=c(0.5,0.61),to=c(0.5,0.54),lwd=2,arr.pos=0.6)
straightarrow(from=c(0.8,0.61),to=c(0.8,0.54),lwd=2,arr.pos=0.6)
lines(c(0.025,0.975),c(0.715,0.715),lwd=2,lty=3)
straightarrow(from=c(0.5,0.46),to=c(0.5,0.365),lwd=2,arr.pos=0.6)
bentarrow(from=c(0.2,0.46),to=c(0.25,0.325),lwd=2,path="V")
bentarrow(from=c(0.8,0.46),to=c(0.75,0.325),lwd=2,path="V")
lines(c(0.025,0.975),c(0.415,0.415),lwd=2,lty=3)
straightarrow(from=c(0.5,0.285),to=c(0.5,0.215),lwd=2,arr.pos=0.6)
straightarrow(from=c(0.5,0.135),to=c(0.5,0.08),lwd=2,arr.pos=0.6)


text("Projection",x=0.025,y=0.175,srt="90",cex=1.5)
text("Conditioning",x=0.025,y=0.56,srt="90",cex=1.5)
text("Pre-Conditioning",x=0.025,y=0.86,srt="90",cex=1.5)

# End of Initiation alternatives------------------------------------------------







