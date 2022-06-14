# The objective here is to explore the description of growth as used in the MSE model,
# How to express summer and winter growth summer = Sep - Feb, winter = Mar - Aug
# include the option of different proportions of that growth by season


library(rutilsMH)
library(AbMSE)

datadir <- "./../../abMSE/abMSE/data-raw/"

data(condDat)
str(condDat)

glb <- condDat$globals

nblock <- condDat$nblock
blkdef <- condDat$blkpop
numpop <- condDat$numpop

midpts <- condDat$midpts
set.seed(condDat$randomseed)
blockNames <- condDat$blockNames
projectionLML <- condDat$projLML
historicalLML <- condDat$histLML
if (condDat$Condition)
  projLML <- historicalLML else projLML <- projectionLML


popdefs <- definepops(nblock,blockI,condDat) # define all pops in one go
popnum <- 1
popdef <- popdefs[popnum,]
pop <- makeabpop(popdef,midpts,projLML)

str(pop)

round(pop$Nt[,1],1)


Nclass <- glb$Nclass

growp <- popdef[1:4]
Nt0 <- numeric(length(midpts))
Nt0[1:2] <- c(1000,1000)



yrs <- 8
part <- c(1,2)
npart <- length(part)
Nt <- matrix(0,nrow=Nclass,ncol=(yrs*npart))
N1 <- matrix(0,nrow=Nclass,ncol=yrs,dimnames=list(midpts,1:yrs))
colnames(Nt) <- sort(c(paste0(1:yrs,part[1]),paste0(1:yrs,part[2])))
rownames(Nt) <- midpts
P1 <- 0.1
P2 <- 0.9
dvar <- 0.78

G1 <- STM(growp,midpts)
gpw <- c(growp[1]*P1,growp[2:3],growp[4]*dvar)
Gw <- STM(gpw,midpts)
gps<- c(growp[1]*P2,growp[2:3],growp[4]*dvar)
Gs <- STM(gps,midpts)

N1[,1] <- G1 %*% Nt0
Nt[,1] <- Gw %*% Nt0
for (yr in 2:yrs) N1[,yr] <- G1 %*% N1[,(yr-1)]
pyr <- 1
for (yrp in 2:(yrs*npart)) {
  if (!(yrp %% 2))  {
     Nt[,yrp] <- Gs %*% Nt[,(yrp-1)]
  } else {
     Nt[,yrp] <- Gw %*% Nt[,(yrp-1)]
  }
}

summer <- c(1,seq(2,16,2))

plotprep(width=7,height=4)
plot1(midpts,N1[,1],limity=c(0,400))
for (yr in 2:yrs) lines(midpts,N1[,yr],lwd=2,col=1)
for (yrp in 2:(yrs*npart))
  if (!(yrp %% 2)) lines(midpts,Nt[,yrp],lwd=2,col=2,lty=2)
ssq <- 0.0
for (yr in 2:yrs) ssq <- ssq + sum((N1[,yr] - Nt[,summer[yr]])^2)
ssq

