

library(MQMF)

data(abdat)
# or read.csv("filename_and_path",header=TRUE)
# format:
# year catch cpue
# 1978  123   97
# 1979  234   98
# ..............
# By all means use abdat to start with but this will be better
# done with your own data. set it up in Excel and save as a
# csv file


# plot the data. A reasonable initial guess at K is something
#like 10x the maximum catch
plotspmdat(abdat)


#cross correlation between cpue and catch cf Fig 7.2
# is there a negative correlation between cpue and catch?
parset()
ccf(x=abdat[,"catch"],y=abdat[,"cpue"],type="correlation",
    ylab="Correlation",plot=TRUE,main="")

cor(abdat[,"catch"],abdat[,"cpue"])

# That's fortunate! Now try an initial fit with guessed
# parameters. 0.25 is clearly too low, try a bigger number

param <- log(c(r=0.3,K=8000,Binit=3000,sigma=0.5))
negatL <- negLL(param,simpspm,abdat,logobs=log(schaef[,"cpue"]),schaefer=FALSE)
negatL
ans <- plotspmmod(inp=param,indat=abdat,schaefer=FALSE,
                  addrmse=TRUE,plotprod=FALSE)


# Now you have an approximate solution you can fit the thing
#
param <- log(c(r=0.3,K=8000,Binit=3000,sigma=0.5))
pnams <- c("r","K","Binit","sigma")
best <- optim(par=param,fn=negLL,funk=simpspm,indat=abdat,
              logobs=log(abdat[,"cpue"]),method="BFGS",schaefer=FALSE)
outfit(best,digits=4,title="Optim",parnames = pnams)
cat("\n")
best2 <- nlm(negLL,best$par,funk=simpspm,indat=abdat,
             logobs=log(abdat[,"cpue"]),schaefer=FALSE)
outfit(best2,digits=4,title="nlm",parnames = pnams)

# Now plot the optimum
ans <- plotspmmod(inp=best2$estimate,indat=abdat,schaefer=FALSE,
                  addrmse=TRUE,plotprod=TRUE)


# Now examine the uncertainties first described in
# Chapter 6 and specifically in chapter 7

data(abdat); fish <- as.matrix(abdat)
param <- log(c(r=0.3,K=11500,Binit=3300,sigma=0.05))
bestmod <- nlm(f=negLL1,p=param,funk=simpspm,schaefer=FALSE,
               logobs=log(fish[,"cpue"]),indat=fish,hessian=TRUE)
optpar <- exp(bestmod$estimate)
ans <- plotspmmod(inp=bestmod$estimate,indat=fish,schaefer=FALSE,
                  target=0.4,addrmse=TRUE, plotprod=FALSE)


# R-chunk 57 Page 312

out <- spm(bestmod$estimate,indat=fish,schaefer=FALSE)
str(out, width=65, strict.width="cut")

library(knitr)


catches <- seq(700,1000,50)   # projyr=10 is the default
projans <- spmprojDet(spmobj=out,projcatch=catches,plotout=TRUE)

### Accounting for Uncertainty
### Using Asymptotic Errors
# R-chunk 60  Page 315
# generate parameter vectors from a multivariate normal
# project dynamics under a constant catch of 900t

library(mvtnorm)
matpar <- parasympt(bestmod,N=1000) #generate parameter vectors
projs <- spmproj(matpar,fish,projyr=10,constC=900)#do dynamics

# R-chunk 61  Page 315
# Fig 7.30  1000 replicate projections asymptotic errors

outp <- plotproj(projs,out,qprob=c(0.1,0.5),refpts=c(0.2,0.4))

### Using Bootstrap Parameter Vectors
# R-chunk 62  Page 316
#bootstrap generation of plausible parameter vectors for Fox

reps <- 1000
boots <- spmboot(bestmod$estimate,fishery=fish,iter=reps,schaefer=FALSE)
matparb <- boots$bootpar[,1:4] #examine using head(matparb,20)

# R-chunk 63  Page 316
#bootstrap projections. Lower case b for boostrap  Fig7.31

projb <- spmproj(matparb,fish,projyr=10,constC=900)
outb <- plotproj(projb,out,qprob=c(0.1,0.5),refpts=c(0.2,0.4))

### Using Samples from a Bayesian Posterior
# R-chunk 64  Pages 317 - 318
#Generate 1000 parameter vectors from Bayesian posterior

param <- log(c(r=0.3,K=11500,Binit=3300,sigma=0.05))
set.seed(444608)
N <- 1000
result <- do_MCMC(chains=1,burnin=100,N=N,thinstep=2048,
                  inpar=param,infunk=negLL,calcpred=simpspmC,
                  calcdat=as.matrix(fish),obsdat=log(fish[,"cpue"]),
                  priorcalc=calcprior,schaefer=FALSE,
                  scales=c(0.065,0.055,0.1,0.475))
parB <- result[[1]][[1]] #capital B for Bayesian
cat("Acceptance Rate = ",result[[2]],"\n")

# R-chunk 65  Page 318
# auto-correlation, or lack of, and the K trace Fig 7.32

oldp <- parset(plots=c(2,1),cex=0.85)
acf(parB[,2],lwd=2)
plot(1:N,parB[,2],type="l",ylab="K",ylim=c(8000,19000),xlab="")
par(oldp)  # return par to old settings; this line not in book

# R-chunk 66  Page 318
#  Fig 7.33

matparB <- as.matrix(parB[,1:4]) # B for Bayesian
projs <- spmproj(matparB,fish,constC=900,projyr=10) # project them
plotproj(projs,out,qprob=c(0.1,0.5),refpts=c(0.2,0.4)) #projections

## Concluding Remarks
## Appendix: The Use of Rcpp to Replace simpspm
# R-chunk 67  Page 321

library(Rcpp)
cppFunction('NumericVector simpspmC(NumericVector pars,
             NumericMatrix indat, LogicalVector schaefer) {
   int nyrs = indat.nrow();
   NumericVector predce(nyrs);
   NumericVector biom(nyrs+1);
   double Bt, qval;
   double sumq = 0.0;
   double p = 0.00000001;
   if (schaefer(0) == TRUE) {
     p = 1.0;
   }
   NumericVector ep = exp(pars);
   biom[0] = ep[2];
   for (int i = 0; i < nyrs; i++) {
      Bt = biom[i];
      biom[(i+1)] = Bt + (ep[0]/p)*Bt*(1 - pow((Bt/ep[1]),p)) -
                          indat(i,1);
      if (biom[(i+1)] < 40.0) biom[(i+1)] = 40.0;
      sumq += log(indat(i,2)/biom[i]);
    }
    qval = exp(sumq/nyrs);
    for (int i = 0; i < nyrs; i++) {
      predce[i] = log(biom[i] * qval);
    }
    return predce;
 }')




library(microbenchmark)


param <- log(c(r=0.3,K=11500,Binit=3300,sigma=0.05))

microbenchmark(
  x <- simpspm(param,as.matrix(fish),schaefer=FALSE),
  xC <- simpspmC(param,as.matrix(fish),schaefer=FALSE)
)







