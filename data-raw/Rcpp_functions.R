
library(Rcpp)

cppFunction('List oneyearcatC(List inpopC, NumericVector inNt, double Nclass,
            double incat, int yr) {
    int nyrs = indat.nrow();
    NumericVector Ne(Nclass);
    NumericVector Cat(Nclass);
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
      biom[(i+1)] = Bt + (ep[0]/p)*Bt*(1 - pow((Bt/ep[1]),p)) - indat(i,1);
      if (biom[(i+1)] < 40.0) biom[(i+1)] = 40.0;
      sumq += log(indat(i,2)/biom[i]);
    }
    qval = exp(sumq/nyrs);
    for (int i = 0; i < nyrs; i++) {
      predce[i] = log(biom[i] * qval);
    }

    return List::create(vect, NaL, catchN, NumNe);
}')


oneyearcat <- function(inpopC,inNt,Nclass,incat,yr) {  #
   # yr=2; pop=2; inpopC=zoneC[[pop]]; inNt=zoneD$Nt[,yr-1,pop];
   # Nclass=glb$Nclass; inH=0.05;
   MatWt <- inpopC$MatWt/1e06
   SelectWt <- inpopC$SelWt[,yr]/1e06
   selyr <- inpopC$Select[,yr]
   Ne <- numeric(Nclass)
   Cat <- numeric(Nclass)
   Os <- exp(-inpopC$Me/2)
   NumNe <- (Os * (inpopC$G %*% inNt))
   midyexpB <- sum(SelectWt * NumNe) #SelectWt=Select*WtL =midyrexB
   estH <- min(incat/midyexpB,0.8) # no more than 0.8 harvest rate
   Fish <- 1-(estH*selyr)
   newNt <- (Os * (Fish * NumNe))
   Cat <- (estH*selyr) * NumNe  #numbers at size in the catch
   ExploitB <- sum(SelectWt * newNt) # end of year exploitable biomass
   avExpB <- (midyexpB + ExploitB)/2.0 #av start and end exploitB
   MatureB <- sum(MatWt*newNt)
   Catch <- sum(inpopC$WtL*Cat)/1e06
   ce <- inpopC$popq * avExpB * 1000.0  #ExploitB
   vect <- c(exploitb=ExploitB,midyexpB=midyexpB,matureb=MatureB,
             catch=Catch,cpue=ce)
   ans <- list(vect=vect,NaL=newNt,catchN=Cat,NumNe=NumNe)
   return(ans)
} # End of oneyearcat





