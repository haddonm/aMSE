


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
})
dropdir <- getDBdir()
prefixdir <- pathtopath(dropdir,"/A_codeUse/aMSETAS/hsargs/")
#prefixdir <- pathtopath(dropdir,"/A_codeUse/aMSEUse/")


postfixdir <- "BCtime"
verbose <- TRUE
rundir <- path.expand(filenametopath(prefixdir,postfixdir))
controlfile <- paste0("control",postfixdir,".csv")
outdir <- "C:/aMSE_scenarios/BC/"
confirmdir(rundir,ask=FALSE)
confirmdir(outdir,ask=FALSE)

load(file=paste0(outdir,"BCtime",".RData"))

prody <- out$production

matB <- prody[,"MatB",37:48]
catch <- prody[,"Catch",37:48]
totY <- rowSums(catch)

plot1(H,totY)

H <- as.numeric(rownames(catch))
for (i in 1:12){
  pick <- which(catch[,i] == max(catch[,i]))
  cat(i, pick,catch[pick,i],H[pick],"\n")


}





































