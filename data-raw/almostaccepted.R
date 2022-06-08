



options("show.signif.stars"=FALSE,
        "stringsAsFactors"=FALSE,
        "max.print"=50000,
        "width"=240)
suppressPackageStartupMessages({
  library(aMSE)
  library(TasHS)
  library(hutils)
  library(hplot)
  library(makehtml)
  library(knitr)
})
dropdir <- getDBdir()
prefixdir <- paste0(dropdir,"A_codeUse/aMSEUse/scenarios/")

startime <- Sys.time()
postfixdir <- "M15h7L75"
verbose <- TRUE
#hsfile <- "TasHS1_Tas.R"
rundir <- filenametopath(prefixdir,postfixdir)
controlfile <- paste0("control",postfixdir,".csv")
outdir <- "C:/aMSE_scenarios/M15h7L75/"
confirmdir(outdir)


# explore zonebiology --------------------------------------

filen <- "zonebiology.csv"

filename <- filenametopath(rundir,filen)
filename

biol <- read.csv(filename,header=TRUE)

plotprep(width=8, height=4.5,newdev=FALSE)
parset()
x <- biol[,"B0"]*biol[,"MSYDepl"]
plot1(x,biol[,"MSY"],type="p",pch=16,cex=1,col=factor(biol[,"SAU"]))
model <- lm(biol[,"MSY"] ~ x - 1)
abline(model,col=1,lwd=2)
summary(model)






















