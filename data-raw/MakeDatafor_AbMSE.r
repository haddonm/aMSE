#rm(list=ls())   # Cleanup the R console if required
# Set up the run ----------


# library(r4cpue)
library(mhutils)
library(aMSE)
library(microbenchmark)


options("show.signif.stars"=FALSE,"stringsAsFactors"=FALSE,
        "max.print"=50000,"width"=240)
#     listFunctions("C:/A_Mal/A_Book/rcode/agestruct/age_utils.r")

# two block abalone zone with 6 populations -----------


# Setup data for aMSE -----------------------------------------------------

datadir <- "./../../rcode2/aMSE/data-raw/"

infile <- paste0(datadir,"twoblock.csv")

#datafile <- datafileTemplate(numblock=3,filename=outfile)
condDat <- readdataFile(infile)
str(condDat,max.level = 2)

save(condDat,file=paste0(datadir,"condDat.RData"))


# check and transfer -------------------------------------------------------



tools::checkRdaFiles(paths=datadir)
tools::resaveRdaFiles(paths=datadir,compress="auto")
tools::checkRdaFiles(paths=datadir)

devtools::use_data(condDat,
                   pkg="C:/A_mal/Rcode/Abalone/AbMSE",
                   internal=FALSE, overwrite=TRUE)

devtools::use_data(condDat3,
                   pkg="C:/A_mal/Rcode/Abalone/AbMSE",
                   internal=FALSE, overwrite=TRUE)

devtools::use_data(abdat,
                   pkg="C:/A_mal/Rcode/Abalone/abspatial",
                   internal=FALSE, overwrite=TRUE)

















