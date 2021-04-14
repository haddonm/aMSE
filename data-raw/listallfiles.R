
 # infile="C:/Users/User/Dropbox/A_Code/aMSE/R/inputfiles.R"; allfuns=allfilesort[,"function"]









if (dir.exists("c:/Users/User/DropBox")) {
  ddir <- "c:/Users/User/DropBox/A_code/"
} else {
  ddir <- "c:/Users/Malcolm/DropBox/A_code/"
}

indir <- paste0(ddir,"aMSE/R/")
files <- c("aMSE_utils.R","defineZone.R","dynamics.R","generate_results.R",
           "getfunctions.R","HS-HCR.R","inputfiles.R","plotfuns.R",
           "projection.R","RcppExports.R")
outfile <- paste0(ddir,"aMSE/data-raw/allfile.csv")

x <- describefunctions(indir,files,outfile=outfile)




