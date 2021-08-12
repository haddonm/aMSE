



library(rutilsMH)

if (dir.exists("c:/Users/User/DropBox")) {
  ddir <- "c:/Users/User/DropBox/A_code/"
} else {
  ddir <- "c:/Users/Malcolm/DropBox/A_code/"
}

indir <- paste0(ddir,"rutilsMH/R/")

files <- dir(indir)

nfiles <- length(files)

total <- 0
comments <- 0
for (i in 1:nfiles) {
  dat <- readLines(paste0(indir,files[i]))
  lenchar <- nchar(dat)
  total <- total + countgtzero(lenchar)
  comments <- comments + length(grep("#\'",dat))
}
{cat("number of non-empty lines:  ",total,"\n")
cat("number of lines of comment: ",comments,"\n")
cat("number of lines of code:    ",total - comments,"\n")}
# dat <- readLines(paste0(indir,files[2]))
# lenchar <- nchar(dat)
# countgtzero(lenchar)
#
# grep("#'",dat)
