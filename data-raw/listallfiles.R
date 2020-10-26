
 # infile="C:/Users/User/Dropbox/A_Code/aMSE/R/inputfiles.R"; allfuns=allfilesort[,"function"]

findfuns <- function(infile,allfuns) {
  numfun <- length(allfuns)
  content <- readLines(con=infile)
  rfun <- tail(unlist(strsplit(infile,"/")),1)
  rfile <- substr(rfun,1,nchar(rfun)-2)
  funLines <- grep("function",content)
  titles <- grep("@title",content)
  testhash <- substr(content[funLines],1,4)
  omit <- grep("#",testhash)
  if (length(omit) > 0) {
    funLines <- funLines[-omit]
    testhash <- testhash[-omit]
  }
  omit2 <- grep("  ",testhash) # remove functions internal to other functions
  if (length(omit2) > 0) funLines <- funLines[-omit2]
  nfun <- length(funLines)
    outf <- as.data.frame(matrix("",nrow=nfun,ncol=1))
    bounds <- matrix(0,nrow=nfun,ncol=2,
                     dimnames=list(paste0(rfile,1:nfun),c("start","end")))
    bounds[,1] <- funLines + 1
    if (nfun > 1) {
      bounds[,2] <- c((titles[2:nfun] - 2),length(content))
  } else {
    bounds[,2] <- length(content)
  }
  for (i in 1:nfun) { # i=1
    funname <- removeEmpty(unlist(strsplit(content[funLines[i]],"<-"))[1])
    funcont <- content[bounds[i,1]:bounds[i,2]]
    testhash <- substr(funcont,1,5)
    omit <- grep("#",testhash)
    if (length(omit) > 0) funcont <- funcont[-omit]
    whichfun <- ", "
    for (j in 1:numfun) {  #  j = 65
      if (allfuns[j] != funname)
        if (length(grep(allfuns[j],funcont)) > 0)
          whichfun <- paste0(whichfun,allfuns[j],", ")
    }
    outf[i,] <- whichfun
  }
  return(outf)
} # end of findfuns


#' @title listfuns produces a listing of all functions in an input R file
#'
#' @description listfuns reads in a given R file and then identifies each
#'     function header within it and pulls out the function name, its syntax,
#'     the line-number in the file, and associates that with the filename.
#'
#' @param infile the R file to be examined
#'
#' @return a data.frame of syntax, function name, line number, and file name
#' @export
#'
#' @examples
#' print("wait for an example")
listfuns <- function(infile) { # infile=filename; console=FALSE
  content <- readLines(con=infile)
  rfun <- tail(unlist(strsplit(infile,"/")),1)
  rfile <- substr(rfun,1,nchar(rfun)-2)
  funLines <- grep("function",content)
  testhash <- substr(content[funLines],1,4)
  omit <- grep("#",testhash)
  if (length(omit) > 0) {
    funLines <- funLines[-omit]
    testhash <- testhash[-omit]
  }
  omit2 <- grep("  ",testhash) # remove functions internal to other functions
  if (length(omit2) > 0) funLines <- funLines[-omit2]
  nLine <- length(funLines)
  delF <- NULL
  for (i in 1:nLine) {
    tmpLine <- gsub(" ","",content[funLines[i]])
    if ((length(grep("function\\(",tmpLine)) == 0) |
        (substr(tmpLine,1,2) == "#'") |
        (length(grep("<-function",tmpLine)) == 0) |
        (length(grep("} #",tmpLine)) > 0)) delF <- c(delF,i)
  }
  ndelF <- length(delF)
  if (ndelF > 0) {
    funLines <- funLines[-delF]
  }
  if (ndelF == nLine) {
    txt <- paste0(infile,"  contained no recognizable functions")
    warning(cat(txt,"\n\n"))
    return(txt)
  } else {
    outlines <- sort(c(funLines))
    out <- content[outlines]
    funnames <- out
    n <- length(out)
    for (i in 1:n) {  # i=1
      out[i] <- gsub(" ","",(unlist(strsplit(out[i],"\\{")))[1])
      funnames[i] <- removeEmpty(unlist(strsplit(out[i],"<-"))[1])
      out[i] <- gsub("<-function","",out[i])
    }
    columns <- c("syntax","function","linenumber","file","references")
    rows <- paste0(rfile,1:n)
    outfuns <- as.data.frame(matrix(NA,nrow=n,ncol=length(columns),
                                    dimnames=list(rows,columns)))
    outfuns[,"syntax"] <- out
    outfuns[,"function"] <- funnames
    outfuns[,"linenumber"] <- funLines
    outfuns[,"file"] <- rfile
  }
  outfuns
  return(outfuns)
} # end of listfuns

indir <- "C:/Users/User/Dropbox/A_Code/aMSE/R/"
files <- c("aMSE_utils.R","defineZone.R","dynamics.R","generate_results.R",
           "getfunctions.R","HS-HCR.R","inputfiles.R","plotfuns.R")
nfiles <- length(files)
allfiles <- NULL

for (i in 1:nfiles)
  allfiles <- rbind(allfiles,listfuns(paste0(indir,files[i])))

allfilesort <- allfiles[order(allfiles[,"function"]),]

allrefs <- NULL
for (i in 1:nfiles)
  allrefs <- rbind(allrefs,findfuns(paste0(indir,files[i]),allfilesort[,"function"]))

allfiles[,"references"] <- allrefs
allfilesort <- allfiles[order(allfiles[,"function"]),]



outfile <- "C:/Users/User/Dropbox/A_Code/aMSE/data-raw/allfile.csv"
write.csv(allfilesort,file = outfile)
