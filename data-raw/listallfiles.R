
listfuns <- function(infile,console=TRUE) { # infile=filename; console=FALSE
  content <- readLines(con=infile)
  rfun <- tail(unlist(strsplit(infile,"/")),1)
  rfile <- substr(rfun,1,nchar(rfun)-2)
  funLines <- grep("function",content)
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
      if (console) {
        if (nchar(out[i]) == 2) {
          cat("\n")
        } else {
          cat(outlines[i],out[i],"\n")
        }
      } # end of if(console)
    }
    columns <- c("syntax","function","linenumber","file")
    rows <- paste0(rfile,1:n)
    outfuns <- as.data.frame(matrix(NA,nrow=n,ncol=4,
                                    dimnames=list(rows,columns)))
    outfuns[,"syntax"] <- out
    outfuns[,"function"] <- funnames
    outfuns[,"linenumber"] <- funLines
    outfuns[,"file"] <- rfile
  }
  outfuns
  return(invisible(outfuns))
}



first <- listfuns("C:/Users/User/Dropbox/rcode2/aMSE/R/inputfiles.R",console=FALSE)
second <- listfuns("C:/Users/User/Dropbox/rcode2/aMSE/R/defineZone.R",console=FALSE)
third <- listfuns("C:/Users/User/Dropbox/rcode2/aMSE/R/aMSE_utils.R",console=FALSE)
fourth <- listfuns("C:/Users/User/Dropbox/rcode2/aMSE/R/makehtml_funs.R",console=FALSE)
fifth <- listfuns("C:/Users/User/Dropbox/rcode2/aMSE/R/plotfuns.R",console=FALSE)
sixth <- listfuns("C:/Users/User/Dropbox/rcode2/aMSE/R/getfunctions.R",console=FALSE)

allfiles <- rbind(first,second,third,fourth,fifth,sixth)

allfilesort <- allfiles[order(allfiles[,"function"]),]


outfile <- "C:/Users/User/Dropbox/rcode2/aMSE/data-raw/allfile.csv"
write.csv(allfilesort,file = outfile)
