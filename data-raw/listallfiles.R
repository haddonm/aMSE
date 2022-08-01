
 # infile="C:/Users/User/Dropbox/A_Code/aMSE/R/inputfiles.R"; allfuns=allfilesort[,"function"]


library(codeutils)

if (dir.exists("c:/Users/User/DropBox")) {
  ddir <- "c:/Users/User/DropBox/A_code/"
} else {
  ddir <- "c:/Users/Malcolm/DropBox/A_code/"
}

indir <- paste0(ddir,"aMSE/R/")
files <- c("aMSE_utils.R","defineZone.R","docondition.R","domse.R","dynamics.R",
           "generate_results.R","getfunctions.R","inputfiles.R","plotfuns.R",
           "projection.R")
outfile <- paste0(ddir,"aMSE/data-raw/allfile.csv")



x <- describefunctions(indir=indir,files=files,outfile=outfile)










infun <- "makeequilzone"

toplevel <- extractpathway(indir,infun,allfuns=x)
subdiv <- length(toplevel)
out <- vector(mode="list",length=subdiv)
out[[1]] <- toplevel
for (i in 2:subdiv) {
  out[[i]] <- extractpathway(indir,infun=toplevel[i],allfuns=x)
}
out[[1]] <- infun

filen = "C:/Users/Malcolm/Dropbox/A_Code/rutilsMH/data-raw/test.Rmd"
filen=""
setuprmd(filen=filen)
nfun <- length(out)
cat(' \n\n','# ',out[[1]],'  \n\n',sep="",file=filen, append=TRUE)
if (nfun > 1) {
  for (i in 2:nfun) { #  i = 5
    vectf <- out[[i]]
    nsub <- length(vectf)
    cat('## ',vectf[1],'  \n\n',sep="",file=filen, append=TRUE)
    if (nsub > 1) {
      for (j in 2:nsub) { # j = 2
        cat('### ',vectf[j],' \n\n',sep="",file=filen,append=TRUE)
        pick <- which(x[,"functions"] == vectf[j])
        if (length(pick) > 0) {
          refs <- removeEmpty(unlist(strsplit(x[pick,"references"],",")))
          nref <- length(refs)
          if (nref > 0) {
            for (k in 1:nref)
              cat('#### ',refs[k],'     \n',sep="",file=filen, append=TRUE)
            cat(' \n\n',sep="",file=filen, append=TRUE)
          }
        } else {
          stop("Unknown function found in allfuns within extractpathway \n")
        }
      }
    } # end of nsub test
  }
} # end of nfun test









# develop a working example-----------------------------------------------------
# this works down to x <- describefunctions(...)
filen <- tempfile("test",fileext=".R")
txt <- c("# this is a comment",
         "@title ...",
         "dummy <- function() {",
         "  out <- anotherdummy()",
         "  return(out)",
         "}",
         "# a possibly confusing use of function",
         "@title ...",
         "anotherdummy <- function() {",
         "  return(NULL)",
         "}")
write(txt,file=filen)
usedir <- paste0(tempdir(),"//")
filename <- tail(unlist(strsplit(filen,"\\",fixed=TRUE)),1)
x <- describefunctions(indir=usedir,files=filename,outfile="")

extractpathway(indir=usedir,"dummy",x)




