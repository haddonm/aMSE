
#' @title addfilename adds a filename to the autoresult csv file
#'
#' @description addfilename is used to facilitate the produciton of
#'     the HTML results summary for a particular run. This depends
#'     upon a csv file containing the names of each file to be
#'     plotted or tabulated. This function adds a filename and the
#'     supporting caption and category, without one needing to
#'     remember the syntax.
#'
#' @param filename the full path and filename for the file being added
#' @param tabfile the file to be added to, which is defined by
#'     setuphtml found in aMSE_utils
#' @param category what HTML tab should it be displayed on?
#' @param type the type of addition, either a plot or table
#' @param caption the caption for the figure or table

#'
#' @return nothing but it does add a line to filen
#' @export
#'
#' @examples
#' indir <- tempdir()
#' plotdir <- filenametopath(indir,"plots")
#' dirExists(plotdir,verbose=FALSE)
#' plottabfile <- setuphtml(plotdir,"example_only")
#' filename <- filenametopath(plotdir,"example.png")
#' png(filename=filename,width=7,height=4,units="in",res=300)
#' plot(runif(100),runif(100),type="p")
#' graphics.off()  # could use dev.off()
#' addfilename(filename=filename,tabfile=plottabfile,"A_category",
#'             type="plot",caption="Example Figure")
#' dir(plotdir)
addfilename <- function(filename,tabfile,category,type,caption) {
  cat(c(filename,category,type,as.character(Sys.time()),caption," \n"),
      file=tabfile,sep=",",append=TRUE)
}

#' @title dirExists: Checks for the existence of a directory
#'
#' @description dirExists: does a directory exist? It uses dir.exists
#'     and reports existence if already present and uses dir.create
#'     it it does not exist, but avoids the warning message is one
#'     already exists. The option of not creating a new directory is
#'     also present.
#'
#' @param indir a character string containing the name of the directory
#'     whose existence is to be checked before it is created if it
#'     does not already exist.
#' @param make if the directory does NOT exist should it be created.
#'     default = TRUE; if make=FALSE and a directory does not exist
#'     a warning will be given to the console.
#' @param verbose default=TRUE, prints directory status to the console,
#'     If make is set to FALSE and a directory does not exist a
#'     warning will always be given.
#'
#' @return a message to the screen if the directory exists or is
#'     created; if make is TRUE then it also creates the directory as
#'     listed in 'indir'.
#' @export
#'
#' @examples
#' indirect <- getwd()
#' dirExists(indirect)
dirExists <- function(indir,make=TRUE,verbose=TRUE) {
  if (dir.exists(indir)) {
    if (verbose) cat(indir,":  exists  \n")
  } else {
    if (make) {
      dir.create(indir)
      if (verbose) cat(indir,":  created  \n")
    } else {
      warning(cat(indir,":  does not exist \n"))
    }
  }
}  # end of dirExists


#' @title filenametopath safely add a filename to a path
#'
#' @description filenametopath add a filename to a path safely, using
#'     pathtype to get the seperator and then checks the end character.
#'     If the separator is nothing or a '/' or a '//' then it adds to
#'     the path appropriately. Without this one can unwittingly include
#'     extra separators or none at all.
#'
#' @param inpath the path to be analysed
#' @param infile the filename to be added to the path 'inpath'
#'
#' @return the completed filename or extended path
#' @export
#'
#' @examples
#' indir <- tempdir()
#' infile <- "control.csv"
#' filenametopath(indir,infile)
filenametopath <- function(inpath,infile) {
  typepath <- pathtype(inpath)
  endpath <- pathend(inpath)
  if (is.na(endpath)) {
    outfile <- paste(inpath,infile,sep=typepath)
  } else { outfile <- paste(inpath,infile,sep="")
  }
  return(outfile)
} # end of filenametopath

#' @title getConst extracts 'nb' numbers from a line of text
#'
#' @description getConst parses a line of text and extracts 'nb' pieces of
#'     text as numbers
#'
#' @param inline text line to be parsed, usually obtained using readLines
#' @param nb the number of numbers to extract
#' @param index which non-empty object to begin extracting from?
#'
#' @return a vector of length 'nb'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Not exported, prefix with AbMSE:::
#'   txtline <- "MaxDL , 32,32,32"
#'   AbMSE:::getConst(txtline,nb=3,index=2)
#' }
getConst <- function(inline,nb,index=2) { # parses lines containing numbers
  ans <- numeric(nb)
  tmp <- unlist(strsplit(inline,","))
  if (length(tmp) < (nb+1))
    warning(paste("possible problem with data",tmp[1],
                  "missing comma?",sep=" "),"\n")
  count <- 0
  for (j in index:(nb+index-1)) {
    count <- count + 1
    ans[count] <- as.numeric(tmp[j])
  }
  return(ans)
}   # end getConst


#' @title getLogical extracts nb logicals from an input line of text
#'
#' @description getLogical obtains nb logicals from an input line
#'
#' @param inline text line to be parsed, usually obtained using readLines
#' @param nb the number of logicals to extract, if nb is longer than the
#'     number of logicals within inline the vector will contain NAs
#'
#' @return a vector of length nb
#' @export
#'
#' @examples
#' \dontrun{
#'  # Not exported, prefix with AbMSE:::
#'  txtline <- "Depleted, TRUE"
#'  AbMSE:::getLogical(txtline,nb=1)
#'  txtline2 <- "calcthis, TRUE, FALSE"
#'  AbMSE:::getLogical(txtline2,nb=2)
#' }
getLogical <- function(inline,nb) {  #inline <- txtline; nb=2
  tmp <- unlist(strsplit(inline,","))
  tmp <- removeEmpty(tmp)
  outtmp <- as.logical(as.character(tmp[2:(nb+1)]))
  return(outtmp)
}

#' @title getsingleNum extracts a single number from an input line of text
#'
#' @description getsingleNum obtains a number from an input line. The
#'     variable context is the a variable declared in each function within
#'     which getsingleNum is to be called. This cvan then be used in a
#'     the stop function.
#'
#' @param varname the name of the variable to get from intxt
#' @param intxt text to be parsed, usually obtained using readLines
#' @param context the function name within which getsingleNum is
#'     being used, defaults to empty string
#'
#' @return a single number
#' @export
#'
#' @examples
#' \dontrun{
#'  # Not exported, prefix with aMSE:::
#'  context = "Function Example"
#'  txtline <- "replicates, 100"
#'  aMSE:::getsingleNum("replicates",txtline,context=context)
#'  aMSE:::getsingleNum("eeplicates",txtline,context=context)
#' }
getsingleNum <- function(varname,intxt,context="") {
  begin <- grep(varname,intxt)
  if (length(begin) > 0) {
    return(as.numeric(getConst(intxt[begin],1)))
  } else {
    stop(paste0("No data for ",varname," within ",context))
  }
}

#' @title getStr obtains a string from an input text line
#'
#' @description  getStr obtains a string from an input text line in
#'     which any parts are separated by ','. Then, after ignoring the
#'     first component, assumed to be a label, it returns the first
#'     nb parts.
#'
#' @param inline input text line with components separated by ','
#' @param nb number of parts to return
#'
#' @return a vector of character string(s)
#' @export
#'
#' @examples
#'   txt <- "runlabel, development_run, label for this particular run"
#'   getStr(txt,1)
getStr <- function(inline,nb) {
  tmp <- unlist(strsplit(inline,","))
  tmp <- removeEmpty(tmp)
  outconst <- as.character(tmp[2:(nb+1)])
  return(outconst)
} # end of getStr

#' @title pathend determines what character is at the end of a path
#'
#' @description pathend determines what character is at the end of a
#'     path uses pathtype to get the seperator and then checks the end
#'     character
#'
#' @param inpath the path to be analysed
#'
#' @return the end character of the path; either NA, '/', or "\\"
#' @export
#'
#' @examples
#'   indir <- "C:/Users/Malcolm/Dropbox/rcode2/aMSE/data-raw"
#'   pathend(indir)
pathend <- function(inpath) {
  lookfor <- pathtype(inpath)
  endpath <- NA
  if (lookfor == "/") {
    if(length(grep("/$",inpath)) > 0) endpath <- "/"
  } else {
    if(length(grep("\\\\$",inpath)) > 0) endpath <- "\\"
  }
  return(endpath)
} # end of pathend

#' @title pathtype finds the type of separator used in a path
#'
#' @description pathtype finds the type of separator used in a path,
#'     this is either a '/' or a '\\'
#'
#' @param inpath - the path to be analysed
#'
#' @return the type of path divider, either a 0 = '\\' or a
#'    1 = '/'
#' @export
#'
#' @examples
#' indir <- "C:/Users/Malcolm/Dropbox/rcode2/aMSE/data-raw"
#' pathtype(indir)
pathtype <- function(inpath) {
  typepath <- "/"
  if (length(grep("\\\\",inpath)) > 0) typepath <- "\\"
  return(typepath)
} # end of pathtype

#' @title setupdirs sets up and checks the directories used in a run
#'
#' @description setupdirs sets up and checks the directories used in
#'     a run, if they do not exist to start with they will be created.
#'
#' @param rundir the principle directory in which all all files
#'     relating to a particular run are to be held
#' @param verbose default=TRUE, prints directory status to the console
#'     if FALSE no status reports will be made to the console
#'
#' @return nothing, it allocates directly to the global environment
#' @export
#'
#' @examples
#' rund <- tempdir()
#' out <- setupdirs(rund)
#' str(out)
setupdirs <- function(rundir, verbose=TRUE) { # rundir=plotdir; runname=runname; verbose=TRUE
  dirExists(rundir,verbose=verbose)
  datadir <- filenametopath(rundir,"data")
  plotdir <- filenametopath(rundir,"plots")
  dirExists(datadir,verbose=verbose)
  dirExists(plotdir,verbose=verbose)
  return(list(datadir=datadir,plotdir=plotdir))
} # end of setupdirs

#' @title setuphtml initiates csv files lsiting results to be included
#'
#' @description setuphtml initiates the csv file used to contain the
#'     filenames, captions, and categories of the plots and tables to
#'     be included in the html results. The format of the csv file
#'     is to have column names of file, caption, category, and
#'     timestamp. Then, each plot and table is included with an entry
#'     for each column.
#'
#' @param pldir full path to the directory to contain the plots
#' @param runname the name of the particular run being summarized.
#'
#' @return full path to the plottabfile. creating the file in plotdir
#' @export
#'
#' @examples
#' indir <- tempdir()
#' plotdir <- filenametopath(indir,"plots")
#' dirExists(plotdir,verbose=FALSE)
#' plottabfile <- setuphtml(plotdir,"example_only")
#' dir(plotdir)
setuphtml <- function(pldir, runname) {  # pldir=plotdir; runname=runname
  plottabfile <- filenametopath(pldir,paste0("plotFileTable_",
                                              runname,".csv"))
  label <- c("file","category","type","timestamp","caption")
  cat(label,"\n",file = plottabfile,sep=",",append=FALSE)
  return(plottabfile)
} # end of setuphtml

