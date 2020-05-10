
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

#' @title htmltable generates the html to print out a table
#'
#' @description htmltable generates the required html code to add a
#'     table to the website of the results. This requires the
#'     logfilename function to include both a category biology,
#'     productivity, etc, and a type, which is currently limited to
#'     plot or table.
#'
#' @param inmat the 2D matrix or data.frame to be printed
#' @param filename the filename to which to add the html, defined by
#'     makehtml
#' @param caption the caption text placed at the top of the table
#'
#' @return nothing but it does add some html to the input filename
#' @export
htmltable <- function(inmat,filename,caption) {
  rows <- rownames(inmat)
  numrow <- length(rows)
  columns <- colnames(inmat)
  numcol <- length(columns)
  cat('<br><br> \n',file=filename,append=TRUE)
  cat('<table> \n',file=filename,append=TRUE)
  cat('<caption> ',caption,'</caption> \n',file=filename,append=TRUE)
  cat('<br> \n',file=filename,append=TRUE)
  cat('<tr> \n',file=filename,append=TRUE)
  cat(' <th>Var</th> \n',file=filename,append=TRUE)
  for (cl in 1:numcol)
    cat(' <th>',columns[cl],'</th> \n',file=filename,append=TRUE)
  for (rw in 1:numrow) {
    if ((rw %% 2) == 0) {
      cat('<tr> \n',file=filename,append=TRUE)
    } else {
      cat('<tr class="odd"> \n',file=filename,append=TRUE)
    }
    cat(' <th>',rows[rw],'</th> \n',file=filename,append=TRUE)
    for (cl in 1:numcol)
      cat(' <th>',inmat[rw,cl],'</th> \n',file=filename,append=TRUE)
    cat('</tr> \n',file=filename,append=TRUE)
  }
  cat('</table>',file=filename,append=TRUE)
} # end of htmltable

#' @title logfilename adds a filename to the autoresult csv file
#'
#' @description logfilename is used to facilitate the produciton of
#'     the HTML results summary for a particular run. This depends
#'     upon a csv file containing the names of each file to be
#'     plotted or tabulated. This function adds a filename and the
#'     supporting caption and category, without one needing to
#'     remember the syntax. If no category is added explicitly then
#'     the local webpage will have an 'any' tab containing these
#'     unloved results.
#'
#' @param filename the full path and filename for the file being added
#' @param resfile the file to be added to, which is defined by
#'     setuphtml found in aMSE_utils
#' @param category what HTML tab should it be added to? default="any"
#' @param caption the caption for the figure or table, default = ""
#'
#' @return nothing but it does add a line to resfile
#' @export
#'
#' @examples
#' indir <- tempdir()
#' resdir <- filenametopath(indir,"result")
#' dirExists(resdir,verbose=FALSE)
#' resfile <- setuphtml(resdir,"example_only")
#' filename <- filenametopath(resdir,"example.png")
#' png(filename=filename,width=7,height=4,units="in",res=300)
#' plot(runif(100),runif(100),type="p")
#' graphics.off()  # could use dev.off()
#' logfilename(filename=filename,resfile=resfile,"A_category",
#'             caption="Example Figure")
#' dir(resdir)
logfilename <- function(filename,resfile,category="any",caption="") {
  type <- getextension(filename)
  cat(c(filename,category,type,as.character(Sys.time()),caption," \n"),
      file=resfile,sep=",",append=TRUE)
}

#' @title make_html create HTML files to view results in a browser.
#'
#' @description make_html writes a set of HTML files with tabbed
#'     navigation, with each tab containing the results relating to a
#'     given set of results or diagnostics from a particular aMSE run.
#'     This code was borrowed from Ian Taylor's r4ss, but has been
#'     extensively modified to improve both the css (see write_css)
#'     and the HTML. By default, this function will look in the
#'     results directory where PNG files were created for a
#'     resultTable_runname.csv file with the name 'runname' written
#'     added by aMSE. HTML files are written to link to these plots
#'     and put in the same directory.
#'
#' @param replist Object created by an aMSE run
#' @param rundir Directory where a particular run's files,
#'     including the region files, the results, as tables and plots,
#'     and any other files, are all held. It will always contain at
#'     least the 'data' and 'results' sub-directories.
#' @param width Width of plots (in pixels). Default = 500
#' @param openfile Automatically open index.html in default browser?
#' @param runnotes Add additional notes to home page.
#' @param verbose Display more info while running this function?
#' @author Originally Ian Taylor, modified by Malcolm Haddon
#'
#' @export
make_html <- function(replist=NULL,
                      rundir=NULL,
                      width=500,
                      openfile=TRUE,
                      runnotes=NULL,
                      verbose=TRUE) {
  # replist=reportlist;rundir=rundir;width=500;openfile=TRUE;runnotes=runnotes;verbose=FALSE
  # Clarify data
  resdir <- filenametopath(rundir,"results")
  if(is.null(resdir)) stop("input 'resdir' required \n")
  write_css(resdir)
  filenames <- dir(resdir)
  filetable <- filenames[grep("resultTable",filenames)]
  if(length(filetable)==0) stop("No resultTable, did the run fail? \n")
  filename <- filenametopath(resdir,filetable)
  tablefile <- read.csv(filename,colClasses = "character")
  if(!is.data.frame(tablefile))
    stop("The list of files to plot needs to be a data.frame \n")
  #tablefile$TimeMade <- as.POSIXlt(tablefile$png_time)
  tablefile$basename <- basename(as.character(tablefile$file))
  tablefile$dirname <- resdir
  # identify the categories and name each html file
  categories <- unique(tablefile$category)
  types <- tablefile$type
  for (icat in 0:length(categories)) { # icat=1
    if(icat==0){
      category <- "Home"
      htmlfile <- paste0(resdir,"/aMSEout.html")
      htmlhome <- htmlfile
      if(verbose) cat("Home HTML file with output will be:\n",htmlhome,'\n')
    }  else{
      category <- categories[icat]
      htmlfile <- paste0(resdir,"/aMSEout_",category,".html")
      if(verbose) cat("tab HTML file with output will be:\n",htmlfile,'\n')
    }
    write_head(htmlfile)
    cat('<body> \n',file=htmlfile, append=TRUE)
    cat('<!-- Site navigation menu -->\n',
        '  <ul id="tabnav">\n',file=htmlfile, append=TRUE)
    for(itab in 0:length(categories)){
      if(itab==0){
        tab <- "Home"
        cat('    <li class="tab1"><a href="aMSEout.html">Home</a></li>\n',sep="",
            file=htmlfile, append=TRUE)
      }else{
        tab <- categories[itab]
        cat('    <li class="tab',itab+1,'"><a href="aMSEout_',tab,'.html">',tab,'</a></li>\n',sep="",
            file=htmlfile, append=TRUE)
      }
    }
    cat('  </ul>\n', file=htmlfile, append=TRUE)

    if (category=="Home") {    # add text on "Home" page
      newcat <- "Run Details"
      cat('\n\n<h2><a name="', category, '">', newcat, '</a></h2>\n', sep="",
          file=htmlfile, append=TRUE)
      MSE_info <- packageDescription("aMSE")
      goodnames <- c("Version", "Date", "Built",grep("Remote", names(MSE_info),
                                                     value=TRUE))
      MSE_info_text <- '<b>aMSE info:</b><br><br>\n'
      for(name in goodnames) {
        MSE_info_text <- c(MSE_info_text,
                           paste0(name, ": ",MSE_info[name], "<br>\n"))
      }
      cat('\n\n<p>',MSE_info_text,'</p>\n',
          '<p><b>Name of Run: </b>',replist$runname,'</p>\n',
          '<p><b>Directory: </b>',rundir,'</p>\n',
          '<p><b>Starting time of model: </b>',
          replist$starttime,'</p>\n',
          '<p><b>End time of model: </b>',
          replist$endtime,'</p>\n\n',
          sep="",file=htmlfile, append=TRUE)

      if (!is.null(runnotes)) {
        for(i in 1:length(runnotes)) {
          cat('<p><b>Notes:</b>\n',paste(runnotes,collapse='</b>\n'),
              '</p>\n\n',sep="", file=htmlfile, append=TRUE)
        }
      } # end of runnotes
    } else {   #  split on category if statement
      plotinfo <- tablefile[tablefile$category==category,]
      cat('\n\n<h2><a name="', category, '">', category, '</a></h2>\n', sep="",
          file=htmlfile, append=TRUE)
      for(i in 1:nrow(plotinfo)){  # i=1
        if (plotinfo$type[i] == "plot") {
          cat("<p align=left><a href='", plotinfo$basename[i],
              "'><img src='", plotinfo$basename[i],
              "' border=0 width=", width, "></a><br>",
              plotinfo$caption[i],
              "<br><i>file: <a href='", plotinfo$basename[i],
              "'>", plotinfo$basename[i], "</a></i></p>\n\n",
              sep="",  file=htmlfile,  append=TRUE)
        }
        if (plotinfo$type[i] == "table") {
          datafile <- filenametopath(resdir,plotinfo$basename[i])
          dat <- read.csv(file=datafile,header=TRUE,row.names=1)
          htmltable(inmat=dat,filename=htmlfile,caption=plotinfo$caption[i])
        }
      }
    } # end of category if else statement

    cat("\n\n</body>\n</html>", file=htmlfile, append=TRUE)
  }
  # open HTML file automatically:
  if(openfile) browseURL(htmlhome)
} # end of make_html2 tabs

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

#' @title setuphtml initiates csv files lsiting results to be included
#'
#' @description setuphtml initiates the csv file used to contain the
#'     filenames, captions, and categories of the plots and tables to
#'     be included in the html results. The format of the csv file
#'     is to have column names of file, caption, category, and
#'     timestamp. Then, each plot and table is included with an entry
#'     for each column.
#'
#' @param resdir full path to the directory to contain the plots
#' @param runname the name of the particular run being summarized.
#'
#' @return full path to the resfile. creating the file in resdir
#' @export
#'
#' @examples
#' indir <- tempdir()
#' resdir <- filenametopath(indir,"results")
#' dirExists(resdir,verbose=FALSE)
#' resfile <- setuphtml(resdir,"example_only")
#' dir(resdir)
setuphtml <- function(resdir, runname) {  # resdir=resdir; runname=runname
  resfile <- filenametopath(resdir,paste0("resultTable_",
                                          runname,".csv"))
  label <- c("file","category","type","timestamp","caption")
  cat(label,"\n",file = resfile,sep=",",append=FALSE)
  return(resfile)
} # end of setuphtml

#' @title write_css generates a CSS file used by all html files
#'
#' @description write_css generates a cascading style sheet that will
#'     be used by each separate html files made for each category of
#'     results. The origin came from Ian Taylor although I have
#'     modified the styles and ensure correct Version 3 CSS
#'
#' @param resdir the directory within the run directory that contains
#'     all the plot results
#'
#' @return nothing but it does generate a .css file in the resdir
write_css <- function(resdir) {
  filename <- filenametopath(resdir,"aMSEout.css")
  cat('    \n',
      '    body {\n',
      '      font-size:  18px; \n',
      '      font-family: Cambria, "Hoefler Text", "Liberation Serif", Times, "Times New Roman", serif;\n',
      '      background-color: #fff;\n',
      '      margin: 50px;\n',
      '    }\n',
      '    \n',
      '    /* begin css tabs */\n',
      '    \n',
      '    #tabnav { /* general settings */\n',
      '      text-align: left; /* set to left, right or center */\n',
      '      margin: 1em 0 1em 0; /* set margins as desired */\n',
      '      font: bold 11px verdana, arial, sans-serif; /* set font as desired */\n',
      '      border-bottom: 1px solid #6c6; /* set border COLOR as desired */\n',
      '      list-style-type: none;\n',
      '      padding: 3px 10px 2px 10px; /* THIRD number must change with respect to padding-top (X) below */\n',
      '    }\n',
      '    \n',
      '    #tabnav li { /* do not change */\n',
      '      display: inline;\n',
      '    }\n',
      '    \n',
      '    #tab1 li.tab1, #tab2 li.tab2, #tab3 li.tab3, #tab4 li.tab4 { /* settings for selected tab */\n',
      '      border-bottom: 1px solid #fff; /* set border color to page background color */\n',
      '      background-color: #fff; /* set background color to match above border color */\n',
      '    }\n',
      '    \n',
      '    #tab1 li.tab1 a, #tab2 li.tab2 a, #tab3 li.tab3 a, #tab4 li.tab4 a { /* settings for selected tab link */\n',
      '      background-color: #fff; /* set selected tab background color as desired */\n',
      '      color: #000; /* set selected tab link color as desired */\n',
      '      position: relative;\n',
      '      top: 1px;\n',
      '      padding-top: 4px; /* must change with respect to padding (X) above and below */\n',
      '    }\n',
      '    \n',
      '    #tabnav li a { /* settings for all tab links */\n',
      '      padding: 2px 4px; /* set padding (tab size) as desired; FIRST number must change with respect to padding-top (X) above */\n',
      '      border: 1px solid #6c6; /* set border COLOR as desired; usually matches border color specified in #tabnav */\n',
      '      background-color: #cfc; /* set unselected tab background color as desired */\n',
      '      color: #666; /* set unselected tab link color as desired */\n',
      '      margin-right: 0px; /* set additional spacing between tabs as desired */\n',
      '      text-decoration: none;\n',
      '      border-bottom: none;\n',
      '    }\n',
      '    \n',
      '    #tabnav a:hover { /* settings for hover effect */\n',
      '      background: #fff; /* set desired hover color */\n',
      '    }\n',
      '    \n',
      '    /* end css tabs */\n',
      '    \n',
      '    \n',
      '    h2 {\n',
      '      font-size: 20px;\n',
      '      color: #1611A7;\n',
      '      padding-top: 1px;\n',
      '      font-weight: bold;\n',
      '      border-bottom-width: 1px;\n',
      '      border-bottom-style: solid;\n',
      '      border-bottom-color: #6c6;\n',
      '      padding-bottom: 2px;\n',
      '      padding-left: 0px;\n',
      '    }\n',
      '    .odd {\n',
      '      background-color: #cfc;\n ',
      '    } \n',
      '    table { \n',
      '      table-layout: fixed; \n',
      '      font-size: 15px; \n',
      '      border: 1px solid black;\n',
      '      border-collapse: collapse; \n',
      '    } \n',
      '    th, td {\n',
      '      padding: 5px; \n',
      '      border: 1px solid black;\n',
      '      border-collapse: collapse; \n',
      '    }\n',
      '    th { \n',
      '       text-align: right; \n',
      '    }\n',
           sep = "", file=filename, append=FALSE)
} # end of write_css

#' @title write_head adds the <head> tag to  each html file in results
#'
#' @description write_head adds the head tag to each html file.
#'     It links to the css file rather than writing the complete
#'     css code to every file. There were missing elements which are
#'     now in place, so the code is now valid HTML5.
#'
#' @param htmlfile the particular html file being worked on. This is
#'     defined within the make_html function.
#'
#' @return nothing but it does add the <head> tag to each html file
write_head <- function(htmlfile) {
  cat('<!DOCTYPE html> \n',
      '<html> \n',
      '  <head>',
      '    <meta charset="utf-8"> \n',
      '    <meta name="format-detection" content="telephone=no"/> \n',
      '    <title>', 'aMSEout', '</title>\n',
      '    <!-- source for text below is http://unraveled.com/publications/css_tabs/ -->\n',
      '    <!-- CSS Tabs is licensed under Creative Commons Attribution 3.0 - http://creativecommons.org/licenses/by/3.0/ -->\n',
      '    <!-- When visiting unraveled.com/publications/css_tabs it appeared to be a toxic website - BE CAREFUL -->\n',
      '    \n',
      '    <link href="aMSEout.css" rel="stylesheet" type="text/css"> \n',
      '    \n',
      '  </head>\n',
      sep = "", file=htmlfile, append=FALSE)
} # end of write_head




