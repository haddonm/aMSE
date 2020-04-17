


#' @title write_css generates a CSS file used by all html files
#'
#' @description write_css generates a cascading style sheet that will
#'     be used by each separate html files made for each category of
#'     results. The origin came from Ian Taylor although I have
#'     modified the styles and ensure correct Version 3 CSS
#'
#' @param plotdir the directory within the run directory that contains
#'     all the plot results
#'
#' @return nothing but it does generate a .css file in the plotdir
write_css <- function(plotdir) {
  filename <- filenametopath(plotdir,"aMSEout.css")
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




#' @title make_html create HTML files to view results in a browser.
#'
#' @description make_html writes a set of HTML files with tabbed
#'     navigation, with each tab containing the results relating to a
#'     given set of results or diagnostics from a particular aMSE run.
#'     This code was borrowed from Ian Taylor's r4ss, but has been
#'     extensively modified to improve both the css (see write_css)
#'     and the HTML. By default, this function will look in the
#'     plots directory where PNG files were created for a
#'     plotFileTable_runname.csv file with the name 'runname' written
#'     added by aMSE. HTML files are written to link to these plots
#'     and put in the same directory.
#'
#' @param replist Object created by an aMSE run
#' @param rundir Directory where a particular run's files,
#'     including the region files, the results, as tables and plots,
#'     and any other files, are all held. It will always contain at
#'     least the 'data' and 'plots' sub-directories.
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
  # Clarify data
  plotdir <- filenametopath(rundir,"plots")
  if(is.null(plotdir)) stop("input 'plotdir' required \n")
  write_css(plotdir)
  filenames <- dir(plotdir)
  filetable <- filenames[grep("plotFileTable",filenames)]
  if(length(filetable)==0) stop("No plotFileTable, did the run fail? \n")
  filename <- filenametopath(plotdir,filetable)
  tablefile <- read.csv(filename,colClasses = "character")
  if(!is.data.frame(tablefile))
    stop("The list of files to plot needs to be a data.frame \n")
  #tablefile$TimeMade <- as.POSIXlt(tablefile$png_time)
  tablefile$basename <- basename(as.character(tablefile$file))
  tablefile$dirname <- plotdir
  # identify the categories and name each html file
  categories <- unique(tablefile$category)
  for (icat in 0:length(categories)) { # icat=1
    if(icat==0){
      category <- "Home"
      htmlfile <- paste0(plotdir,"/aMSEout.html")
      htmlhome <- htmlfile
      if(verbose) cat("Home HTML file with output will be:\n",htmlhome,'\n')
    }  else{
      category <- categories[icat]
      htmlfile <- paste0(plotdir,"/aMSEout_",category,".html")
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
      cat('\n\n<h2><a name="', category, '">', category, '</a></h2>\n', sep="",
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
        cat("<p align=left><a href='", plotinfo$basename[i],
            "'><img src='", plotinfo$basename[i],
            "' border=0 width=", width, "></a><br>",
            plotinfo$caption[i],
            "<br><i>file: <a href='", plotinfo$basename[i],
            "'>", plotinfo$basename[i], "</a></i></p>\n\n",
            sep="",  file=htmlfile,  append=TRUE)
      }
    } # end of category if else statement

    cat("\n\n</body>\n</html>", file=htmlfile, append=TRUE)
  }
  # open HTML file automatically:
  if(openfile) browseURL(htmlhome)
} # end of make_html2 tabs


