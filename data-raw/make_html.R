

add_css <- function(htmlfile) {
    cat('<html><head><title>', 'AbMSEout', '</title>\n',
      '    <!-- source for text below is http://unraveled.com/publications/css_tabs/ -->\n',
      '    <!-- CSS Tabs is licensed under Creative Commons Attribution 3.0 - http://creativecommons.org/licenses/by/3.0/ -->\n',
      '    \n',
      '    <style type="text/css">\n',
      '    \n',
      '    body {\n',
      '    font: 100% verdana, arial, sans-serif;\n',
      '    background-color: #fff;\n',
      '    margin: 50px;\n',
      '    }\n',
      '    \n',
 
      '    /* begin css tabs */\n',
      '    \n',
      '    ul#tabnav { /* general settings */\n',
      '    text-align: left; /* set to left, right or center */\n',
      '    margin: 1em 0 1em 0; /* set margins as desired */\n',
      '    font: bold 11px verdana, arial, sans-serif; /* set font as desired */\n',
      '    border-bottom: 1px solid #6c6; /* set border COLOR as desired */\n',
      '    list-style-type: none;\n',
      '    padding: 3px 10px 2px 10px; /* THIRD number must change with respect to padding-top (X) below */\n',
      '    }\n',
      '    \n',
      '    ul#tabnav li { /* do not change */\n',
      '    display: inline;\n',
      '    }\n',
      '    \n',
      '    body#tab1 li.tab1, body#tab2 li.tab2, body#tab3 li.tab3, body#tab4 li.tab4 { /* settings for selected tab */\n',
      '    border-bottom: 1px solid #fff; /* set border color to page background color */\n',
      '    background-color: #fff; /* set background color to match above border color */\n',
      '    }\n',
      '    \n',
      '    body#tab1 li.tab1 a, body#tab2 li.tab2 a, body#tab3 li.tab3 a, body#tab4 li.tab4 a { /* settings for selected tab link */\n',
      '    background-color: #fff; /* set selected tab background color as desired */\n',
      '    color: #000; /* set selected tab link color as desired */\n',
      '    position: relative;\n',
      '    top: 1px;\n',
      '    padding-top: 4px; /* must change with respect to padding (X) above and below */\n',
      '    }\n',
      '    \n',
      '    ul#tabnav li a { /* settings for all tab links */\n',
      '    padding: 2px 4px; /* set padding (tab size) as desired; FIRST number must change with respect to padding-top (X) above */\n',
      '    border: 1px solid #6c6; /* set border COLOR as desired; usually matches border color specified in #tabnav */\n',
      '    background-color: #cfc; /* set unselected tab background color as desired */\n',
      '    color: #666; /* set unselected tab link color as desired */\n',
      '    margin-right: 0px; /* set additional spacing between tabs as desired */\n',
      '    text-decoration: none;\n',
      '    border-bottom: none;\n',
      '    }\n',
      '    \n',
      '    ul#tabnav a:hover { /* settings for hover effect */\n',
      '    background: #fff; /* set desired hover color */\n',
      '    }\n',
      '    \n',
      '    /* end css tabs */\n',
      '    \n',
      '    \n',
      '    h2 {\n',
      '    font-size: 20px;\n',
      '    color: #4c994c;\n',
      '    padding-top: 1px;\n',
      '    font-weight: bold;\n',
      '    border-bottom-width: 1px;\n',
      '    border-bottom-style: solid;\n',
      '    border-bottom-color: #6c6;\n',
      '    padding-bottom: 2px;\n',
      '    padding-left: 0px;\n',
      '    }\n',
      '    </style>',
      '</head>\n',
  sep = "", file=htmlfile, append=FALSE)
} # end of add_css


#' @title make_html create a HTML files to view figures in a browser.
#' 
#' @description make_html erites a set of HTML files with tabbed navigation 
#'     between them, with each tab containing output relating to a particular
#'     set of results or diagnostics from a particular AbMSE run. This code 
#'     was borrowed and extensively modified from Ian Taylors r4ss. By default,
#'     this function will look in the directory where PNG files were created 
#'     for a plotFileTable_runname.csv file with the name 'runname' written by
#'     AbMSE. HTML files are written to link to these plots and put in the same
#'     directory.
#' 
#' @param replist Object created by an AbMSE run
#' @param plotdir Directory where .png files are saved
#' @param plotFileTable CSV file with info on PNG files. The plotdir will
#'     be sequentially searched for filenames included in this csv file,
#'     which includes a figure legend and a category = tab
#' @param width Width of plots (in pixels).
#' @param openfile Automatically open index.html in default browser?
#' @param runnotes Add additional notes to home page.
#' @param verbose Display more info while running this function?
#' @note 
#' @author Originally Ian Taylor, modified by Malcolm Haddon
#' 
#' @export
make_html <- function(replist=NULL,
                      plotdir=NULL,
                      plotFileTable="plotFileTable",
                      width=500,
                      openfile=TRUE,
                      runnotes=NULL,
                      verbose=TRUE) {
# Clarify data  
  if(is.null(plotdir)) stop("input 'plotdir' required \n")
  filenames <- dir(plotdir) 
  filetable <- filenames[grep(plotFileTable,filenames)]
  if(length(filetable)==0) stop("No plotFileTable, did the run fail? \n")
  filename <- file.path(plotdir,filetable)
  tablefile <- read.csv(filename,colClasses = "character")
  if(!is.data.frame(tablefile))
    stop("The list of files to plot needs to be a data.frame \n")  
  tablefile$png_time <- as.POSIXlt(tablefile$png_time)
  tablefile$basename <- basename(as.character(tablefile$file))
  tablefile$dirname <- plotdir
# identify the categories and name each html file
  categories <- unique(tablefile$category)
  for (icat in 0:length(categories)) { # icat=1
     if(icat==0){
       category <- "Home"
       htmlfile <- paste0(plotdir,"/AbMSEout.html")
       htmlhome <- htmlfile
       if(verbose) cat("Home HTML file with output will be:\n",htmlhome,'\n')
      }  else{
       category <- categories[icat]
       htmlfile <- paste0(plotdir,"/AbMSEout_",category,".html")
       if(verbose) cat("tab HTML file with output will be:\n",htmlfile,'\n')
     }
     
    add_css(htmlfile=htmlfile)
    
    cat('<!-- Site navigation menu -->\n',
        '  <ul id="tabnav">\n',
        file=htmlfile, append=TRUE)
    for(itab in 0:length(categories)){
      if(itab==0){
        tab <- "Home"
        cat('    <li class="tab1"><a href="AbMSEout.html">Home</a></li>\n',sep="",
            file=htmlfile, append=TRUE)
      }else{
        tab <- categories[itab]
        cat('    <li class="tab',itab+1,'"><a href="AbMSEout_',tab,'.html">',tab,'</a></li>\n',sep="",
            file=htmlfile, append=TRUE)
      }
    }
    cat('  </ul>\n', file=htmlfile, append=TRUE)
    
    if (category=="Home") {    # add text on "Home" page
      cat('\n\n<h2><a name="', category, '">', category, '</a></h2>\n', sep="",
          file=htmlfile, append=TRUE)
      MSE_info <- packageDescription("AbMSE")
      goodnames <- c("Version", "Date", "Built",grep("Remote", names(MSE_info),
                     value=TRUE))
      MSE_info_text <- '<b>AbMSE info:</b><br><br>\n'
      for(name in goodnames) {
          MSE_info_text <- c(MSE_info_text,
                              paste0(name, ": ",MSE_info[name], "<br>\n"))
      }
      cat('\n\n<p>',MSE_info_text,'</p>\n',
          '<p><b>Name of Run: </b>',replist$runname,'</p>\n',
          '<p><b>Starting time of model: </b>',
      substring(replist$starttime,12),'</p>\n\n',sep="", 
                file=htmlfile, append=TRUE)
      if (!is.null(runnotes)) {
         for(i in 1:length(runnotes)) {
            cat('<p><b>Notes:</b>\n',paste(runnotes,collapse='</b>\n'),
                '</p>\n\n',sep="", file=htmlfile, append=TRUE)
          }
      } # end of runnotes
    } else {   #  split on category if statement
      plotinfo <- tablefile[tablefile$category==category,]
      cat('\n\n<h2><a name="', category, '">', category, '</h2>\n', sep="",
          file=htmlfile, append=TRUE)
      for(i in 1:nrow(plotinfo)){  # i=1
        cat("<p align=left><a href='", plotinfo$basename[i],
            "'><img src='", plotinfo$basename[i],
            "' border=0 width=", width, "></a><br>",
            plotinfo$caption[i],
            "<br><i><small>file: <a href='", plotinfo$basename[i],
            "'>", plotinfo$basename[i], "</a></small></i>\n",
            sep="",  file=htmlfile,  append=TRUE)
      }
    } # end of category if else statement

    cat("\n\n</body>\n</html>", file=htmlfile, append=TRUE)
  }
  # open HTML file automatically:
  if(openfile) browseURL(htmlhome)
  
} # end of make_html tabs
  
  
  
  
  
  abmselist <- list(runname="Test_HTML",starttime=Sys.time(),
                    Nwarnings=NA,runnotes=NULL)
  plotdir="C:/Users/Malcolm/Dropbox/AbMSE/AbMSEuse/plots"

  
  
    # Need some indication of a successful AbMSE run.
  
  make_html(replist=abmselist,plotdir=plotdir,width=750)
  
  replist=abmselist
  plotdir=plotdir
  title=title
  plotFileTable = "plotFileTable"
  width = 500
  openfile=TRUE
  runnotes = NULL
  verbose = TRUE
  
  
  
  