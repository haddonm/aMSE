

#' @title checkctrldat checks datadir contains the required csv files
#'
#' @description checkctrldat checks datadir contains the required csv
#'     files including the named control file, which then contains
#'     the names of the region data file, and the population data
#'     file. The run stops if any are not present or are misnamed.
#'
#' @param datadir the directory in which all all files relating to a
#'     particular run are to be held.
#' @param ctrlfile default="control.csv", the filename of the control
#'     file present in datadir containing information regarding the
#'     run.
#'
#' @return the control list for the run
#' @export
#'
#' @examples
#' resdir <- tempdir()
#' ctrlfiletemplate(resdir)
#' zonefiletemplate(resdir)
#' datafiletemplate(6,resdir,filename="zone1sau2pop6.csv")
#' ctrl <- checkctrldat(resdir)
#' ctrl
checkctrldat <- function(datadir,ctrlfile="control.csv") { # resdir=resdir; ctrlfile="control.csv"
  filenames <- dir(datadir)
  if (length(grep(ctrlfile,filenames)) != 1)
    stop(cat(ctrlfile," not found in resdir \n"))
  ctrol <- readctrlfile(datadir,infile=ctrlfile)
  if (length(grep(ctrol$zonefile,filenames)) != 1)
      stop("zone data file not found \n")
  if (length(grep(ctrol$datafile,filenames)) != 1)
    stop("population data file not found \n")
  cat("All required files appear to be present \n")
  return(ctrol)
} # end of checctrldat


