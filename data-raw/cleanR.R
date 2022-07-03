




#' @title alldirExists Checks the existence of both a run and data directory
#'
#' @description alldirExists answers the questions 'do both a rundir and a
#'     datadir directory exist?' It uses dir.exists and reports existence if
#'     already present and provides a stop warning if either does not exist. Of
#'     course it can also be used to determine whether a single directory exists.
#'     By default the second directory, indir2 = indir1. This allows for the
#'     datadir = rundir within aMSE, but also allows for a separate datadir.
#'     was in aMSE_utils.R
#'
#' @param indir1 a character string containing the name of the first directory
#'     whose existence is to be checked before it is created if it
#'     does not already exist.
#' @param indir2 a character string containing the name of a second directory
#'     whose existence is to be checked, by default this has the same value as
#'     indir1
#' @param make if the directory does NOT exist should it be created.
#'     default = FALSE; if make=FALSE and a directory does not exist
#'     a warning will be given to the console.
#' @param verbose default=TRUE, prints directory status to the console,
#'     If make is set to FALSE and a directory does not exist a
#'     warning will always be given.
#'
#' @return a message to the screen if the directory exists or is
#'     created; if make is TRUE then it also creates the directory as
#'     listed in 'indir1'.
#' @export
#'
#' @examples
#' indirect <- getwd()
#' alldirExists(indirect)
alldirExists <- function(indir1,indir2=indir1,make=FALSE,verbose=TRUE) {
  if (dir.exists(indir1)) {
    if (verbose) cat("rundir, ",indir1,":  exists  \n")
  } else {
    if (make) {
      dir.create(indir1, recursive = TRUE)
      if (verbose) cat(indir1,":  created  \n")
    } else {
      warning(cat(indir1,":  does not exist \n"))
    }
  }
  if (indir2 != indir1) {
    if (dir.exists(indir2)) {
      if (verbose) cat("datadir, ",indir2,":  exists  \n")
    } else {
      if (make) {
        dir.create(indir2, recursive = TRUE)
        if (verbose) cat("datadir, ",indir2,":  created  \n")
      } else {
        warning(cat("datadir, ",indir2,":  does not exist \n"))
      }
    }
  } # end of datadir != rundir
}  # end of alldirExists

































