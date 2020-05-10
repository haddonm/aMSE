

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
setupdirs <- function(rundir, verbose=TRUE) { # rundir=resdir; runname=runname; verbose=TRUE
  dirExists(rundir,verbose=verbose)
  datadir <- filenametopath(rundir,"data")
  resdir <- filenametopath(rundir,"results")
  dirExists(datadir,verbose=verbose)
  dirExists(resdir,verbose=verbose)
  return(list(datadir=datadir,resdir=resdir))
} # end of setupdirs
