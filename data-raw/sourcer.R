



#' Title
#'
#' @param plotdir
#' @param tabledir
#' @param runname
#'
#' @return
#' @export
#'
#' @examples
setuphtml <- function(pldir, tabdir,rname) {
  plottabfile <<- filenametopath(pldir,paste0("plotFileTable_",rname,".csv"))
  tabletabfile <<- filenametopath(tabdir,paste0("TableFileTable_",rname,".csv"))
  label <- c("file","caption","category","TimeMade")
  cat(label,"\n",file = plottabfile,sep=",")
  cat(label,"\n",file = tabletabfile,sep=",")
} # end of setuphtml





