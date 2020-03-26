




#' Title
#'
#' @param pop
#' @param res
#' @param location
#'
#' @return
#' @export
#'
#' @examples
findF1 <- function(pop,res,location=TRUE) {
  harv <- results[,"AnnH",pop]
  catch <- results[,"Catch",pop]
  grad <- numeric(nH-1)
  for (i in 1:(nH-1)) {
    divisor <- harv[i+1] - harv[i]
    numerator <- catch[i+1] - catch[i]
    grad[i] <- numerator/divisor
  }
  pickF1 <-  which.closest(0.1,grad/grad[1])
  if (location) {
    return(pickF1)
  } else {
    return(grad)
  }
}








