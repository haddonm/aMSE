

# rutilsMH::listFunctions("C:/Users/User/Dropbox/A_code/aMSE/aMSE_utils.R")

#' @title summarizeprod generates a summary of the productivity properties
#'
#' @description summarizeprod generates a summary of the productivity properties
#'     by examining the productivity matrix for each SAU/population and
#'     extracting the Bmsy, annualH, MSY, Depletion, and RelCE at the maximum
#'     catch level (which approximates the MSY). It summarizes the total zone
#'     by summing all the productivity matrices and search for the largest
#'     catch again. It generates estimates of the annualH, depletion and RelCE
#'     by using a weighted average of those values from teh separate SAU or
#'     populations, where the weighting is the proportion of the sum of the
#'     MSYs taken in each sau or population. This latter is only an
#'     approximation but provides at least an indication.
#'
#' @param product The productivity array from doproduction containing the
#'     range of imposed harvest rates, and the resulting outputs for each
#'     population
#' @param saunames the names of the different SAU
#'
#' @return a matrix containing the approximate productivity matrix for the zone
#' @export
#'
#' @examples
#' \dontrun{
#' data(zone)
#' product <- zone$product
#' zoneprod <- summarizeprod(product,saunames=zone$zone1$SAUnames)
#' round(zoneprod,3)
#' }
summarizeprod <- function(product,saunames) { # product=zone$product; saunames=zone$zone1$SAUnames
  numrow <- dim(product)[1]
  numpop <- dim(product)[3]
  prodrows <- rownames(product[,,1])
  prodcols <- colnames(product[,,1])
 # nSAU <- length()  # here numpop = nSAU
  columns <- c("Bmsy","AnnH","MSY","Deplet","RelCE")
  rows <- c(saunames,"zone")
  nrows <- length(rows)
  ans <- matrix(0,nrow=nrows,ncol=length(columns),dimnames=list(rows,columns))
  for (sau in 1:numpop) {
    sauprod <- product[,,sau]
    pick <- which(sauprod[,"Catch"] == max(sauprod[,"Catch"],na.rm=TRUE))
    ans[sau,] <- sauprod[pick,2:6]
  }
  zoneprod <- matrix(0,nrow=numrow,ncol=length(prodcols),
                     dimnames=list(prodrows,prodcols))
  for (i in 1:numpop) zoneprod <- zoneprod + product[,,i]
  pick <-  which(zoneprod[,"Catch"] == max(zoneprod[,"Catch"],na.rm=TRUE))
  ans[nrows,] <- zoneprod[pick,2:6]
  pmsy <- ans[1:numpop,"MSY"]/sum(ans[1:numpop,"MSY"],na.rm=TRUE)
  ans[nrows,"AnnH"] <- sum(ans[1:numpop,"AnnH"] * pmsy,na.rm=TRUE)
  ans[nrows,"Deplet"] <- sum(ans[1:numpop,"Deplet"] * pmsy,na.rm=TRUE)
  ans[nrows,"RelCE"] <- sum(ans[1:numpop,"RelCE"] * pmsy,na.rm=TRUE)
  return(ans)
} # end of summarizeprod


#' @title wtedmean calculates the weighted mean of a set of values and weights
#'
#' @description wtedmean solves the problem of calculating a weighted mean
#'     value from a set of values with different weights. Within the aMSE this
#'     is common when trying to summarize across populations within an SAU or
#'     summarize SAU within a zone by finding a mean value weighted by the
#'     respective catch from each related population or SAU.
#'
#' @param x the values whose weighted mean is wanted
#' @param wts the weights to use, often a set of catches
#'
#' @return a single real number
#' @export
#'
#' @examples
#' saucpue <- c(91.0,85.5,88.4,95.2)
#' saucatch <- c(42.0,102.3,75.0,112.0)
#' wtedmean(saucpue,saucatch)
#' saucatch/sum(saucatch)  # the relative weights
wtedmean <- function(x,wts) {
  pwts <- wts/sum(wts,na.rm=TRUE)
  ans <- sum((x * pwts),na.rm=TRUE)
  return(ans)
} # end of wtedmean

