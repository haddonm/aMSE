

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
#' data(product)
#' zoneprod <- summarizeprod(product)
#' round(zoneprod,3)
summarizeprod <- function(product,saunames) { # product=out$product; saunames=zone1$SAUnames
  numrow <- dim(product)[1]
  prodrows <- rownames(product[,,1])
  prodcols <- colnames(product[,,1])
  numpop <- dim(product)[3]  # here numpop = nSAU
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
