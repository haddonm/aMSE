

# rutilsMH::listFunctions("C:/Users/User/Dropbox/A_code/aMSE/aMSE_utils.R")

#' @title alltosau sums population properties to form SAU totals
#'
#' @description alltosau sums the projected values for matureB, exploitB,
#'     acatch, catch, and recruit, by population, into totals by SAU, using
#'     sauindex from glb
#'
#' @param zoneDR the dynamics projeciton object
#' @param glb the global constants object
#'
#' @return a list of some SAU properties
#' @export
#'
#' @examples
#' print("wait on appropriate internal data sets")
alltosau <- function(zoneDR,glb) {  # zoneDR=zoneDP; glb=glb
  saunames <- glb$saunames
  matB <- poptosau(zoneDR$matureB,glb=glb)
  exB <- poptosau(zoneDR$exploitB,glb=glb)
  recruit <- poptosau(zoneDR$recruit,glb=glb)
  catsau=zoneDR$catsau
  dimnames(catsau)[[2]] <- saunames
  acatch=zoneDR$acatch
  dimnames(acatch)[[2]] <- saunames
  cesau=zoneDR$cesau
  dimnames(cesau)[[2]] <- saunames
  sumsau <- list(matureB=matB,exploitB=exB,catsau=catsau,acatch=acatch,
                 cesau=cesau,recruit=recruit)
  return(sumsau)
} # end of alltosau


#' @title poptosau converts projected population dynamics to SAU scale results
#'
#' @description poptosau the MSE dynamics are run at the population level but
#'     all management decisions are made at the SAU level, hence the results of
#'     the projections need to be translated into results at the SAU level.
#'     poptosau uses the sauindex to sum the variables that can be summed, which
#'     include matureB, exploitB, catch, and recruit, and combines their
#'     population results to form their respective SAU results. Thus, for an
#'     input array of dimensions [projyrs,numpop,reps] one receives an array of
#'     [projyrs, nSAU,reps]
#'
#' @param invar either zoneDP$ matureB, exploitB, catch, or recruit
#' @param glb the global constants object
#'
#' @return a results array of dimnesion [projyrs, nSAU,reps]
#' @export
#'
#' @examples
#' print("wait on appropriate internal datasets")
poptosau <- function(invar,glb) {  # invar=zoneDP$matureB; glb=glb
  numpop <- glb$numpop
  nsau <- glb$nSAU
  nyrs <- dim(invar)[1]
  reps <- dim(invar)[3]
  sauindex <- glb$sauindex
  saunames <- glb$saunames
  result <- array(0,dim=c(nyrs,nsau,reps),
                  dimnames=list(1:nyrs,saunames,1:reps)) #aspirational catches
  for (iter in 1:reps)
    for (yr in 1:nyrs)
      result[yr,,iter] <- tapply(invar[yr,,iter],sauindex,sum,na.rm=TRUE)
  return(result)
} # end of poptosa


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

#' @title sumpop2sau gathers population data into sau data using sauindex
#'
#' @param invect a vector of population values for a given variable
#' @param sauindex the indices of each sau for each population
#'
#' @return a vector of length nsau containing the sum of population values
#'     for each sau
#' @export
#'
#' @examples
#' vect <- c(5.8,6.2,13.2,23.8,3.3,3.7,29.7,26.3,38.9,9.1)
#' sauind <- c(1,1,2,2,3,3,4,4,5,5)
#' sumpop2sau(vect,sauind)   # should be 12 37 7 56 48
sumpop2sau <- function(invect,sauindex) {
  return(tapply(invect,sauindex,sum,na.rm=TRUE))
} # end of sumpop2sau

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

