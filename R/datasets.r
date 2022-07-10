
#' @title lfs is a 3-D array of size-composition data for use in examples
#'
#' @description lfs is a 38 x 31 x 8 array of values for size-composition data
#'     from the west coast of Tasmania made up of 38 size classes, 31 years,
#'     from 1990 - 2020, and 8 statistical blocks. This is purely for use with
#'     examples within the aMSE R package and requires quality control checks
#'     before use anywhere else.
#'
#' @name lfs
#'
#' @docType data
#'
#' @examples
#'  data(lfs)
#'  lfs[,10:20,6]
NULL

#' @title midg is an abalone tagging data-set from the Actaeons
#'
#' @description midg is a tagging data-set for blacklip abalone
#'     (\emph{Haliotis rubra}) from the middle ground in the Actaeons
#'     in Tasmania's Block 13. All individuals were recaptured in 2003,
#'     the site number was 478, at Latitude -43.54 longitude 146.99,
#'     and there are 347 observations. All Dt = 1 year.
#'
#' @name midg
#'
#' @docType data
#'
#' @format A data.frame of abalone tagging data
#' \describe{
#'   \item{RecapL}{the length at recapture}
#'   \item{Lt}{the length at tagging}
#'   \item{Dt}{The time interval between tagging and recapture, in this
#'       instance they are all listed as 1 year}
#'   \item{DL}{the growth increment in mm}
#' }
#'
#' @section Subjects:
#'  \itemize{
#'    \item growth curves
#'    \item inverse logistic, von Bertalanffy, Gompertz
#'    \item Static model fitting
#'  }
#'
#' @source Thanks to the Institute of Marine and Antarctic Science,
#'     which is part of the University of Tasmania, and especially to
#'     Dr Craig Mundy, leader of the Abalone Group, for permission to use
#'     this data collected in 2003.
#'
#' @examples
#'  data(midg)
#'  head(midg,20)
#'  oldpar <- par(no.readonly=TRUE)
#'  plot(midg$Lt,midg$DL,type="p",pch=16,cex=1.0,xlim=c(5,180))
#'  abline(h=0,col=1)
#'  par(oldpar)
NULL

#' @title saudat is a matrix of constants read by readsaudatafile
#'
#' @description saudat is a 32 x 2 matrix of values for biological variables
#'     relating to growth, weight-at-size, maturity-at-size, average unfished
#'     recruitment, selectivity, and other details, including the variation
#'     expected for each main variable. This information is used to define the
#'     populations within each SAU.
#'
#' @name saudat
#'
#' @docType data
#'
#' @format A data.frame of biological and fishery parameters
#' \describe{
#'   \item{DLMax}{maximum growth increment}
#'   \item{L50}{size at half the maximum growth increment}
#'   \item{L50inc}{L95 - L50, size where growth increment only 5percent of max}
#'   \item{Wtb}{exponent of weight-at-size relationship}
#'   \item{Wta}{intercept of weight-at-size relationship}
#' }
#'
#'
#' @examples
#'  data(saudat)
#'  round(saudat,5)
NULL

#' @title tasab is a matrix of abalone maturity-at-length data
#'
#' @description tasab is a 715 x 4 matrix of maturity-at-length data
#'     for blacklip abalone (\emph{Haliotis rubra}) from two sites
#'     along the Tasmanian west coast. All data was collected in
#'     February 1995, but details, such as site name, accurate
#'     location, statistical block, year, month, and other
#'     details have been omitted for brevity.
#'
#' @name tasab
#'
#' @docType data
#'
#' @format A data.frame of maturity-at-length data
#' \describe{
#'   \item{site}{an identifier for the two different sites sampled}
#'   \item{sex}{I = immature, M = male, F = female}
#'   \item{length}{the shell length in mm}
#'   \item{mature}{was the animal mature = 1 or not = 0}
#' }
#'
#' @section Subjects:
#'  \itemize{
#'    \item maturity ogives or logistic curves
#'    \item Binomial likelihoods
#'  }
#'
#' @source Thanks to the Institute of Marine and Antarctic Science,
#'     which is part of the University of Tasmania, and especially to
#'     Dr Craig Mundy, leader of the Abalone Group, for permission to use
#'     this data collected in February 1995.
#'
#' @examples
#'  data(tasab)
#'  head(tasab,20)
#'  table(tasab$site,tasab$sex)
NULL


#' @title zone a list of 8 major components from makeequilzone function
#'
#' @description zone is the output from the function makeequilzone and contains
#'     8 R objects including zoneC, zoneD, glb, constants, saudat, product,
#'     ctrl, and zone1. This is a big object of about 1 Mb.
#'
#' @name zone
#'
#' @docType data
#'
#' @examples
#'  data(zone)
#'  str(zone,max.level=1)
NULL




