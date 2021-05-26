
#' @title blockE13 is an abalone data-set for testing performance measures
#'
#' @description blockE13 is a fishery data-set for blacklip abalone
#'     (\emph{Haliotis rubra}) from block 13 in the eastern zone, this is
#'     Tasmania's Block 13. It constitutes three time-series of the same cpue
#'     data in different formats (see below). It is for use when testing the
#'     performance measures within Tasmania's MCDA, although it could be used
#'     for other purposes, such as illustrating the typical linear relationship
#'     between catch and CPUE. Note there will only ever be one less PM value
#'     than there are cpue data in the time-series.
#'
#' @name blockE13
#'
#' @docType data
#'
#' @format A data.frame of abalone fishery data
#' \describe{
#'   \item{year}{the year of fishing}
#'   \item{coeff}{the back-transformed coefficients from a standardization}
#'   \item{scaled}{the coefficients scaled to a mean of one}
#'   \item{cpue}{the coefficients scaled to the nominal geometric mean of the
#'       time series, to place it on the nominal scale}
#'   \item{catch}{the related catch from the eastern parts of block 13}
#' }
#'
#' @section Subjects:
#'  \itemize{
#'    \item performance measure estimation
#'    \item relationship between catch and CPUE
#'  }
#'  @export
#'
#' @source Mundy, C. and J.M. McAllister (2020) Tasmanian Abalone Assessment
#'     2019. IMAS, University of Tasmania.
#'
#' @examples
#' data(blockE13)
#' blockE13
NULL

#' @title lf10 contains three years of length-composition data for block 10
#'
#' @description lf10 is a data.frame of length-composition of commercial catch
#'     from block 10 on the west coast of Tasmania in a longdat format as
#'     derived from the commlf function makelongdat.This makes it useful for
#'     illustrating the use of makewidedat.
#'
#' @name lf10
#'
#' @docType data
#'
#' @format A data.frame of abalone size-composition data
#' \describe{
#'   \item{year}{the year of sampling}
#'   \item{sau}{the block or SAU}
#'   \item{length}{the measured size in mm, with 2mm size classes}
#'   \item{counts}{the counts at each length class}
#'   \item{propcounts}{the proportion of each length class of the total}
#' }
#'
#' @section Subjects:
#'  \itemize{
#'    \item selectivity
#'    \item size-based stock assessment model fitting
#'    \item size-distributions of commercial catches
#'  }
#'  @export
#'
#' @source Thanks to the Institute of Marine and Antarctic Science,
#'     which is part of the University of Tasmania, and especially to
#'     Dr Craig Mundy, leader of the Abalone Group, for permission to use
#'     this data.
#'
#' @examples
#'  data(lf10)
#'  head(lf10,20)
#'  mids <- seq(138,210,2)
#'  answer <- makewidedat(lf10,mids)
#'  answer[1:20,]
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
#'  @export
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

#' @title zone the primary object obtained from the function makeequilzone
#'
#' @description zone contains seven objects, including 5 lists, a matrix, and
#'     an array. This is the
#'
#' @name zone
#'
#' @docType data
#'
#' @format A list of objects plus a matrix and array that make up the initial
#'     equilibrium zone
#' \describe{
#'   \item{zoneC}{a list of the constants for each population}
#'   \item{zoneD}{a list of the dynamic parts of the populations of a zone}
#'   \item{glb}{a list of global constants, containing numpop,nSAU,midpts,
#'       Nclass, Nyrs}
#'   \item{constants}{a matrix of biological properties for each population in
#'       the zone, derived from the datafile}
#'   \item{product}{the productivity array from doproduction}
#'   \item{ctrl}{the list containing control information for the run, including
#'       the datafile for the constants, the reps, the variation to be included
#'       when projecting}
#'   \item{zone1}{a list of objects used in the MSE}
#' }
#'
#' @examples
#'  data(zone)
#'  str(zone,max.level=1)
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
#'  @export
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


