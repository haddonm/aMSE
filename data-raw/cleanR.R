

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




























