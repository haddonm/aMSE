

#' @title condDat is an abalone data-set for conditioning the MSE
#'
#' @description condDat is an abalone data-set for blacklip abalone
#'     (\emph{Haliotis rubra}) that can be used to condition the
#'     operating model as an example when running the aMSE function
#'     examples. It describes a region containing two hypothetical
#'     Tasmanian blocks (the level of SMU) with a total of six
#'     populations. If this were a data-set in a file it would require
#'     the readdatafile function to read the file. The function
#'     datafiletemplate will generate a template of the required
#'     format to be read by readdatafile, which can be edited
#'     to suit a given fishery.
#'
#' @name condDat
#'
#' @docType data
#'
#' @format A list of hypothetical conditioning data for the MSE
#' \describe{
#'   \item{constants}{the PDFs used to describe the biological parameters}
#'   \item{Lt}{the length at tagging}
#'   \item{Dt}{The time interval between tagging and recapture, in this
#'       instance they are all listed as 1 year}
#'   \item{DL}{the growth increment in mm}
#' }
#'
#' @examples
#'  data(condDat)
#'  str(condDat)
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

