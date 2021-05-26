
#' @importFrom grDevices palette dev.off
#' @importFrom graphics grid hist lines mtext par plot points title
#' @importFrom graphics legend abline text
#' @importFrom stats dnorm pnorm rlnorm rnorm quantile sd lm coef median
#' @importFrom utils object.size browseURL packageDescription read.csv
#' @importFrom utils write.table head tail
#' @importFrom rutilsMH getmax getmin setpalette plotprep removeEmpty
#' @importFrom rutilsMH splitDate which.closest parset getparplots
#' @importFrom makehtml addplot addtable dirExists filenametopath
#' @importFrom makehtml getextension htmltable logfilename make_html pathend
#' @importFrom makehtml pathtype setuphtml write_css write_head
NULL

#' @title aMSE functions for Conditioning and Running an Invertebrate MSE
#'
#' @description The aMSE package provides functions to facilitate
#'     the conditioning and running of an invertebrate Management
#'     Strategy Evaluation system. The Operating model dynamics are
#'     based around using size-based dynamics rather than age-based
#'     dynamics (because most invertebrates, like abalone, are difficult
#'     to age accurately and consistently). The growth dynamics are
#'     described using the inverse logistic curve (see Haddon et al.
#'     2008). Earlier versions of this invertebrate MSE are described
#'     on A development version is #useDynLib aMSE, registration = TRUE
#'     #importom Rcpp evalCpp
#'     available on GitHub at github.com/haddonm/aMSE.
#'
#' @references Haddon, M., Mundy, C., and D. Tarbath (2008) Using an
#'     inverse-logistic model to describe growth increments of blacklip
#'     abalone (\emph{Haliotis rubra}) in Tasmania. \emph{Fishery Bulletin}
#'     106: 58-71.
#'
#'     Haddon, M., Mayfield, S., Helidoniotis, F., Chick, R. and C.
#'     Mundy (2013) \emph{Identification and Evaluation of Performance
#'     Indicators for Abalone Fisheries}. FRDC Final Report 2007/020.
#'     CSIRO Oceans and Atmosphere and Fisheries Research Development
#'     Corporation. 295 p.
#'
#'     Haddon, M. and C. Mundy (2016) \emph{Testing abalone empirical
#'     harvest strategies for setting TACs and associated LMLs, which
#'     include the use of novel spatially explicit performance measures}.
#'     FRDC Final Report 2011/028. CSIRO Oceans and Atmosphere and
#'     Fisheries Research Development Corporation. Hobart 182 p
#'
#' @section Data sets:
#' \describe{
#'   \item{constants}{ is conditioning data for 6 populations}
#'   \item{ctrl}{control file for a 2 MSU 6 population MSE run}
#'   \item{midg}{an abalone tagging data-set from the Actaeons}
#'   \item{product}{productivity matrix from a 2 MSU 6 population
#'       example}
#'   \item{zone1}{the constants common to a zone}
#'   \item{tasab}{Abalone maturity data from two sites on the
#'       south-west of Tasmania.}
#'   \item{taszoneC}{a zone list made up of 6 equilibrium populations}
#'   \item{taszoneD}{a list of 8 matrices and 2 arrays defining the
#'       dynamics of a zone}
#' }
#' @docType package
#' @name aMSE
NULL




