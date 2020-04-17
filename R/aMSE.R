
#' @importFrom grDevices palette
#' @importFrom graphics grid hist lines mtext par plot points title
#' @importFrom stats dnorm pnorm rlnorm rnorm quantile sd
#' @importFrom utils object.size browseURL packageDescription read.csv
#' @importFrom rutilsMH getmax getmin setpalette plotprep removeEmpty
#' @importFrom rutilsMH splitDate which.closest parset
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
#'     on A development version is
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
#'   \item{midg}{a tagging data-set from abalone taken from the Middle
#'       Ground in the Actaeons.}
#'   \item{tasab}{Abalone maturity data from two sites on the
#'       south-west of Tasmania.}
#'   \item{condDat}{a two block, six population data set for conditioning
#'       the operating model to be used in the function examples.}
#' }
#' @docType package
#' @name aMSE
NULL




