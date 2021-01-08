

#' @title product is the productivity curve matrix from doproduction
#'
#' @description product is the productivity curve matrix from
#'     doproduction when the example zone is generated using the
#'     inbuilt datasets ctrl, zone1, and constants. The slowest
#'     part of building the whole is to use the modregC function
#'     to adjust the zoneC and generate the production array. To
#'     save that time in the examples (to avoid time limits on
#'     examples should this package go to CRAN), then this dataset can
#'     be used instead. This is a three dimensional array of
#'     productivity variables.
#'
#' @name product
#'
#' @docType data
#'
#' @section contents:
#' \itemize{
#'   \item harvestrate the initial harvest rates applied
#'   \item productivity variables ExB, MatB, AnnH, Catch, Deplet, RelCE
#'   \item population the index of each population
#' }
#'
#' @examples
#'  data(product)
#'  product[1:20,,1]
NULL

#' @title testzoneC is a zone list made up of 6 equilibrium populations
#'
#' @description testzoneC is a zone list made up of 6 equilibrium
#'     populations. These have been run with a laral dispersal rate of
#'     0.03 so the change from B0 to effB0 is not great, but still
#'     required for an initial equilibrium. This is here to simplify
#'     the internal testing of funcitons that require a completed
#'     zone starting at equilibrium. Its name is to avoid conflict
#'     with any actual use of zoneC. use str(testzoneC, max.level=1)
#'     to see its format. It can be expected to be used with testzoneD
#'
#' @name testzoneC
#'
#' @docType data
#'
#' @section Subjects:
#'  \itemize{
#'    \item testing of functions that require a full zone
#'    \item initial equilibrium
#'  }
#'  @export
#'
#' @examples
#'  data(testzoneC)
#'  data(testzoneD)
#'  data(zone1)
#'  glb <- zone1$globals
#'  r0 <- getvar(testzoneC,"R0")
#'  move <- makemove(glb$numpop,r0,glb$larvdisp)
#'  glb$move <- move
#'  ans <- testequil(testzoneC, testzoneD, glb)
#'  str(testzoneC[[1]])
NULL

#' @title testzoneD is a list of 8 matrices and 2 arrays defining the dynamics of a zone
#'
#' @description testzoneD is a list of 8 matrices and 2 arrays defining
#'     the dynamics of a zone. These have been run with a larval
#'     dispersal rate of 0.03 to achieve an initial equilibrium. This
#'     is here to simplify the internal testing of functions that
#'     require a completed zone starting at equilibrium. Its name is
#'     to avoid conflict with any actual use of zoneD. use
#'     str(testzoneD, max.level=1) to see its format. It can be
#'     expected to be used with testzoneC.
#'
#' @name testzoneD
#'
#' @docType data
#'
#' @section Subjects:
#'  \itemize{
#'    \item testing of functions that require a full zone
#'    \item initial equilibrium
#'  }
#'  @export
#'
#' @examples
#'  data(testzoneC)
#'  data(testzoneD)
#'  data(zone1)
#'  glb <- zone1$globals
#'  r0 <- getvar(testzoneC,"R0")
#'  move <- makemove(glb$numpop,r0,glb$larvdisp)
#'  glb$move <- move
#'  ans <- testequil(testzoneC, testzoneD, glb)
#'  str(testzoneD)
NULL

#' @title getgrad1 calculates the one year gradient score
#'
#' @description getgrad1 calculates the one year gradient score for a
#'    vector of CPUE. The equation used is [CE(y) / CE(y-1)] - 1,
#'    which provides the annual proportional change in CPUE.
#'
#' @param vectce vector of cpue for a given spatial scale
#'
#' @return a vector of gradient1 for each year, starting from year 2
#' @export
#'
#' @examples
#' data(blockE13)
#' nyr <- length(blockE13$year)
#' grad1 <- getgrad1(blockE13$cpue)
#' score1 <- getscore(grad1)
#' cbind(blockE13$year[2:nyr],grad1,score1)
getgrad1 <- function(vectce) {
  nyr <- length(vectce)
  grad1 <- (vectce[2:nyr]/vectce[1:(nyr-1)])-1
  return(grad1)
} # end of getgrad1




#' @title getgrad4 applies a linear regression in steps of wid to input
#'
#' @description getgrad4 takes an input vector of cpue and, in chunks
#'     of length wid, converts them to proportional changes by
#'     dividing through by the first value of the short series, then
#'     applies a linear regression keeping only the gradient.
#'
#' @param vectce the input vector of cpue
#' @param wid the number of years of the input data to use in each
#'     regression
#'
#' @return a vector of length (wid-1) shorter than the input vector
#' @export
#'
#' @examples
#' x <- c(0.0169,0.1953,0.1102,0.1511,-0.0403,-0.0247,-0.0255,-0.1089,
#'        -0.1458,-0.2082,0.0289,-0.0267)
#' grad4 <- getgrad4(x,wid=4)
#' grad3 <- getgrad4(x,wid=3)
#' cbind(c(NA,grad4),grad3)
getgrad4 <- function(vectce,wid=4) { # vectce=ab$cpue; wid=4
  nyr <- length(vectce)
  inc <- wid-1
  num <- nyr-wid+1
  x <- 1:wid
  grad4 <- numeric(num)
  for (i in 1:num) {
    propce <- vectce[i:(i+wid-1)]/vectce[i]
    grad4[i] <- coef(lm(propce ~ x))[2]
  }
  return(grad4)
} # end of getgrad4
