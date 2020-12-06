

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

