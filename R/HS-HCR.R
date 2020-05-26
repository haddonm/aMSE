
#' @title getaav calculates annual absolute variation in catch
#'
#' @description getaav calculates the annual absolute change in catch
#'     for an input vector of catches, which could be across a series
#'     of years or even across different spatial units for a single
#'     year (an unusual use).
#'     The equation used is aav = 100 x sum(|Ct - Ct-1|)/(sum(Ct).
#'
#' @param invect a vector of catches
#'
#' @return a single scalar value the AAV of the input catches
#' @export
#'
#' @examples
#'   catch <- c(1,2,3,4,5,4,3,2,1)
#'   getaav(catch)  # should equal 0.32
getaav <- function(invect) { # invect=x
  nyr <- length(invect)
  totC <- sum(invect,na.rm=T)
  aac <- sum(abs(invect[2:nyr] - invect[1:(nyr-1)]))
  aav <- 0.0
  if (totC > 0.0) aav <- aac/totC
  return(aav)
} # end of getaav

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
#' score1 <- getscore1(grad1)
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

#' @title getscore calculates the scores for the grad1 PM
#'
#' @description getscore calculates the scores for the grad1 and grad4
#'     performance measures. It does this by re-scaling the range of
#'     the PM values and then fitting separate linear regressions to
#'     the values above and below zero. These enable it to calculate
#'     the predicted scores.
#'
#' @param grad14 the values derived from the function getgrad1 or
#'     getgrad4
#'
#' @return a vector of scores to be included in the MCDA
#' @export
#'
#' @examples
#' data(blockE13)
#' nyr <- length(blockE13$year)
#' grad1 <- getgrad1(blockE13$cpue)
#' score1 <- getscore(grad1)
#' cbind(blockE13$year[2:nyr],grad1,score1)
getscore <- function(grad14) {
  bounds <- round((range(grad14,na.rm=TRUE) * 1.1),2)
  low <- seq(bounds[1],0.0,length=6)
  high <- seq(0.0,bounds[2],length=6)
  xax <- c(low[1:5],high)
  model1 <- lm(0:5 ~ xax[1:6])
  vars <- coef(model1)
  score <- grad14 # just to get a vector of the correct length
  pickl0 <- which(grad14 <= 0)
  score[pickl0] <- grad14[pickl0]*vars[2] + vars[1]
  model2 <- lm(5:10 ~ xax[6:11])
  vars2 <- coef(model2)
  pickg0 <- which(grad14 >= 0)
  score[pickg0] <- grad14[pickg0]*vars2[2] + vars2[1]
  return(score)
} # end of getscore
