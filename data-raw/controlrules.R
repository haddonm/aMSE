



## BLOCK BASED HCR
## scores need to be rbinded rbind(gradScore,targScore,rate1Score, etc)
## in the same order as the weights for the different HCR
## scores <- rbind(blkGrad,blkTarg); weight=mcdaWts

#' @title blockMCDA generates the TAC adjustment from the scores and weights
#'
#' @description blockMCDA - generates the TAC adjustment from the scores and
#'     weights for however many performance measures yolu have HCR for. This
#'     will need some careful alterations to enable the generation of new
#'     HCR that can be added easily to the dynamics - or old ones removed.
#' @param scores - a vector of scores for each perforamcne measure
#' @param weight - a vector of relative weigths to assign to each of the
#'    scores. Must sum to one or the run is stopped with an error.
#' @param schedule - identifies the TAC adjustment schedule 1 = TACadj
#'    schedule 2 = TACadj2, both of which are defined in the ctrlfile
#' @return a vector of TAC multipliers, one for each block being simulated.
#' @export
#'
#' @examples
#' \dontrun{
#'  TACadj <- c(0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.25)
#'  TACadj2 <- c(0.25,0.8,0.85,0.9,1,1,1.05,1.1,1.15,1.2,1.2)
#'  scores <- rbind(c(4,3.5,3.5,4),c(2.7,2.8,3.2,3.7),c(5,4,4,4.8))
#'  print(scores)
#'  mcdaWts <- c(0.4,0.5,0.1)
#'  blockMCDA(scores,mcdaWts,schedule=1)
#'  blockMCDA(scores,mcdaWts,schedule=2)
#' }
blockMCDA <- function(scores, weight, schedule=1) {
  if (sum(weight) != 1.0) stop("FATAL ERROR: Invalid weights in blockMCDA")
  nblock <- length(scores[1,])
  pickIndex <- numeric(nblock) # to identify the TACadj value to be used
  lookupindex <- 1:11
  finalScore <- as.numeric(weight %*% scores)
  if (schedule == 1) {
    pickIndex <- trunc(finalScore) + 1
    multiplyTAC <- TACadj[pickIndex]
  } else {
    for (blk in 1:nblock) pickIndex[blk] <- which.closest(finalScore[blk],lookupindex)
    multiplyTAC <- TACadj2[pickIndex]
  }
  return(multiplyTAC)   # a vector of length nblock
}  # end of MCDA





#' @title draftHCR generates a brief outline for a new HCR
#'
#' @description draftHCR generates a brief outline for a new HCR which
#'    is designed to calculate the scores for a new performance measure
#'
#' @return returns nothing but prints an outline draftHCR function to
#'    the screen
#' @export
#'
#' @examples
#' # library(AbMSE)
#' draftHCR()
draftHCR <- function() {
  cat("draftHCR <- function(pmbyBlock,scorescale) { \n")
  cat("  #pmbyBlock is a matrix nrow=Nyrs,ncol=(nblock * reps) \n")
  cat("  score <- numeric((nblock * reps)) \n")
  cat("  #Calculate the score using the performance measure \n")
  cat("  #the two lines here are just in case but should not be reached \n")
  cat("  #so adjust your calculations appropriately - perhaps using \n")
  cat("  #scorescale \n")
  cat("  score[score > 10.0] <- 10.0 \n")
  cat("  score[score < 0.0] <- 0.0 \n")
  cat("  return(score) \n")
  cat("} # end of draftHCR \n")
} # end of draftHCR



#' @title gradblockHCR - calculate the gradient4 performance measure
#'
#' @description gradblockHCR - calculate the gradient4 performance
#'     measure for each block using cpueBlock
#'
#' @param incpueBlock - the matrix of cpue by block by year
#' @param maxGradient - the scale of gradients used.
#'
#' @return a vector of scores for Grad4
#' @export gradblockHCR
#'
#' @examples
#' \dontrun{
#'  print("Need to use a dataset for an example")
#' }
gradblockHCR <- function(incpueBlock, maxGradient=maxGrad4) {
  score <- rep(5,nblock)
  tmp <- rep(5,nblock)
  cePeriod <- length(incpueBlock[,1]) # includes implementE
  yrs <- seq(1,cePeriod,1)
  percCE <- apply(incpueBlock,2,(function(x) x/x[1]))
  for (blk in 1:nblock) {
    model <- lm(percCE[,blk]~yrs)
    grad <- model$coeff[2]
    trial <- (5/maxGradient) * grad + 5  # A REAL NUMBER
    if (trial > 10) trial <- 10
    if (trial < 0) trial <- 0
    score[blk] <- trial
  }
  return(score)
} # end of gradblockHCR


#' @title targceHCR - calculates SMU scores for the targetCE
#'
#' @description targceHCR - calculates SMU scores for the targetCE
#'
#' @param incpueSMU the matrix of cpue by SMU by year
#' @param targetCE - a vector of the targetCE for each of the SMUs
#' @param modifyTarg - the constant that sets the range of CPUE in the
#'    scoring function
#'
#' @return a vector of scores relating to the targetCE PM
#' @export
#'
#' @examples
#'  print("Need to use a dataset for an example")
targceHCR <- function(incpueBlock, targetCE, modifyTarg=deltaCE) {
  delCE <- 5.0/modifyTarg
  score <- (delCE * incpueBlock) + 5.0 - (delCE * targetCE)
  score[score > 10.0] <- 10.0
  score[score < 0.0] <- 0.0
  return(score)  # not yet an integer
}  # end of targblockHCR



