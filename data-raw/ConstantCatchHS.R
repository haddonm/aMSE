
#' @title makeouthcr creates an object to contain the replicate HCR outputs
#'
#' @description makeouthcr creates an objects required to hold the cumulative
#'     values for the grad4 pm. Otherwise the full run from 1992 - endyr must be
#'     recalculated each time. In this case it is a 3-D array
#'     set up to hold the values used on each SAUs previous years
#'     aspirational catch. Copies of the output array could also be used to
#'     store the scores developed for the different components of the HCR.
#'
#' @param glb the globals object for the MSE
#' @param hsargs the harvest strategy arguments, in this case a vector of
#'     aspirational catches for each sau
#'
#' @return in this case it merely returns the vector of catches in hsargs
#' @export
#'
#' @examples
#' makeoutconst(glb=NULL,hsargs=c(10,20,30,40,50,60,70,80))
makeoutconst <- function(glb,hsargs) { # glb=glb; hsargs=hsargs
  return(hsargs)
} # end of makeoutconst

#' @title consthcr conducts the MCDA and returns the TAC multiplier and details
#'
#' @description consthcr conducts the MCDA on the basis of a vector of cpue and
#'      other details prescribed in the function's arguments. It returns the
#'      TAC multiplier, the combined score, and all the details of the
#'      calculation. The hcr is a vector of 1:10 values where the cell index
#'      is used to allocate each combined score to a multiplier. The score is
#'      rounded up to the nearest integer and that is the index within hcr.
#'      Thus a score of <=1 points to the first cell, >1 and <=2 points to the
#'      second cell, and so on, up to a score between >9 and <=10 which points
#'      to the last cell. All the arguments used within the mcda are brought in,
#'      and hence are flexible, inside a list names hcrargs. This list includes
#'      the 'wid', the number of years to use in grad4, default=4; the
#'      'targqnt', the quantile for the cpue target level, default=0.55, the
#'      'pmwts', performance measure weights, default = 0.65, 0.25, 0.1 for the
#'      target, grad4, and grad1 PMs respectively, so the order matters; the
#'      'hcr', the harvest control rule scales that transform the combined score
#'      into a TAC multiplier. A vector of 1 - 10 where each cell index
#'      represents the upper limit of the the combined score.
#'
#' @param indat a list of the array of cpue across all years that forms the
#'     basis of the assessment, the yearnames for the extent of years making up
#'     the cpue data, the aspirational catches from the previous year (that
#'     will be multiplied to give the coming year's aspirational catches), and
#'     a NULL fis and NaS value
#' @param hsargs a vector setting catches by SAU. They should sum to the TAC
#'
#' @return a list of the acatch and TAC
#' @export
#'
#' @examples
#' print("wait on data and time")
#' # hcrdata,hsargs,saunames=glb$saunames
#' # indat=constdata; hsargs=hsargs;saunames=glb$saunames
consthcr <- function(indat,hsargs,saunames) {
  acatch <- hsargs
  TAC <- sum(acatch,na.rm=TRUE)
  out <- list(acatch=acatch,TAC=TAC)
  return(out)
} # end of mcdahcr



#' @title constCPUE prepares the projected cpue data ready for the constHCR
#'
#' @describe constCPUE prepares the projected cpue data and the yearnames from
#'    zoneDP ready for their use within the constHCR. The yearnames identifies
#'    the extent of the years of data being included in the calculation of the
#'    cpue performance measures.
#'
#' @param cesau the projected cpue by SAU for a particular year and iteration
#'     from zoneDP
#' @param iter the specific iteration being considered
#'
#' @return a list of the array of cpue and the yearnames vector
#' @export
#'
#' @examples
#' print("wait on suitable data sets")
#' # cesau=zoneDP$cesau; iter=1; year=59
constCPUE <- function(cesau,year) {
  yrnames <- as.numeric(rownames(cesau))
  begin <- which(yrnames == 1990)
  arrce=cesau[begin:(year-1),]
  return(list(arrce=arrce,yearnames=yrnames[begin:(year-1)]))
} # end of constCPUE

#' @title constdata processes the projections to produce the data for the HCR
#'
#' @description constdata takes the projection results from zoneDP, and histCE
#'     from 'otherdata', and appropriately samples them ready for the HCR
#'
#' @param constCPUE a function that generates the CPUE statistics
#' @param constFIS a function that generates the FIS statistics
#' @param constNaS a function that generates the Numbers-at-size samples
#' @param sauCPUE the cpue by sau object from the dynamics object from the
#'     replicate model runs
#' @param sauacatch aspirational catches for each SAU from year-1
#' @param year the specific year being considered
#' @param iter the specific iteration among the replicates being considered
#' @param yrnames the yearnames from glb use to name array dimensions
#'
#' @return a list of arrce, yearnames, acatch, fis, and nas
#' @export
#'
#' @examples
#' print("wait on suitable data sets")
#' # constCPUE=sampleCE;constFIS=sampleFIS;constNAS=sampleNaS;sauCPUE=zoneDP$cesau;year=year;iter=iter;
#' # sauacatch=zoneDP$acatch
constdata <- function(constCPUE, constFIS, constNaS, sauCPUE, sauacatch, year, iter) {
  return(NULL)
} # end of constdata


#' @title constFIS calculates required FIS data
#'
#' @description constFIS no FIS data are currently required by the Tasmanian
#'     HCR so this function merely returns NULL
#'
#' @param x a dummy variable not used
#'
#' @return NULL
#' @export
#'
#' @examples
#' x <- constFIS(NA)
#' print(x)
constFIS <- function(x) {
  return(NULL)
}

#' @title constNas calculates required Numbers-at-Size data
#'
#' @description constNas no numbers-at-size data are currently required by the
#'     Tasmanian HCR so this function merely returns NULL
#'
#' @param x a dummy variable not used
#'
#' @return NULL
#' @export
#'
#' @examples
#' x <- constNaS(NA)
#' print(x)
constNaS <- function(x) {
  return(NULL)
}

