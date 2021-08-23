
# mult the multiplier on the performance measure bounds to expand them both
#      upwards and downwards. default value = 1.1 = 10 percent increase.
# wid  the number of years over which to calculate the gradient, default value
#      = 4, meaning four years
# targqnts what quantile of the distribution of cpue to use as the target,
#      default value = 0.55
# pmwts what weights to give to each of the performance measures. Their order
#      is targetCE, grad4, and grad1 with default values = c(0.65, 0.25,0.1)
# hcr  is used to translate the overall score between 0 - 10, into a
#      multiplier for the previous aspirational catch


hsargs <- list(mult=0.05,
               wid = 4,
               targqnt = 0.55,
               maxtarg = c(30,30,30,30,30,30),
               pmwts = c(0.65,0.25,0.1),
               hcr = c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2),
               startCE = 1990
               )


#' @title calcexpectpopC generates population catches from TAC or acatch
#'
#' @description calcexpectpopC converts the TAC into expected sau catches,
#'     which are then translated into population level catches while
#'     introducing management implementation error. Here this
#'     error is implemented as Log-Normal errors on diver intuitions
#'     concerning the relative abundance in each population. The error is
#'     imposed separately on the populations in each SAU. The option of TAC for
#'     a zone or aspirational catches for SAU is given to allow for different
#'     processes between jurisdictions. If only a TAC is produced by an HCR then
#'     calcexpectpopC needs to be able to generate catches by SAU prior to
#'     subdividing those SAU catches among populations.
#'
#' @param hcrout in Tasmania this contains the aspirational catches by SAU
#'     along with many other details.
#' @param exb the exploitable biomass at the end of the previous year. In the
#'     first year of the projections this would be the last year of the
#'     conditioning.
#' @param sauindex the SAU index for each population
#' @param sigmab the Log-Normal standard deviation of implementation error. The
#'     default value = 1e-08, which effectively means no errors.
#'
#' @return a vector of population catches for the year to be imposed after the
#'     estimation of the mid-year exploitable biomass (midyexpB).
#' @export
#'
#' @examples
#' print("wait on suitable internal data sets")
calcexpectpopC <- function(hcrout,exb,sauindex,sigmab=0.1) {
  # exb=zoneDD$exploitB[endyr,]; suaindex=glb$sauindex; sigmaB=0.1;
  acatch <- hcrout$acatch
  npop <- length(exb)
  nSAU <- length(acatch)
  sauexb <- tapply(exb,sauindex,sum,na.rm=TRUE)# * exp(rnorm(nSAU,mean=0,sd=sigmab))
  popC <- acatch[sauindex] * (exb/sauexb[sauindex])
  return(popC)
} # end of calcexpectpopC

#' @title catchbysau calculates the catch by sau with error for Tasmania
#'
#' @description catchbysau calculates the catch by sau with error for the
#'     Tasmanian HS. This takes the sum of the aspirational catches as the TAC,
#'     which is multiplied by the proportion of exploitable biomass for each SAU,
#'     which has log-Normal error included with sd=sigmaB.
#'
#' @param inexpB a particular years' exploitable biomass by population
#' @param sauindex a vector containing the sau index number for each population
#' @param TAC the sum of the aspirational catches for a given year
#' @param sigmab the sd of the log-normal errors included in the estimates of
#'     the exploitable biomass by SAU
#'
#' @return a vector of the catches by sau
#' @export
#'
#' @examples
#' print("wait on suitable internal data")
catchbysau <- function(inexpB,sauindex,TAC,sigmab) {
  #  iter=1;year=2
  # inexpB=zoneDP$exploitB[(year - 1),,iter];sauindex=sauindex;TAC=origTAC[iter]; sigmab=0.05
  totexB <- sum(inexpB,na.rm=TRUE)
  sauexpB <- tapply(inexpB,sauindex,sum,na.rm=TRUE)
  nSAU <- length(sauexpB)
  npop <- length(inexpB)
  sauexpB <- sauexpB * exp(rnorm(nSAU,mean=0,sd=sigmab))
  totsau <- sum(sauexpB,na.rm=TRUE)
  sauexpB <- sauexpB * (totexB/totsau)
  catbysau <- TAC * sauexpB/sum(sauexpB)
  return(catbysau)
} # end of catchbysau

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
#'
#' @return in this case a list three dimensional array designed to hold
#' @export
#'
#' @examples
makeouthcr <- function(glb) { # glb=glb; hsargs=hsargs
  begin <- hsargs
  grad4val <- array(0,dim=c(glb$pyrs,glb$nSAU,glb$reps),
                           dimnames=list(glb$pyrnames,glb$saunames,1:glb$reps))
  hcrin <- list(grad4val=grad4val)
  return(hcrin)
} # end of makeouthcr

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


#' @title targscore generates the HCR score for the target PM
#'
#' @description targscore takes in a vector of cpue x years and defines the targ
#'     from the complete series, and the limit reference point. This differs
#'     from the two getgradone and getgradwid in needing multiple years at once.
#'     One meta-rule is that the target cpue should not rise above
#'
#' @param vectce the vector of cpue to be assessed
#' @param qnt the quantile of the input vector selected as the target,
#'     default = 0.55
#' @param mult the multiplier on the bounds to expand them upwards and
#'     downwards. default value = 1.1 = 10 percent increase
#' @param maxtarg a meta-rule that sets an upper limit on the target cpue, the
#'     default=150kg/hr
#'
#' @return a list of the final year's score, the internals to the calculations,
#'     and he target and limit reference points
#' @export
#'
#' @examples
#' print("wait on suitable data")
#' # vectce=cpue[30:88,1]; qnt=0.55;mult=0.1; maxtarg=150.0
targscore <- function(vectce,qnt=0.55,mult=0.1,maxtarg=150.0) { #
  nyr <- length(vectce)
  targ <- quantile(vectce,probs=c(qnt),na.rm=TRUE)
  if (targ > maxtarg) targ <- maxtarg
  bounds <- round(extendrange(vectce,range(vectce,na.rm=TRUE),f=mult),3)
  low <- seq(bounds[1],targ,length=6)
  high <- seq(targ,bounds[2],length=6)
  xax <- c(low[1:5],high)
  vars <- getlmcoef(0:5,xax[1:6])
  score <- numeric(length(vectce))
  pickl <- which(vectce <= targ)
  score[pickl] <- vectce[pickl]*vars[2] + vars[1]
  vars2 <- getlmcoef(5:10,xax[6:11])
  pickh <- which(vectce >= targ)
  score[pickh] <- vectce[pickh]*vars2[2] + vars2[1]
  rp <- c(bounds[1],targ,bounds[2]); names(rp)=c("lower","target","upper")
  result <- tail(score,1)
  return(list(result=result,scores=score,rp=rp))
} # end of targscore


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

