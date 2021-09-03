
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
  # hcrout= list(acatch=acatch)
  acatch <- hcrout$acatch
  TAC <- sum(acatch,na.rm=TRUE)  # means TAC estimation irrelevant in TAS
  totexb <- sum(exb,na.rm=TRUE)
  npop <- length(exb)
  nSAU <- length(acatch)
  sauexb <- tapply(exb,sauindex,sum,na.rm=TRUE) * exp(rnorm(nSAU,mean=0,sd=sigmab))
  sauexb <- sauexb * (totexb/sum(sauexb,na.rm=TRUE))
  popC <- acatch[sauindex] * (exb/sauexb[sauindex])
  popC <- popC * TAC/sum(popC)
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


#' @title getgrad1 calculates the one year gradient score for all years of cpue
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
  grad1 <- c(NA,(vectce[2:nyr]/vectce[1:(nyr-1)])-1)
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
  grad4 <- numeric(nyr)
  grad4[1:inc] <- NA
  for (i in 1:num) { # i=1
    propce <- vectce[i:(i+wid-1)]/vectce[i]
    grad4[i+inc] <- getlmcoef(propce,1:wid)[2]
    #grad4[i+inc] <- coef(lm(propce ~ x))[2]
  }
  return(grad4)
} # end of getgrad4


#' @title getgradone calculates the one year gradient score
#'
#' @description getgradone calculates the one year gradient score for a
#'    vector of CPUE. The equation used is [CE(y)/CE(y-1)] - 1,
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
#' grad1 <- getgradone(blockE13$cpue,yr = 10)
#' print(grad1)
getgradone <- function(vectce) {
  yr <- length(vectce)
  grad1 <- (vectce[yr]/vectce[yr-1])-1
  return(grad1)
} # end of getgradone

#' @title getgradwid returns a single gradient across a cpue vector
#'
#' @description getgradwid calculates the regression gradient of a time-series
#'     of cpue values of length wid, using the analytical formulae rather than
#'     lm (for speed). It always uses the (lastyear-wid-1):lastyear to make
#'     a single gradient estimate.
#'
#' @param vectce the vector of cpue. it must be at least wid yr long
#' @param wid the number of years across which the gradient is to be estimated
#'
#' @return a single gradient
#' @export
#'
#' @examples
#' vectcpue <- c(0.867,0.885,0.877,1.008,1.169,1.231,1.226,1.237,1.238,1.156)
#' getgradwid(vectcpue,wid=4)   #  -0.01704731
getgradwid <- function(vectce,wid=4) {
  lastyr <- length(vectce)
  inc <- wid-1
  yrrge <- (lastyr-inc):lastyr
  x <- 1:wid
  xres <- (x-mean(x,na.rm=TRUE))
  propce <- vectce[yrrge]/vectce[lastyr-inc]
  yres <- (propce - mean(propce,na.rm=TRUE))
  x2 <- sum(xres^2)
  xy <- sum(xres * yres)
  return(xy/x2)
} #end getgradwid

#' @title getlmcoef is a greatly simplified way of gasining regression coefs
#'
#' @description getlmcoef is a simplified replacement for the lm function that
#'     is much faster than lm but only returns the linear regression
#'     coefficients. This is limited to a simple y ~ x model. When in use be
#'     sure to give it the values of y required.
#'
#' @param y the dependent variable
#' @param x the independent variable
#'
#' @return a vector of the intercept and gradient, in that order.
#' @export
#'
#' @examples
#' y <- c(1.226,1.237,1.238,1.156)
#' x <- 1:4
#' yp <- y/y[1]  # to make the cpue relative to the first in the series
#' getlmcoef(yp,x)
#' getgradwid(yp,4)
getlmcoef <- function(y,x) {
  mx <- mean(x,na.rm=TRUE)
  xres <- (x-mx)
  my <- mean(y,na.rm=TRUE)
  yres <- (y - my)
  x2 <- sum(xres^2)
  xy <- sum(xres * yres)
  grad <- xy/x2
  inter <- my - grad*mx
  return(c(inter,grad))
} # end of getlmcoef

#' @title getscore calculates the scores for the grad1 and grad4 PMs
#'
#' @description getscore calculates the scores for the grad1 and grad4
#'     performance measures. It does this by re-scaling the range of
#'     the PM values and then fitting separate linear regressions to
#'     the values above and below zero. These enable it to calculate
#'     the predicted scores. This function now contain a filter that
#'     allows for gaps within the catch and cpue data. These generated
#'     NAs, NaNs, Infs, and -1 values, all of which lead to aberrant
#'     behaviour. This means that the HS being used is applied to all
#'     data series, even those with large gaps, which may not be what
#'     is most desirable (but at least now this does'nt fail).
#'
#' @param pm the raw performance measure values derived from the functions
#'     getgradone,  getgradwid, or targscore
#' @param mult the multiplier on the bounds to expand them upwards and
#'     downwards. default value = 1.1 = 10 percent increase; from hsargs$mult
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
getscore <- function(pm,mult=0.1) { # pm=grad1val[,sau]; mult=hsargs$mult
  pickb <- which((pm == (-1)) | (is.nan(pm)) | (pm == Inf) | (is.na(pm)))
  bounds <- round(extendrange(pm[-pickb],range(pm[-pickb],na.rm=TRUE),f=mult),2)
  low <- seq(bounds[1],0.0,length=6)
  high <- seq(0.0,bounds[2],length=6)
  xax <- c(low[1:5],high)
  vars <- getlmcoef(0:5,xax[1:6])
  score <- numeric(length(pm)) # just to get a vector of the correct length
  pickl0 <- which(pm <= 0)
  score[pickl0] <- pm[pickl0]*vars[2] + vars[1]
  vars2 <- getlmcoef(5:10,xax[6:11])
  pickg0 <- which(pm >= 0)
  score[pickg0] <- pm[pickg0]*vars2[2] + vars2[1]
  return(score)
} # end of getscore


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

#' @title mcdahcr conducts the MCDA and returns the TAC multiplier and details
#'
#' @description mcdahcr conducts the MCDA on the basis of a vector of cpue and
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
#' @param hsargs the arguments used within the Harvest strategy's HCR. See the
#'     description for details.
#' @param yearnames the vector of years used as names for the gradients
#' @param saunames a vecotr of the SAU names, again for labelling matrix edges
#' @param acatches the aspirational catches for the previous year which are to
#'     be multiplied by the multTAC to give the next year's aspirational
#'     catches. They are summed to give the next year's TAC
#'
#' @return a list of the acatch, TAC, TAC multiplier, the score, all the
#'     details and a matrix of the reference points
#' @export
#'
#' @examples
#' print("wait on data and time")
#' # hcrdata,hsargs,saunames=glb$saunames
#' # indat=hcrdata; hsargs=hsargs;saunames=glb$saunames
mcdahcr <- function(indat,hsargs,saunames) {
  arrce <- indat$arrce
  nsau <- ncol(arrce)
  yrce <- nrow(arrce)
  pmwts <- hsargs$pmwts
  yearnames <- indat$yearnames
  acatches <- indat$acatches
  # define storage matrices
  grad1val <- matrix(0,nrow=yrce,ncol=nsau,dimnames=list(yearnames,saunames))
  grad4val <- targval <- score1 <- score4 <- scoret <- scoretot <- grad1val
  multTAC <- grad1val
  refpts <- matrix(0,nrow=nsau,ncol=3,dimnames=list(saunames,c("low","trp","high")))
  for (sau in 1:nsau) {  #  sau=6
    pickce <- which(!is.na(arrce[,sau]))
    tmp <- getgrad1(arrce[pickce,sau])                     # grad1
    nec <- length(tmp)
    if (nec < yrce) tmp <- c(rep(NA,(yrce-nec)),tmp)
    grad1val[,sau] <- tmp
    score1[,sau] <- getscore(grad1val[,sau],mult=hsargs$mult)
    tmp2 <- getgrad4(arrce[pickce,sau],wid=hsargs$wid)      # grad4
    nec2 <- length(tmp2)
    if (nec2 < yrce) tmp2 <- c(rep(NA,(yrce-nec2)),tmp2) # allow for 6 and 13
    grad4val[,sau] <- tmp2
    score4[,sau] <- getscore(grad4val[,sau],mult=hsargs$mult)
    tmp3 <- targscore(arrce[pickce,sau],qnt=hsargs$targqnt,mult=hsargs$mult,
                      maxtarg=hsargs$maxtarg[sau]) #targ
    nec3 <- length(pickce)
    scrs <- tmp3$scores
    if (nec3 < yrce) scrs <-  c(rep(NA,(yrce-nec3)),scrs)
    scoret[,sau] <- scrs
    scoretot[,sau] <- pmwts[1]*scoret[,sau] + pmwts[2]*score4[,sau] +
      pmwts[3]*score1[,sau]
    pick <- floor(scoretot[,sau]) + 1 # add one to get correct hcr[index]
    pick[pick < 1] <- 1
    pick[pick > 10] <- 10
    multTAC[,sau] <- hsargs$hcr[pick] # index implies scores less than index
    refpts[sau,] <- tmp3$rp
  }
  acatch <- acatches * multTAC[yrce,]
  TAC <- sum(acatch,na.rm=TRUE)
  details <- list(grad4=grad4val,score4=score4,grad1=grad1val,score1=score1,
                  scoret=scoret,scoretot=scoretot)
  out <- list(acatch=acatch,TAC=TAC,multTAC=multTAC,details=details,
              refpts=refpts)
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


#' @title tasCPUE prepares the projected cpue data ready for the tasHCR
#'
#' @describe tasCPUE prepares the projected cpue data and the yearnames from
#'    zoneDP ready for their use within the tasHCR. The yearnames identifies
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
tasCPUE <- function(cesau,year) {
  yrnames <- as.numeric(rownames(cesau))
  begin <- which(yrnames == 1990)
  arrce=cesau[begin:(year-1),]
  return(list(arrce=arrce,yearnames=yrnames[begin:(year-1)]))
} # end of tasCPUE

#' @title tasdata processes the projections to produce the data for the HCR
#'
#' @description tasdata takes the projection results from zoneDP, and histCE
#'     from 'otherdata', and appropriately samples them ready for the HCR
#'
#' @param tasCPUE a function that generates the CPUE statistics
#' @param tasFIS a function that generates the FIS statistics
#' @param tasNaS a function that generates the Numbers-at-size samples
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
#' # tasCPUE=sampleCE;tasFIS=sampleFIS;tasNAS=sampleNaS;sauCPUE=zoneDP$cesau;year=year;iter=iter;
#' # sauacatch=zoneDP$acatch
tasdata <- function(tasCPUE, tasFIS, tasNaS, sauCPUE, sauacatch, year, iter) {
  outce <- tasCPUE(sauCPUE[,,iter],year)
  yearnames <- outce$yearnames   # omit empty first year
  arrce <- outce$arrce
  acatches=sauacatch[year-1,,iter]
  fis <- tasFIS(NA)
  nas <- tasNaS(NA)
  ans <- list(arrce=arrce,yearnames=yearnames,acatches=acatches,fis=fis, nas=nas)
  return(ans)
} # end of tasdata


#' @title tasFIS calculates required FIS data
#'
#' @description tasFIS no FIS data are currently required by the Tasmanian
#'     HCR so this function merely returns NULL
#'
#' @param x a dummy variable not used
#'
#' @return NULL
#' @export
#'
#' @examples
#' x <- tasFIS(NA)
#' print(x)
tasFIS <- function(x) {
  return(NULL)
}

#' @title tasNas calculates required Numbers-at-Size data
#'
#' @description tasNas no numbers-at-size data are currently required by the
#'     Tasmanian HCR so this function merely returns NULL
#'
#' @param x a dummy variable not used
#'
#' @return NULL
#' @export
#'
#' @examples
#' x <- tasNaS(NA)
#' print(x)
tasNaS <- function(x) {
  return(NULL)
}

