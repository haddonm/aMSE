
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


hsargs <- list(mult=0.1,
               wid = 4,
               targqnt = 0.55,
               pmwts = c(0.65, 0.25,0.1),
               hcr = c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2)
               )


#' @title calibrateHCR uses historical CPUE to calibrate the MCDA
#'
#' @description calibrateHCR uses the condC$histCE to apply the getgrad1,
#'     getgrad4, and targscore functions so as to provide vectors of
#'     performance measure scores to enhance the calibration of the MCDA before
#'     any projections are done. This is not absolutely required but improves
#'     the reality of the process. It also means that the productivity of each
#'     SAU needs to be scaled to the current cpue (see Conditioning the MSE).
#'     This uses getgrad1 and getgrad4 instead of getgradone and getgradwid as
#'     it needs to process whole vectors of cpue per sau rather than used a
#'     row of cpue per sau.
#'
#' @param histCE this is a matrix of nominal scale CPUE from the fishery of
#'     interest. Usually standardized. Obtainable from the projC object but
#'     also the condC object
#' @param saunames the names of each SAU (from zone1$SAUnames)
#' @param hsargs the arguments used within the Harvest strategy's HCR. See the
#'     description for details.
#' @param pyrs the number of years of projection
#'
#' @return a list of the list of the PMs, the targce refpts, and then a list of
#'     all the scores, and the multTAC values.
#' @export
#'
#' @examples
#' print("wait on example data being available")
calibrateHCR <- function(histCE,saunames,hsargs) {
  # histCE=condC$histCE; saunames=zone$zone1$SAUnames; hsargs=hsargs;
  #  sauindex=glb$sauindex;
  #  DEPRECATED
  nSAU <- length(saunames)
  yearCE <- as.numeric(rownames(histCE))
  yrce <- length(yearCE)
  # define storage matrices
  grad1val <- matrix(0,nrow=yrce,ncol=nSAU,dimnames=list(yearCE,saunames))
  grad4val <- targval <- score1 <- score4 <- scoret <- scoretot <- grad1val
  multTAC <- grad1val
  refpts <- matrix(0,nrow=nSAU,ncol=2,dimnames=list(saunames,c("trp","lrp")))
  pmwts <- hsargs$pmwts
  for (sau in 1:nSAU) { # sau=1
    pickce <- which(!is.na(histCE[,sau]))
    tmp <- getgrad1(histCE[pickce,sau])
    nec <- length(tmp)
    if (nec < yrce) tmp <- c(rep(NA,(yrce-nec)),tmp)
    grad1val[,sau] <- tmp
    score1[,sau] <- getscore(grad1val[,sau],mult=hsargs$mult)
    tmp2 <- getgrad4(histCE[pickce,sau],wid=hsargs$wid)
    nec2 <- length(tmp2)
    if (nec2 < yrce) tmp2 <- c(rep(NA,(yrce-nec2)),tmp2)
    grad4val[,sau] <- tmp2
    score4[,sau] <- getscore(grad4val[,sau],mult=hsargs$mult)
    tmp3 <- targscore(histCE[pickce,sau],qnt=hsargs$targqnt,mult=hsargs$mult)
    values <- tmp3$ans[,"targval"]
    scrs <- tmp3$ans[,"targsc"]
    nec3 <- length(values)
    if (nec3 < yrce) {
      values <- c(rep(NA,(yrce-nec3)),values)
      scrs <-  c(rep(NA,(yrce-nec3)),scrs)
    }
    targval[,sau] <- values
    scoret[,sau] <- scrs
    scoretot[,sau] <- pmwts[1]*scoret[,sau] + pmwts[2]*score4[,sau] +
      pmwts[3]*score1[,sau]
    pick <- ceiling(scoretot[,sau])
    pick[pick < 1] <- 1
    pick[pick > 10] <- 10
    multTAC[,sau] <- hsargs$hcr[pick]
    refpts[sau,] <- tmp3$details
  }
  pms <- list(grad1val=grad1val,grad4val=grad4val,targval=targval,
              refpts=refpts)
  scores <- list(score1=score1,score4=score4,scoret=scoret,scoretot=scoretot)
  return(list(pms=pms,scores=scores,yrmultTAC=multTAC))
} # end of calibrateHCR





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
#'    vector of CPUE. The equation used is [CE(y) / CE(y-1)] - 1,
#'    which provides the annual proportional change in CPUE.
#'
#' @param vectce vector of cpue for a given spatial scale
#' @param yr the year for which a score is required
#'
#' @return a vector of gradient1 for each year, starting from year 2
#' @export
#'
#' @examples
#' data(blockE13)
#' nyr <- length(blockE13$year)
#' grad1 <- getgradone(blockE13$cpue,yr = 10)
#' print(grad1)
getgradone <- function(vectce,yr) {
  grad1 <- (vectce[yr]/vectce[yr-1])-1
  return(grad1)
} # end of getgradone

#' @title getgradwid returns a single gradient across a cpue vector
#'
#' @description getgradwid calculates the regression gradient of a time-series
#'     of cpue values of length wid, using the analytical formulae rather than
#'     lm (for speed)
#'
#' @param vectce the vector of cpue. it must be at least yr long
#' @param yr the last yr of the wid long time-series whose gradient is wanted
#' @param wid the number of years across which the gradient is to be estimated
#'
#' @return a single gradient
#' @export
#'
#' @examples
#' vectcpue <- c(0.867,0.885,0.877,1.008,1.169,1.231,1.226,1.237,1.238,1.156)
#' getgradwid(vectcpue,yr=10)   #  -0.01704731
getgradwid <- function(vectce,yr,wid=4) {
  inc <- wid-1
  yrrge <- (yr-inc):yr
  x <- 1:wid
  xres <- (x-mean(x,na.rm=TRUE))
  propce <- vectce[yrrge]/vectce[yr-inc]
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
}

#' @title getscore calculates the scores for the grad1 and grad4 PMs
#'
#' @description getscore calculates the scores for the grad1 and grad4
#'     performance measures. It does this by re-scaling the range of
#'     the PM values and then fitting separate linear regressions to
#'     the values above and below zero. These enable it to calculate
#'     the predicted scores.
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
getscore <- function(pm,mult=0.1) { # pm=scoret$ans[,"targval"];mult=hsargs$mult
  bounds <- round(extendrange(pm,range(pm,na.rm=TRUE),f=mult),2)
  if (bounds[1] > -0.05) bounds[1] <- -0.05
  if (bounds[2] <= 0) bounds[2] <- 0.05
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

#' @title makeprojpm creates arrays to hold the projection PMs and HCR results
#'
#' @description makeprojpm creates the arrays required to hold the three
#'     performance measures, grad1val, grad4val, and targval, used in MCDA
#'     calculations. It also creates arrays to hold the multipliers for the
#'     aspirational catches in each SAU, as well as those to hold the updated
#'     values of the target and limit cpue. These are then all output ready to
#'     be used as an argument in doprojection
#'
#' @param histpms the historical performance measure estimates from
#'     calibrateMCDA
#' @param reps the number of replicate projections
#' @param projyrs the number of years of projection
#'
#' @return a list of arrays for the three PMs, the result of the MCDA, and the
#'     updated limit and target reference points for cpue. Finally the number
#'     of years of historical values used.
#' @export
#'
#' @examples
#' print("Wait on suitable internal zone data")
makeprojpm <- function(histpms,reps,projyrs) {
# glb=glb; histpms=cmcda$pms; reps=ctrl$reps; projyrs=projC$projyrs
  saunames <- colnames(histpms$grad1val)
  nSAU <- length(saunames)
  yearCE <- as.numeric(rownames(histpms$grad1val))
  yrce <- length(yearCE)
  allyrs <- yrce + projyrs
  projname <- (yearCE[yrce]+1):(yearCE[yrce]+projyrs)
  yearnames <- c(yearCE,projname)
  grad1val <- array(0,dim=c(allyrs,nSAU,reps),
                    dimnames=list(yearnames,saunames,1:reps))
  grad4val <- targval <- grad1val
  grad1val[1:yrce,,] <- histpms$grad1val
  grad4val[1:yrce,,] <- histpms$grad4val
  targval[1:yrce,,] <- histpms$targval
  multTAC <- array(0,dim=c(projyrs,nSAU,reps),
                   dimnames=list(projname,saunames,1:reps))
  trpCE <- lrpCE <- multTAC
  ans <- list(grad1val=grad1val,grad4val=grad4val,targval=targval,
              multTAC=multTAC,trpCE=trpCE,lrpCE=lrpCE,histyr=yrce)
  return(ans)
} # end of makeprojpm

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
#' @param arrce the array of cpue across all years that forms the basis of the
#'     assessment
#' @param yr the year for which a score is required
#' @param hsargs the arguments used within the Harvest strategy's HCR. See the
#'     description for details.
#'
#' @return a list of the TAC multiplier, the score, and all the details
#' @export
#'
#' @examples
#' print("wait on data and time")
#' #  arrce=condC$histCE; yr=28;
mcdahcr <- function(arrce,hsargs,yearnames,saunames) {
  # arrce=rbind(condC$histCE,zoneDP$cesau[year-1,,iter]); hsargs=hsargs; year=2
  # oldyrs=as.numeric(rownames(condC$histCE)); yearnames=c(oldyrs,tail(oldyrs,1)+1:(year-1))
  # saunames=sort(unique(glb$SAUpop))
  nsau <- ncol(arrce)
  yrce <- nrow(arrce)
  pmwts <- hsargs$pmwts
  # define storage matrices
  grad1val <- matrix(0,nrow=yrce,ncol=nsau,dimnames=list(yearnames,saunames))
  grad4val <- targval <- score1 <- score4 <- scoret <- scoretot <- grad1val
  multTAC <- grad1val
  refpts <- matrix(0,nrow=nsau,ncol=2,dimnames=list(saunames,c("trp","lrp")))
  for (sau in 1:nsau) {  #  sau=1
    pickce <- which(!is.na(arrce[,sau]))
    tmp <- getgrad1(arrce[pickce,sau])                      # grad1
    nec <- length(tmp)
    if (nec < yrce) tmp <- c(rep(NA,(yrce-nec)),tmp)
    grad1val[,sau] <- tmp
    score1[,sau] <- getscore(grad1val[,sau],mult=hsargs$mult)
    tmp2 <- getgrad4(arrce[pickce,sau],wid=hsargs$wid)      # grad4
    nec2 <- length(tmp2)
    if (nec2 < yrce) tmp2 <- c(rep(NA,(yrce-nec2)),tmp2)
    grad4val[,sau] <- tmp2
    score4[,sau] <- getscore(grad4val[,sau],mult=hsargs$mult)
    tmp3 <- targscore(arrce[pickce,sau],qnt=hsargs$targqnt,mult=hsargs$mult) #targ
    values <- tmp3$ans[,"targval"]
    scrs <- tmp3$ans[,"targsc"]
    nec3 <- length(values)
    if (nec3 < yrce) {
      values <- c(rep(NA,(yrce-nec3)),values)
      scrs <-  c(rep(NA,(yrce-nec3)),scrs)
    }
    targval[,sau] <- values
    scoret[,sau] <- scrs
    scoretot[,sau] <- pmwts[1]*scoret[,sau] + pmwts[2]*score4[,sau] +
      pmwts[3]*score1[,sau]
    pick <- ceiling(scoretot[,sau])
    pick[pick < 1] <- 1
    pick[pick > 10] <- 10
    multTAC[,sau] <- hsargs$hcr[pick]
    refpts[sau,] <- tmp3$details
  }
  details <- list(grad4=grad4val,score4=score4,grad1=grad1val,score1=score1,
                  scoret=scoret)
  out <- list(multTAC=multTAC,details=details,refpts=refpts)
  return(out)
} # end of mcdahcr

#' @title targscore generates the HCR score for the target PM
#'
#' @description targscore takes in a vector of cpue x years and defines the targ
#'     from the complete series, and the limit reference point. This differs
#'     from the two getgradone and getgradwid in needing multiple years at once.
#'
#' @param vectce the vector of cpue to be assessed
#' @param qnt the quantile of the input vector selected as the target,
#'     default = 0.55
#' @param mult the multiplier on the bounds to expand them upwards and
#'     downwards. default value = 1.1 = 10 percent increase
#'
#' @return a list of the final year's score, the internals to the calculations,
#'     and he target and limit reference points
#' @export
#'
#' @examples
#' print("wait on suitable data")
targscore <- function(vectce,qnt=0.55,mult=1.1) { # vectce=arrce[,2];qnt=0.55
  nyr <- length(vectce)
  targ <- quantile(vectce,probs=c(qnt),na.rm=TRUE)
  limrp <- min(vectce,na.rm=TRUE)*0.9
  targval <- vectce - targ
  targsc <- getscore(targval,mult=mult)
  ans <- as.matrix(cbind(vectce,targval,targsc))
  dimnames(ans) <- list(1:length(vectce),c("cpue","targval","targsc"))
  details <- c(targ,limrp); names(details)=c("target","limit")
  result <- tail(targsc,1)
  return(list(result=result,ans=ans,details=details))
} # end of targscore

