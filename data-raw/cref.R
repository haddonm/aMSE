constantrefhcr <- function(indat,hsargs,saunames,curryear=NULL) {
  arrce <- indat$arrce
  nsau <- ncol(arrce)
  yrce <- nrow(arrce)
  pmwts <- hsargs$pmwts
  stablewts <- hsargs$stablewts
  pmwtSwitch <- hsargs$pmwtSwitch
  metRunder <- hsargs$metRunder
  metRover <- hsargs$metRover
  yearnames <- indat$yearnames
  acatches <- indat$acatches
  # define storage matrices
  grad1val <- matrix(0,nrow=yrce,ncol=nsau,dimnames=list(yearnames,saunames))
  grad4val <- targval <- score1 <- score4 <- scoret <- scoretot <- grad1val
  multTAC <- grad1val
  refpts <- matrix(0,nrow=nsau,ncol=4,
                   dimnames=list(saunames,c("low","trp","high","realtrp")))
  pickceyrs <- match(hsargs$startCE:hsargs$endCE,indat$yearnames)
  actualtarg <- apply(arrce[pickceyrs,],2,quantile,probs=hsargs$targqnt,
                      na.rm=TRUE)
  for (sau in 1:nsau) {  #  sau=1
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
    # Now estimate target from fixed period
    tmp3 <- targscoreconstref(arrce[pickce,sau],actualtarg[sau],
                              mult=hsargs$mult,maxtarg=hsargs$maxtarg[sau])
    nec3 <- length(pickce)
    scrs <- tmp3$scores
    if (nec3 < yrce) scrs <-  c(rep(NA,(yrce-nec3)),scrs)
    scoret[,sau] <- scrs
    targval[,sau] <- arrce[pickce,sau]
    scoretot[, sau] <- pmwts[1] * scoret[, sau] + pmwts[2] * score4[, sau] +
      pmwts[3] * score1[, sau]
    mr3fire <- FALSE
    if (pmwtSwitch > 0) { # meta rule 3  pmwt
      # Metarule to switch weights & update scoretot
      if (all(tail(arrce[pickce, sau], pmwtSwitch) > refpts[sau, 2])) {
        mr3fire <- TRUE
        scoretot[, sau] <- stablewts[1] * scoret[, sau] +
          stablewts[2] * score4[, sau] + stablewts[3] * score1[, sau]
      }
    }
    pick <- floor(scoretot[,sau]) + 1 # add one to get correct hcr[index]
    pick[pick < 1] <- 1
    pick[pick > 10] <- 10
    multTAC[,sau] <- hsargs$hcr[pick] # index implies scores less than index
    refpts[sau,] <- tmp3$rp # the only reference point is the target cpue
    ## Meta Rules section
    if ((metRover > 0) & (!mr3fire)) {     # MR1  over
      if (all(tail(arrce[, sau], metRover) > refpts[sau, 2])) { # CPUE above Target for 2 years
        tmp5 <- diff(arrce[, sau])
        if (all(tail(tmp5, metRover) > 0)) {
          multTAC[yrce, sau] <- hsargs$hcr[pick[yrce]] # already in place
        } else {
          multTAC[yrce, sau] <- 1 # if above but not increasing assign 1 to multTAC for no change
        }
      }
    }  # end of meta rule 1
    if (metRunder > 0) {      # MR2  under
      if (all(tail(arrce[, sau], metRunder) < refpts[sau, 2])) {
        # CPUE above Target for 2 years
        tmp6 <- diff(arrce[, sau])
        if (all(tail(tmp6, metRunder) >= 0))  {
          multTAC[yrce, sau] <-  1 # if below and increasing 2 years assign 1 to multTAC for no change
        } else {
          multTAC[yrce, sau] <- hsargs$hcr[pick[yrce]] # already in place
        }
      }
    } # end meta rule 2
  } # end of sau loop
  acatch <- round(acatches * multTAC[yrce,],0) # to give whole numbers
  TAC <- sum(acatch,na.rm=TRUE)
  details <- list(grad4=grad4val,score4=score4,grad1=grad1val,score1=score1,
                  targval=targval,scoret=scoret,scoretot=scoretot)
  out <- list(acatch=acatch,TAC=TAC,multTAC=multTAC,details=details,
              refpts=refpts)
  return(out)
} # end of constantrefhcr
