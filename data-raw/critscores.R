#Loop through years
for (i in period:numyears) { # edited 20/09/2015
  # for (i in period+1:numyears) {
  if (!is.na(subdat$stndcpue[i])) {

    yr1 <- i-(period-1)
    activeyear <- i
    yrSeq <- seq(yr1:activeyear)

    ## * determine score for CPUE Target PM ####
    subdat$CEtarg[i] <- cpuex[2]
    subdat$CEval[i] <- subdat$stndcpue[activeyear]
    if(!is.na(subdat$stndcpue[activeyear])) {
      subdat$CritCPUE[i] <-
        approx(cpuex, cpuey, xout = subdat$stndcpue[activeyear], yleft=0,yright=10)$y
    }

    ## * determine score for CPUE Gradient4 PM ####
    percCE <- subdat$stndcpue[yr1:activeyear]/subdat$stndcpue[yr1]

    possibleError <- tryCatch(
      slopeModel <- lm(percCE ~ yrSeq),
      error=function(e) e
    ) # End tryCatch catch

    if(inherits(possibleError, "error")) {next} # skip to next value in for loop if error found

    grad <- slopeModel$coeff[2]
    subdat$grad4[i] <- grad
    subdat$CritCPUEGrad[i] <- approx(cpuegradx, cpuegrady, xout = grad , yleft=0,yright=10)$y

    ## * determine score for CPUE Gradient1 PM ####
    pc <- (subdat$stndcpue[activeyear]- subdat$stndcpue[activeyear-1])/subdat$stndcpue[activeyear-1]
    subdat$grad1[i] <- pc
    subdat$CritCPUEChng[i] <- approx(cpuechngx, cpuechngy, xout = pc, yleft=0,yright=10)$y

    ## Calculate combined score from individual PM score * weighting ####
    subdat$ComScore[i] <- subdat$CritCPUE[i]*subdat$wghtCEtarg[i] +
      subdat$CritCPUEGrad[i]*subdat$wghtgrad4[i] +
      subdat$CritCPUEChng[i]*subdat$wghtgrad1[i]


    ## Pass score through to the Control Rule
    possibleError <- tryCatch(
      {
        pick <- which(subdat$ComScore[i] >= HCR$LowerCut & subdat$ComScore[i] < HCR$UpperCut)
        subdat$TACAdjust[i] <- HCR$Action[pick]

      },
      error=function(e) e
    ) # End tryCatch HCR

    if(inherits(possibleError, "error")) {next} # skip to next value in for loop if error found


    ## caluculate TACC for following Fishing Year ####
    if (i >= period & i <= numyears ) {
      subdat$TAC[i+1] <- subdat$TACAdjust[i] * subdat$Catch[i]
    }

  } #End IF is.na loop

}
