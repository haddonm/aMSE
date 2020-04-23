



cat("recthreshold, 0.0000001  \n",file=filename, append=TRUE)
cat("hcrLabel,  mcdaHCR  \n",file=filename, append=TRUE)
cat("mcdaHCR,  TRUE   \n",file=filename, append=TRUE)
cat("ConstC,   FALSE  \n",file=filename, append=TRUE)
cat("ConstH,   FALSE  \n",file=filename, append=TRUE)
TACadj <- c(0.75,0.8,0.85,0.9,0.95,1.0,1.05,1.10,1.15,1.2,1.25)
cat(paste0("TACadj,  ",makeLabel(TACadj,insep=",")," \n"),file=filename,
    append=TRUE)
mcdaWeights <- c(0.4,0.5,0.1)
cat(paste0("mcdaWts,  ",makeLabel(mcdaWeights,insep=",")," \n"),
    file=filename, append=TRUE)
postmcdaWeights <- mcdaWeights  # wts after target
cat(paste0("postmcdaWts, ",makeLabel(postmcdaWeights,insep=",")," \n"),
    file=filename,
    append=TRUE)
cat(paste0("cpuePeriod,  4   \n"),file=filename, append=TRUE)
cat(paste0("maxGrad4,   0.2   \n"),file=filename, append=TRUE)
cat(paste0("maxRate1,   0.4   \n"),file=filename, append=TRUE)
# Target Based HCR
targCE <- rep(80.0,4)
cat(paste0("CETarg,  ",makeLabel(targCE,insep=",")," \n"),file=filename,
    append=TRUE)   # if using rep then same target in each block
deltCE <- rep(45,4)
cat(paste0("deltaCE,  ",makeLabel(deltCE,insep=",")," \n"),file=filename,
    append=TRUE)
cat(paste0("implementE, 0.0  \n"),file=filename, append=TRUE) # 0 = no delay, 1 = 1 year delay, used in all HCR
cat(paste0("LRPTAC,  TRUE   \n"),file=filename, append=TRUE)
cat(paste0("TACLower, 500.0   \n"),file=filename, append=TRUE)# TAC can only decline to >= TACLower * origTAC
cat(paste0("TACUpper, 1000.0  \n"),file=filename, append=TRUE)
cat(paste0("refyr,  20     \n"),file=filename, append=TRUE)
