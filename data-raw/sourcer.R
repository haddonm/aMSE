

#' Title
#'
#' @param regC
#' @param regD
#' @param glob
#'
#' @return
#' @export
#'
#' @examples
findunfished <- function(regC,regD,glob) {
#  regC=regionC; regD=regionD;  glob=glb
  numpop <- glb$numpop
  inH <- rep(0.0,numpop)
  regD <- runthreeH(regC,regD,glb,inH)
  for (pop in 1:numpop) {
    regC[[pop]]$effB0 <- regD$matureB[1,pop]
    regC[[pop]]$effExB0 <- regD$exploitB[1,pop]
    regC[[pop]]$popq <- regC[[pop]]$popdef["MaxCE"]/regD$exploitB[1,pop]
    regD$cpue[1,pop] <- NA
    regD$deplsB[1,pop] <- 1.0
    regD$depleB[1,pop] <- 1.0
  }
  return(list(regionC=regC,regionD=regD))
} # end of findunfished


#' Title
#'
#' @param regC
#' @param regD
#' @param glob
#' @param catch
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
testequil <- function(regC,regD,glob,inH,verbose=TRUE) {
  Nyrs <- glob$Nyrs
  Nclass <- glob$Nclass
  npop <- glob$numpop
  larvdisp <- glob$larvdisp
  for (yr in 2:Nyrs)
    regD <- oneyearD(regC=regC,regD=regD,Ncl=Nclass,
                     inHt=inH,year=yr,sigmar=1e-08,npop=npop,
                     deltarec=larvdisp)
  if (verbose) {
    if (all(trunc(regD$matureB[1,],3) == trunc(regD$matureB[Nyrs,],3))) {
      print("matureB OK",quote=FALSE)
    } else {
      print("matureB varies",quote=FALSE)
    }
    if (all(trunc(regD$exploitB[1,],3) == trunc(regD$exploitB[Nyrs,],3))) {
      print("exploitB OK",quote=FALSE)
    } else {
      print("exploitB varies",quote=FALSE)
    }
    if (all(trunc(regD$recruit[1,],3) == trunc(regD$recruit[Nyrs,],3))) {
      print("recruitment OK",quote=FALSE)
    } else {
      print("recruitment varies",quote=FALSE)
    }
    if (all(round(regD$deplsB[1,],3) == round(regD$deplsB[Nyrs,],3))) {
      print("spawning depletion OK",quote=FALSE)
    } else {
      print("spawning depletion varies",quote=FALSE)
    }
  }
  return(regD)
} # end of testequil


findF1 <- function(pop,res,location=TRUE) {
  harv <- results[,"AnnH",pop]
  catch <- results[,"Catch",pop]
  grad <- numeric(nH-1)
  for (i in 1:(nH-1)) {
    divisor <- harv[i+1] - harv[i]
    numerator <- catch[i+1] - catch[i]
    grad[i] <- numerator/divisor
  }
  pickF1 <-  which.closest(0.1,grad/grad[1])
  if (location) {
    return(pickF1)
  } else {
    return(grad)
  }
}








