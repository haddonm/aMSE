

#' @title doprojection applies the Tasmanian MCDA to a simulated zone
#'
#' @description doprojection applies the Tasmanian MCDA to a simulated zone. It
#'     is still in need of optimization as 1000 iterations would currently take
#'     about 20 minutes.
#'
#' @param zoneCP the modified zoneC object ready for projections
#' @param zoneDP the modified zoneD object ready for projections
#' @param glob the global variables object
#' @param ctrl the ctrl object
#' @param projyrs the number of years of projection
#' @param applyHS the function describing the harvest strategy to be used
#' @param hsargs a list of any arguments used by applyHS
#' @param projpms the matrices of performance measures from calibrateHS
#' @param inityrs the number of years kept from the period of conditioning,
#' @param ... the ellipsis should contain any additional arguments
#'
#' @return the filled in zoneDP list
#' @export
#'
#' @examples
#' print("wait on data files")
#' # zoneCP=zoneCP;zoneDP=zoneDR;glob=glb;ctrl=ctrl;projyrs=projC$projyrs;inityrs=projC$inityrs
#' # wid=4;targqnt=0.55;pmwts=c(0.65, 0.25,0.1);hcr = c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2)
doprojection <- function(zoneCP,zoneDP,glob,ctrl,projyrs,applyHS,hsargs,
                         projpms,inityrs,...) {
  # get  important constants
  sigmaR <- ctrl$withsigR # needed to add recruitment variation
  npop <- glob$numpop
  nsau <- glob$nSAU
  Ncl <- glob$Nclass
  nyrs <- projyrs
  movem <- glob$move
  reps <- ctrl$reps
  sauindex <- glob$sauindex
  matb <- numeric(npop)                   # use the same initial TAC for all reps
  origTAC <- mean(colSums(zoneDP$catch[1,,])) # mean sum of catches in last year

  saucatch <- array(0,dim=c(nyrs,nsau,reps))
  saucpue <- saucatch

  # for (iter in 1:reps) { # generate SAU total catches and catch-weighted cpue iter=1
  #   for (yr in 1:inityrs) { # iter=1; yr=4
  #     saucatch[yr,,iter] <- tapply(zoneDP$catch[yr,,iter],sauindex,sum,na.rm=TRUE)
  #     wts <- zoneDP$catch[yr,,iter]/(saucatch[yr,sauindex,iter])
  #     saucpue[yr,,iter] <- tapply((zoneDP$cpue[yr,,iter] * wts),sauindex,sum,na.rm=TRUE)
  #     if (yr >= hsargs$wid) {
  #       grad1[yr,,iter] <- apply(saucpue[1:yr,,iter],2,getgradone,yr=yr)
  #       grad4[yr,,iter] <- apply(saucpue[1:yr,,iter],2,getgradwid,yr=yr,wid=hsargs$wid)
  #     }
  #   }
  #   for (sau in 1:nsau)
  #     targetce[sau] <- targscore(saucpue[1:inityrs,sau,iter],qnt=hsargs$targqnt)
  #   targsc[1:inityrs,,iter] <- targetce
  # }
  # now do replicates, updating saucatch and saucpue each year
  if(ctrl$randseedP > 0) set.seed(ctrl$randseedP) # set random seed if desired
  for (iter in 1:reps) {
    TAC <- origTAC  # use the same original TAC for each replicate
    for (year in (inityrs+1):nyrs) { # iter=1; year=11
      #  catpop <- colSums(zoneDP$catch[1:(year - 1),,iter])
      inexpB <- zoneDP$exploitB[(year - 1),,iter]
      sauexpB <- tapply(inexpB,sauindex,sum,na.rm=TRUE)
      catbysau <- TAC * sauexpB/sum(sauexpB)  # no error initially
      multh <- apply(saucpue[1:(year-1),,1],2,applyHS,yr=(year-1)) # apply mcdahcr
      TAC <- sum(catbysau * multh)
      divererr <- sauexpB * exp(rnorm(nsau,mean=0,sd=ctrl$withsigB))
      catbysau <- TAC * (divererr/sum(divererr)) # currently no error on TAC
      catbypop <- catbysau[sauindex] * (inexpB/sauexpB[sauindex]) # no error on pops
      for (popn in 1:npop) { # year=11; iter=1; pop=1
        out <- oneyearcat(inpopC=zoneCP[[popn]],inNt=zoneDP$Nt[,year-1,popn,iter],
                          Nclass=Ncl,incat=catbypop[popn],yr=year)
        zoneDP$exploitB[year,popn,iter] <- out$ExploitB
        zoneDP$matureB[year,popn,iter] <- out$MatureB
        zoneDP$catch[year,popn,iter] <- out$Catch
        zoneDP$harvestR[year,popn,iter] <- out$Harvest
        zoneDP$cpue[year,popn,iter] <- out$ce
        zoneDP$Nt[,year,popn,iter] <- out$Nt
        zoneDP$catchN[,year,popn,iter] <- out$CatchN
        matb[popn] <- out$MatureB
      } # pop
      steep <- getvect(zoneCP,"steeph")
      r0 <- sapply(zoneCP,"[[","R0")
      b0 <- sapply(zoneCP,"[[","B0")
      recs <- oneyearrec(steep,r0,b0,matb,sigR=ctrl$withsigR)
      newrecs <- movem %*% recs
      zoneDP$recruit[year,,iter] <- newrecs
      zoneDP$Nt[1,year,,iter] <- newrecs
      zoneDP$deplsB[year,,iter] <- zoneDP$matureB[year,,iter]/b0
      zoneDP$depleB[year,,iter] <- zoneDP$exploitB[year,,iter]/sapply(zoneCP,"[[","ExB0")
      saucatch[year,,iter] <- tapply(zoneDP$catch[year,,iter],sauindex,sum,na.rm=TRUE)
      wts <- zoneDP$catch[year,,iter]/(saucatch[year,sauindex,iter])
      saucpue[year,,iter] <- tapply((zoneDP$cpue[year,,iter] * wts),sauindex,sum,na.rm=TRUE)
    }   # year loop        zoneDR$matureB[,,1]
  }     # rep loop
  zoneDP$cesau <- saucpue
  return(zoneDP=zoneDP)
} # end of doprojection

#' @title product is the productivity curve matrix from doproduction
#'
#' @description product is the productivity curve matrix from
#'     doproduction when the example zone is generated using the
#'     inbuilt datasets ctrl, zone1, and constants. The slowest
#'     part of building the whole is to use the modregC function
#'     to adjust the zoneC and generate the production array. To
#'     save that time in the examples (to avoid time limits on
#'     examples should this package go to CRAN), then this dataset can
#'     be used instead. This is a three dimensional array of
#'     productivity variables.
#'
#' @name product
#'
#' @docType data
#'
#' @section contents:
#' \itemize{
#'   \item harvestrate the initial harvest rates applied
#'   \item productivity variables ExB, MatB, AnnH, Catch, Deplet, RelCE
#'   \item population the index of each population
#' }
#'
#' @examples
#'  data(product)
#'  product[1:20,,1]
NULL

#' @title testzoneC is a zone list made up of 6 equilibrium populations
#'
#' @description testzoneC is a zone list made up of 6 equilibrium
#'     populations. These have been run with a laral dispersal rate of
#'     0.03 so the change from B0 to effB0 is not great, but still
#'     required for an initial equilibrium. This is here to simplify
#'     the internal testing of funcitons that require a completed
#'     zone starting at equilibrium. Its name is to avoid conflict
#'     with any actual use of zoneC. use str(testzoneC, max.level=1)
#'     to see its format. It can be expected to be used with testzoneD
#'
#' @name testzoneC
#'
#' @docType data
#'
#' @section Subjects:
#'  \itemize{
#'    \item testing of functions that require a full zone
#'    \item initial equilibrium
#'  }
#'  @export
#'
#' @examples
#'  data(testzoneC)
#'  data(testzoneD)
#'  data(zone1)
#'  glb <- zone1$globals
#'  r0 <- getvar(testzoneC,"R0")
#'  move <- makemove(glb$numpop,r0,glb$larvdisp)
#'  glb$move <- move
#'  ans <- testequil(testzoneC, testzoneD, glb)
#'  str(testzoneC[[1]])
NULL

#' @title testzoneD is a list of 8 matrices and 2 arrays defining the dynamics of a zone
#'
#' @description testzoneD is a list of 8 matrices and 2 arrays defining
#'     the dynamics of a zone. These have been run with a larval
#'     dispersal rate of 0.03 to achieve an initial equilibrium. This
#'     is here to simplify the internal testing of functions that
#'     require a completed zone starting at equilibrium. Its name is
#'     to avoid conflict with any actual use of zoneD. use
#'     str(testzoneD, max.level=1) to see its format. It can be
#'     expected to be used with testzoneC.
#'
#' @name testzoneD
#'
#' @docType data
#'
#' @section Subjects:
#'  \itemize{
#'    \item testing of functions that require a full zone
#'    \item initial equilibrium
#'  }
#'  @export
#'
#' @examples
#'  data(testzoneC)
#'  data(testzoneD)
#'  data(zone1)
#'  glb <- zone1$globals
#'  r0 <- getvar(testzoneC,"R0")
#'  move <- makemove(glb$numpop,r0,glb$larvdisp)
#'  glb$move <- move
#'  ans <- testequil(testzoneC, testzoneD, glb)
#'  str(testzoneD)
NULL





doprojection <- function(zoneCP,zoneDP,glob,ctrl,projyrs,applyHS,hsargs,
                         projpms) {
  zoneCP=zoneCP
  zoneDP=zoneDP
  glob=glb
  ctrl=ctrl
  projyrs=projC$projyrs
  applyHS=mcdahcr
  hsargs=hsargs
  projpms=cmcda
  # get important constants
  sigmaR <- ctrl$withsigR # needed to add recruitment variation
  sigmaB <- ctrl$withsigB
  npop <- glob$numpop
  nsau <- glob$nSAU
  Ncl <- glob$Nclass
  nyrs <- projyrs
  movem <- glob$move
  reps <- ctrl$reps
  sauindex <- glob$sauindex
  matb <- numeric(npop)
  origTAC <- colSums(zoneDP$catch[1,,]) # mean sum of catches in last year
  # now do replicates, updating saucatch and saucpue each year
  if (ctrl$randseedP > 0) set.seed(ctrl$randseedP) # set random seed if desired
  for (iter in 1:reps) {
    TAC <- origTAC[iter]
    for (year in 2:nyrs) { # iter=1; year=11
      catbysau <- catchbysau(inexpB=zoneDP$exploitB[(year - 1),,iter],sauindex,TAC) # no error initially
      multh <- apply(saucpue[1:(year-1),,1],2,applyHS,yr=(year-1)) # apply mcdahcr
      TAC <- sum(catbysau * multh)
      divererr <- sauexpB * exp(rnorm(nsau,mean=0,sd=ctrl$withsigB))
      catbysau <- TAC * (divererr/sum(divererr)) # currently no error on TAC
      catbypop <- catbysau[sauindex] * (inexpB/sauexpB[sauindex]) # no error on pops
      for (popn in 1:npop) { # year=11; iter=1; pop=1
        out <- oneyearcat(inpopC=zoneCP[[popn]],inNt=zoneDP$Nt[,year-1,popn,iter],
                          Nclass=Ncl,incat=catbypop[popn],yr=year)
        zoneDP$exploitB[year,popn,iter] <- out$ExploitB
        zoneDP$matureB[year,popn,iter] <- out$MatureB
        zoneDP$catch[year,popn,iter] <- out$Catch
        zoneDP$harvestR[year,popn,iter] <- out$Harvest
        zoneDP$cpue[year,popn,iter] <- out$ce
        zoneDP$Nt[,year,popn,iter] <- out$Nt
        zoneDP$catchN[,year,popn,iter] <- out$CatchN
        matb[popn] <- out$MatureB
      } # pop
      steep <- getvect(zoneCP,"steeph")
      r0 <- sapply(zoneCP,"[[","R0")
      b0 <- sapply(zoneCP,"[[","B0")
      recs <- oneyearrec(steep,r0,b0,matb,sigR=sigmar)
      newrecs <- movem %*% recs
      zoneDP$recruit[year,,iter] <- newrecs
      zoneDP$Nt[1,year,,iter] <- newrecs
      zoneDP$deplsB[year,,iter] <- zoneDP$matureB[year,,iter]/b0
      zoneDP$depleB[year,,iter] <- zoneDP$exploitB[year,,iter]/sapply(zoneCP,"[[","ExB0")
      saucatch[year,,iter] <- tapply(zoneDP$catch[year,,iter],sauindex,sum,na.rm=TRUE)
      wts <- zoneDP$catch[year,,iter]/(saucatch[year,sauindex,iter])
      saucpue[year,,iter] <- tapply((zoneDP$cpue[year,,iter] * wts),sauindex,sum,na.rm=TRUE)
    }   # year loop        zoneDR$matureB[,,1]
  }     # rep loop
  zoneDP$cesau <- saucpue
  return(zoneDP=zoneDP)
} # end of doprojection




doproj <- function(zoneC,zoneDD,glb,ctrl,projC,applyHS=mcdahcr,HSargs,
                   histpms) {
  ans <- calcprojsel(zoneC,projC,glb) # calculate selectivity for projections
  sel <- ans$projSel
  selwt <- ans$projSelWt
  reps <- ctrl$reps
  pyrs <- projC$projyrs
  zoneDP <- makezoneDP(pyrs,reps,glb,zoneDD) # make object to hold projections
  sigmar <- ctrl$withsigR # needed to add recruitment variation
  npop <- glb$numpop
  nsau <- glb$nSAU
  Ncl <- glb$Nclass

  lastyr <- glb$Nyrs
  histyrs <- projC$inityrs:lastyr
  movem <- glb$move
  reps <- ctrl$reps
  matb <- numeric(npop)
  origTAC <- sum(zoneDD$catch[lastyr,]) # Assumes the TAC is always caught!
  sauindex <- glb$sauindex
  saucatch <- array(0,dim=c(pyrs,nsau,reps))
  saucpue <- saucatch
  grad4 <- array(0,dim=c(pyrs,nsau,reps)) # declare arrays to store PMs
  grad1 <- grad4
  targsc <- grad4
  targetce <- numeric(nsau)
  for (iter in 1:reps) {
    for (year in 1:pyrs) {  # iter=1; year = 1
      ans <- getmcdadatap(zoneDD$catch,zoneDD$cpue,
                          zoneDP$catch[,,iter],zoneDP$cpue[,,iter],
                          histyrs,year)


    }
  }
  #calibrate the HCR on zoneDD before the projections
  for (yr in 1:inityrs) { # iter=1; yr=4
    saucatch[yr,,iter] <- tapply(zoneDP$catch[yr,,iter],sauindex,sum,na.rm=TRUE)
    wts <- zoneDP$catch[yr,,iter]/(saucatch[yr,sauindex,iter])
    saucpue[yr,,iter] <- tapply((zoneDP$cpue[yr,,iter] * wts),sauindex,sum,na.rm=TRUE)
    if (yr >= wid) {
      grad1[yr,,iter] <- apply(saucpue[1:yr,,iter],2,getgradone,yr=yr)
      grad4[yr,,iter] <- apply(saucpue[1:yr,,iter],2,getgradwid,yr=yr)
    }
  }
  for (sau in 1:nsau)
    targetce[sau] <- targscore(saucpue[1:inityrs,sau,iter],qnt=targqnt)$result
  targsc[1:inityrs,,iter] <- targetce

}



