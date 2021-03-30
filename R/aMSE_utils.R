

# rutilsMH::listFunctions("C:/Users/User/Dropbox/A_code/aMSE/aMSE_utils.R")


#' @title poptosau converts projected population dynamics to SAU scale results
#'
#' @description poptosau the MSE dynamics are run at the population level but
#'     all management decisions are made at the SAU level, hence the results of
#'     the projections need to be translated into results at the SAU level.
#'     poptosau uses the sauindex to sum the variables that can be summed, which
#'     include matureB, exploitB, catch, and recruit, and combines their
#'     population results to form their respective SAU results. Thus, for an
#'     input array of dimensions [projyrs,numpop,reps] one receives back an
#'     array of [projyrs, nSAU,reps]
#'
#' @param invar either zoneDP$ matureB, exploitB, catch, or recruit
#' @param glb the global constants object
#'
#' @return a results array of dimension [projyrs, nSAU,reps]
#' @export
#'
#' @examples
#' print("wait on appropriate internal datasets")
poptosau <- function(invar,glb) {  # invar=zoneDP$matureB; glb=glb
  numpop <- glb$numpop
  nsau <- glb$nSAU
  nyrs <- dim(invar)[1]
  reps <- dim(invar)[3]
  sauindex <- glb$sauindex
  saunames <- glb$saunames
  result <- array(0,dim=c(nyrs,nsau,reps),
                  dimnames=list(1:nyrs,saunames,1:reps)) #aspirational catches
  for (iter in 1:reps)
    for (yr in 1:nyrs)
      result[yr,,iter] <- tapply(invar[yr,,iter],sauindex,sum,na.rm=TRUE)
  return(result)
} # end of poptosa

#' @title prepareDDNt converts the population based historical Nt to SAU-based
#'
#' @description prepareDDNt processes the zoneDD, following the conditioning on
#'     the historical fishery data, so that the population numbers-at-size are
#'     converted to SAU-based numbers-at-size. It does this for both Nt and
#'     catchN. This is to allow analysis, tabulation, and plotting at the SAU
#'     scale. Input data is 3D  Nclass x years x numpop
#'
#' @param inNt the 3 dimensional array of population numbers-at-size
#' @param incatchN the 3 dimensional array of catch numbers-at-size
#' @param glb the global constants object
#'
#' @return a list of Nt and catchN x SAU
#' @export
#'
#' @examples
#' print("wait on suitable internal data-sets")
prepareDDNt <- function(inNt,incatchN,glb) { # Nt=zoneDD$Nt; catchN=zoneDD$catchN; glb=glb
  sauindex <- glb$sauindex
  nsau <- glb$nSAU
  nyrs <- glb$Nyrs
  Nc <- glb$Nclass
  catchN <- array(data=0,dim=c(Nc,nyrs,nsau), # define some arrays
                  dimnames=list(glb$midpts,1:nyrs,glb$saunames))
  Nt <- array(data=0,dim=c(Nc,nyrs,nsau),
              dimnames=list(glb$midpts,1:nyrs,glb$saunames))
  for (sau in 1:nsau) { #  yr=1; sau = 1
    pick <- which(sauindex %in% sau)
    for (yr in 1:nyrs) {
      catchN[,yr,sau] <- rowSums(incatchN[,yr,pick])
      Nt[,yr,sau] <- rowSums(inNt[,yr,pick])
    }
  }
  return(list(Nt=Nt,catchN=catchN))
} # end of prepareDDNt

#' @title summarizeprod generates a summary of the productivity properties
#'
#' @description summarizeprod generates a summary of the productivity properties
#'     by examining the productivity matrix for each SAU/population and
#'     extracting the Bmsy, annualH, MSY, Depletion, and RelCE at the maximum
#'     catch level (which approximates the MSY). It summarizes the total zone
#'     by summing all the productivity matrices and search for the largest
#'     catch again. It generates estimates of the annualH, depletion and RelCE
#'     by using a weighted average of those values from teh separate SAU or
#'     populations, where the weighting is the proportion of the sum of the
#'     MSYs taken in each sau or population. This latter is only an
#'     approximation but provides at least an indication.
#'
#' @param product The productivity array from doproduction containing the
#'     range of imposed harvest rates, and the resulting outputs for each
#'     population
#' @param saunames the names of the different SAU
#'
#' @return a matrix containing the approximate productivity matrix for the zone
#' @export
#'
#' @examples
#' \dontrun{
#' data(zone)
#' product <- zone$product
#' zoneprod <- summarizeprod(product,saunames=zone$zone1$SAUnames)
#' round(zoneprod,3)
#' }
summarizeprod <- function(product,saunames) { # product=zone$product; saunames=zone$zone1$SAUnames
  numrow <- dim(product)[1]
  numpop <- dim(product)[3]
  prodrows <- rownames(product[,,1])
  prodcols <- colnames(product[,,1])
 # nSAU <- length()  # here numpop = nSAU
  columns <- c("Bmsy","AnnH","MSY","Deplet","RelCE")
  rows <- c(saunames,"zone")
  nrows <- length(rows)
  ans <- matrix(0,nrow=nrows,ncol=length(columns),dimnames=list(rows,columns))
  for (sau in 1:numpop) {
    sauprod <- product[,,sau]
    pick <- which(sauprod[,"Catch"] == max(sauprod[,"Catch"],na.rm=TRUE))
    ans[sau,] <- sauprod[pick,2:6]
  }
  zoneprod <- matrix(0,nrow=numrow,ncol=length(prodcols),
                     dimnames=list(prodrows,prodcols))
  for (i in 1:numpop) zoneprod <- zoneprod + product[,,i]
  pick <-  which(zoneprod[,"Catch"] == max(zoneprod[,"Catch"],na.rm=TRUE))
  ans[nrows,] <- zoneprod[pick,2:6]
  pmsy <- ans[1:numpop,"MSY"]/sum(ans[1:numpop,"MSY"],na.rm=TRUE)
  ans[nrows,"AnnH"] <- sum(ans[1:numpop,"AnnH"] * pmsy,na.rm=TRUE)
  ans[nrows,"Deplet"] <- sum(ans[1:numpop,"Deplet"] * pmsy,na.rm=TRUE)
  ans[nrows,"RelCE"] <- sum(ans[1:numpop,"RelCE"] * pmsy,na.rm=TRUE)
  return(ans)
} # end of summarizeprod

#' @title sumpop2sau gathers population data into sau data using sauindex
#'
#' @param invect a vector of population values for a given variable
#' @param sauindex the indices of each sau for each population
#'
#' @return a vector of length nsau containing the sum of population values
#'     for each sau
#' @export
#'
#' @examples
#' vect <- c(5.8,6.2,13.2,23.8,3.3,3.7,29.7,26.3,38.9,9.1)
#' sauind <- c(1,1,2,2,3,3,4,4,5,5)
#' sumpop2sau(vect,sauind)   # should be 12 37 7 56 48
sumpop2sau <- function(invect,sauindex) {
  return(tapply(invect,sauindex,sum,na.rm=TRUE))
} # end of sumpop2sau

#' @title zonetosau translates the zonexpop objects to zonexsau objects
#'
#' @description zonetosau combines the dynamic results for each variable so
#'     that results by population become results by SAU. matureB, exploitB,
#'     catch, recruit, catchN, and Nt are simple summations of the totals for
#'     each population into their respective SAU. The harvest rate would be
#'     end of year or beginning of year estimates derived from dividing the
#'     catchxsau by the exploitable biomass x sau. Similarly the deplsB and
#'     depleB are the end of year matureB and exploitB divided by their
#'     respective unfished estimated by SAU obtained using getvar(zoneC,"B0").
#'
#' @param inzone one of the zone dynamics objects made up of populations
#' @param glb the object containng the global constants
#' @param B0 the estimate B0 for each population use getvar(zoneC,"B0")
#' @param ExB0 the estimate ExB0 for each population use getvar(zoneC,"ExB0")
#'
#' @return a list of dynamics variables by SAU
#' @export
#'
#' @examples
#' print("wait on suitable internal data sets ")
zonetosau <- function(inzone,glb, B0, ExB0) { # inzone=zoneDP; glb=glb; B0=B0; ExB0=ExB0
  nsau <- glb$nSAU
  sauindex <- glb$sauindex
  saunames <- glb$saunames
  N <- glb$Nclass
  invar <- inzone$matureB
  nyrs <- dim(invar)[1]
  reps <- dim(invar)[3]
  sauB0 <- tapply(B0,sauindex,sum)
  sauExB0 <- tapply(ExB0,sauindex,sum)
  matureB <- poptosau(invar,glb)
  exploitB <- poptosau(inzone$exploitB,glb)
  recruit <- poptosau(inzone$recruit,glb)
  catchN <- array(data=0,dim=c(N,nyrs,nsau,reps), # define some arrays
                  dimnames=list(glb$midpts,1:nyrs,saunames,1:reps))
  Nt <- array(data=0,dim=c(N,nyrs,nsau,reps),
              dimnames=list(glb$midpts,1:nyrs,saunames,1:reps))
  catch <- inzone$acatch
  cpue <- catch
  harvestR <- catch
  deplsB <- catch
  depleB <- catch
  popcat <- inzone$catch
  popce <- inzone$cpue
  for (iter in 1:reps) {
    for (yr in 1:nyrs) {
      saudyn <- poptosauCE(popcat[yr,,iter],popce[yr,,iter],sauindex)
      cpue[yr,,iter] <- saudyn$saucpue
      catch[yr,,iter] <- saudyn$saucatch
      deplsB[yr,,iter] <- matureB[yr,,iter]/sauB0
      depleB[yr,,iter] <- exploitB[yr,,iter]/sauExB0
      for (sau in 1:nsau) { #  iter=1; yr=1; sau = 1
        pick <- which(sauindex %in% sau)
        catchN[,yr,sau,iter] <- rowSums(inzone$catchN[,yr,pick,iter])
        Nt[,yr,sau,iter] <- rowSums(inzone$Nt[,yr,pick,iter])
      }
    }
    harvestR[,,iter] <- catch[,,iter]/exploitB[,,iter]
  }
  outsau <- list(matureB=matureB,exploitB=exploitB,catch=catch,
                 acatch=inzone$acatch,harvestR=harvestR,cpue=cpue,
                 recruit=recruit,deplsB=deplsB,depleB=depleB,
                 catchN=catchN,Nt=Nt)
  return(outsau)
} # end of zonetosau



