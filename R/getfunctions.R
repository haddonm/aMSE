

##  inzone <- zone1; indexVar <- "Nt"
## gets all LF data from a zone across pops and years
#' @title getregionLF extracts all LF data from a zone across pops and years
#'
#' @description getregionLF extracts all LF data from a region across
#'     all populations for each year. Thus an Nclass x Nyrs X numpop
#'     matrix is compressed into a Nclass x Nyrs matrix
#'
#' @param regD the dynamic part of the region
#' @param glb the global constants
#'
#' @return an Nclass x Nyrs matrix containing LF data across all
#'     populations by year
#' @export
#'
#' @examples
#' \dontrun{
#' print("An example needs to be written")
#' }
getregionLF <- function(regD,glb) { # need to define years
  numpop <- glb$numpop
  Nyrs <- glb$Nyrs
  storeLF <- matrix(0,nrow=glb$Nclass,ncol=Nyrs,
                    dimnames=list(glb$midpts,1:Nyrs))
  for (yr in 1:Nyrs)
    storeLF[,yr] <- rowSums(regD$Nt[,yr,])
  return(storeLF)
} # end of getregionLF

#' @title getlistvar extracts a vector or matrix from regionC
#'
#' @description getlistvar extracts a vector or matrix from regionC.
#'    If a vector of scalars, the names relate to populations, if a
#'    matrix the columns relate to populations. Only Me, R0, B0, effB0,
#'    ExB0, effExB0, MSY, MSYDepl, bLML, popq, SaM, SMUname, popdef,
#'    LML, Maturity, WtL, Emergent, and MatWt are currently valid
#'    choices. The indexvar = popdef would generate a listing of all
#'    the constants. If you only want a single constant from popdefs
#'    then use indexvar2.
#'
#' @param regC the constants components of the simulated region
#' @param indexvar the name of the variable to be extracted; character
#' @param indexvar2 the name of the variable within popdef to extract
#'
#' @return either a vector or matrix of values depending on variable
#' @export
#'
#' @examples
#'  data(region1)
#'  data(constants)
#'  ans <- makeregionC(region=region1,const=constants)
#'  regionC <- ans$regionC
#'  getlistvar(regionC,"MSY")
#'  getlistvar(regionC,"B0")
#'  getlistvar(regionC,"popdef","AvRec")
getlistvar <- function(regC,indexvar,indexvar2="") {
  if (is.character(indexvar)) {
    fields <- c("Me","R0","B0","effB0","ExB0","effExB0","MSY",
                "MSYDepl","bLML","popq","SaM","SMUname",
                "popdef","LML","Maturity","WtL","Emergent","MatWt")
    vects<- c("popdef","LML","Maturity","WtL","Emergent","MatWt")
    pick <- match(indexvar,fields)
    numpop <- length(regC)
    if (is.na(pick)) {
      warning(paste0(fields,"  "))
      stop("Invalid variable name attempted in getlistvar")
    }
    x <- sapply(regC,"[[",indexvar)
    if (fields[pick] %in% vects) {
      colnames(x) <- paste0("p",1:numpop)
    } else {
      names(x) <- paste0("p",1:numpop)
    }
    if ((indexvar == "popdef") & (nchar(indexvar2) > 0)) {
      x <- sapply(regC,"[[",indexvar)[indexvar2,]
      names(x) <- paste0("p",1:numpop)
    }
  } else {
    stop("Non-character variable in indexvar within getlistvar")
  }
  return(x)
} # End of getlistvar



#' @title getunFished - extracts all data relating to year 1; unfished.
#'
#' @description getunFished - extracts all data relating to year 1;
#'    unfished. Of course, this will only work as long as movezoneYear
#'    hasn't been called, which would disrupt the contents of the first
#'    year. But if applied to zone rather than zone1 etc, this should be
#'    fine. The outputs include: R0, A0, B0, B0crypt, popdef, MSY,
#'    MSYDepl, LML, ExB0, Nt, and Ncrypt
#'
#' @param inzone the zone whose unfished state is wanted
#' @param glb contains the global variables used everywhere
#'
#' @return a list of multiple components relating to the unfished stock
#' @export
#' @examples
#' \dontrun{
#'  data(condDat)
#'  glb <- condDat$globals
#'  out <- makeZone(condDat,uplim=0.4)
#'  zone <- out$zone
#'  unfish <- getunFished(zone,glb)
#'  str(unfish,max.level=1)
#' }
getunFished <- function(inzone,glb) {  # inzone=zone
  if((inzone[[1]]$B0 - inzone[[1]]$MatureB[1]) > 0.001) {
    warning("Unfished characters taken from a depleted zone in getunFished \n")
  }
   numpop <- glb$numpop
  nofished <- vector("list",(numpop+1))
  label <- c("R0","B0","B0crypt","popdef","MSY","MSYDepl","LML","ExB0",
             "Nt")
  for (pop in 1:numpop) {
    Nt <- inzone[[pop]]$Nt[,1]

    R0 <- inzone[[pop]]$R0
   # A0 <- inzone[[pop]]$A0
    B0 <- inzone[[pop]]$B0
    B0crypt <- B0 -
            sum(inzone[[pop]]$Emergent * Nt * inzone[[pop]]$MatWt)/1e06
    popdef <- inzone[[pop]]$popdef
    MSY <- inzone[[pop]]$MSY
    MSYDepl <- inzone[[pop]]$MSYDepl
    LML <- inzone[[pop]]$LML[1]
    ExB0 <- inzone[[pop]]$ExB0
    ans <- list(R0,B0,B0crypt,popdef,MSY,MSYDepl,LML,ExB0,Nt)
    names(ans) <- label
    nofished[[pop]] <- ans
  }
  R0 <- 0
 # A0 <- 0
  B0 <- 0
  B0crypt <- 0
  MSY <- 0
  ExB0 <- 0
  Nt <- numeric(glb$Nclass)
  #Ncrypt <- numeric(glb$Nclass)
  for (pop in 1:numpop) { #sums across all pops to go into place numpop+1
    R0 <- R0 + nofished[[pop]]$R0
  #  A0 <- A0 + nofished[[pop]]$A0
    B0 <- B0 + nofished[[pop]]$B0  ## or use: sum(sapply(zone,"[[","B0"))
    B0crypt <- B0crypt + nofished[[pop]]$B0crypt
    MSY <- MSY + nofished[[pop]]$MSY
    ExB0 <- ExB0 + nofished[[pop]]$ExB0
    Nt <- Nt + nofished[[pop]]$Nt
   # Ncrypt <- Ncrypt + nofished[[pop]]$Ncrypt
  }
  zoneT <- list(R0,B0,B0crypt,MSY,ExB0,Nt)
  names(zoneT) <- c("R0","B0","B0crypt","MSY","ExB0","Nt")
  nofished[[numpop+1]] <- zoneT
  return(nofished)
}  # End of getUnfished

#' @title getzoneProps extracts the depletion level for a given year
#'
#' @description getzoneProps extracts the depletion level along with an
#'     array of other statistics for a given year.
#'
#' @param regC the constants for the region being simulated
#' @param glb the global constants
#' @param year which years information is wanted
#'
#' @return a list of size objects
#' @export
#'
#' @examples
#' \dontrun{
#'  print("wait on internal data sets")
#' }
getzoneProps <- function(regC,glb,year=1) { # inzone=zone; glb=glb; year=1
  numpop <- glb$numpop
  Nclass <- glb$Nclass
  totB0 <- sum(sapply(regC,"[[","B0"))
  ExB0 <- sum(sapply(regC,"[[","ExB0"))
  matB <- 0.0
  ExB <- 0.0
  legalMatB <- 0.0
  ce <- 0.0
  for (pop in 1:numpop) {
    ExB <- ExB + regC[[pop]]$ExploitB[year]
    matB <- matB + regC[[pop]]$MatureB[year]
    ce <- ce + log(regC[[pop]]$cpue[year])
    MatWt <- regC[[pop]]$MatWt
    lml <- regC[[pop]]$LML[year]
    pick <- which(glb$midpts >= lml)
    first <- pick[1]
    legalMatB <- legalMatB + (sum(MatWt[first:Nclass]*
                         regC[[pop]]$Nt[first:Nclass,year])/1000000.0)
  }
  deplet <- matB/totB0
  depletEx <- ExB/ExB0
  expCE <- exp(ce/numpop)
  ans <- c(deplet,depletEx,ExB,expCE,(legalMatB/totB0))
  names(ans) <- c("SpBDepl","ExBDepl","ExploitB","ExpCE","LegalMatB")
  return(ans)
}  # end of getzoneProps
