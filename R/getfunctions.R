
#' @title getextension gets the 3 letter file extension from a filename
#'
#' @description getextension is a utility function that is used with
#'     logfilename to determine whether one is dealing with a csv or a
#'     png file. depending on which it will return either a 'table' or
#'     'plot' value. Any other extension type will throw an error and
#'     stop execution.
#'
#' @param filename the filename, with path attached, whose extension
#'     needs to be determined.
#'
#' @return a character string of either 'table' or 'plot'
#' @export
#'
#' @examples
#' filen <- "not_a_real_file.png"
#' getextension(filen)
#' filen <- "another_unreal_file.csv"
#' getextension(filen)
getextension <- function(filename) {
  lenc <- nchar(filename)
  exten <- substr(filename,(lenc-2),lenc)
  type <- ""
  if (exten == "png") type <- "plot"
  if (exten == "csv") type <- "table"
  if (nchar(type) == 0)
    stop("A file added to results must be either 'png' or 'csv'")
  return(type)
} # end of getextension

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

#' @title getregionprops extracts the depletion level for a given year
#'
#' @description getregionprops extracts a set of properties for each
#'     population and summarizes them for the region as well. These
#'     properties include effB0, matureB, legalmatB, propprot, MSY,
#'     effexB0, exploitB, SpBDepl, ExBDepl, legalDepl, MSYDepl, LML,
#'     and bLML. See the model documentation for the meaning of each
#'     of these.
#'
#' @param regC the constants for the region being simulated
#' @param regD the dynamic part of the region being simulated
#' @param glb the global constants
#' @param year which years information is wanted, default=1
#'
#' @return a matrix of population properties with a column of totals
#' @export
#'
#' @examples
#' data(testregC)
#' data(testregD)
#' data(region1)
#' glb <- region1$globals
#' round(getregionprops(testregC,testregD,glb),4)
getregionprops <- function(regC,regD,glb,year=1) { # regC=regionC; regD=regionD; glb=glb; year=1
  numpop <- glb$numpop
  Nclass <- glb$Nclass
  effB0 <- sapply(regC,"[[","effB0")
  effExB0 <- sapply(regC,"[[","effExB0")
  blml <- sapply(regC,"[[","bLML")
  msy <- sapply(regC,"[[","MSY")
  msydepl <- sapply(regC,"[[","MSYDepl")
  matB <- regD$matureB[year,]
  ExB <- regD$exploitB[year,]
  harvestR <- regD$harvestR[year,]
  catch <- regD$catch[year,]
  legalmatB <- numeric(numpop); names(legalmatB) <- 1:numpop
  lml <- sapply(regC,"[[","LML")[year,]
  for (pop in 1:numpop) { #   pop=1
    MatWt <- regC[[pop]]$MatWt
    first <- which(glb$midpts >= lml[pop])[1]
    legalmatB[pop] <- sum(regC[[pop]]$MatWt[first:Nclass] *
                            regD$Nt[first:Nclass,year,pop])/1000000.0
  }
  deplet <- matB/effB0
  depletEx <- ExB/effExB0
  legaldepl <- legalmatB/effB0
  propprot <- (matB - legalmatB)/matB
  label <- c("effB0","matureB","legalmatB","propprot","MSY",
             "effexB0","exploitB","SpBDepl","ExBDepl",
             "legalDepl","MSYDepl","LML","bLML","harvestR","catch")
  ans <- rbind(effB0,matB,legalmatB,propprot,msy,effExB0,ExB,deplet,
               depletEx,legaldepl,msydepl,lml,blml,harvestR,catch)
  rownames(ans) <- label
  tot <- numeric(length(rownames(ans))); names(tot) <- label
  tot[c(1:3,5:7)] <- rowSums(ans)[c(1:3,5:7)]
  tot["propprot"] <- (tot["matureB"] - tot["legalmatB"])/tot["matureB"]
  wgts <- ans["effB0",]/tot["effB0"]
  tot["SpBDepl"] <- sum(wgts * ans["SpBDepl",])
  tot["ExBDepl"] <- sum((ans["effexB0",]/tot["effexB0"]) * ans["ExBDepl",])
  tot["legalDepl"] <- sum(wgts * ans["legalDepl",])
  tot["MSYDepl"] <- sum(wgts * ans["MSYDepl",])
  tot["LML"] <- mean(ans["LML",])
  tot["bLML"] <- sum(wgts * ans["bLML",])
  tot["harvestR"] <- mean(ans["harvestR",])
  tot["catch"] <- sum(ans["catch",])
  ans <- cbind(ans,tot)
  colnames(ans) <- c(paste0("p",1:numpop),"total")
  return(ans)
}  # end of getregionprops

