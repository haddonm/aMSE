
#' @title getConst extracts 'nb' numbers from a line of text
#'
#' @description getConst parses a line of text and extracts 'nb' pieces of
#'     text as numbers
#'
#' @param inline text line to be parsed, usually obtained using readLines
#' @param nb the number of numbers to extract
#' @param index which non-empty object to begin extracting from?
#'
#' @return a vector of length 'nb'
#' @export
#'
#' @examples
#'   txtline <- "MaxDL , 32,32,32"
#'   getConst(txtline,nb=3,index=2)
getConst <- function(inline,nb,index=2) { # parses lines containing numbers
  ans <- numeric(nb)
  tmp <- unlist(strsplit(inline,","))
  if (length(tmp) < (nb+1))
    warning(paste("possible problem with data",tmp[1],
                  "missing comma?",sep=" "),"\n")
  count <- 0
  for (j in index:(nb+index-1)) {
    count <- count + 1
    ans[count] <- as.numeric(tmp[j])
  }
  return(ans)
}   # end getConst



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
#'    ExB0, effExB0, MSY, MSYDepl, bLML, popq, SaM, SMU, popdef,
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
                "MSYDepl","bLML","popq","SaM","SMU",
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

#' @title getLogical extracts nb logicals from an input line of text
#'
#' @description getLogical obtains nb logicals from an input line
#'
#' @param inline text line to be parsed, usually obtained using readLines
#' @param nb the number of logicals to extract, if nb is longer than the
#'     number of logicals within inline the vector will contain NAs
#'
#' @return a vector of length nb
#' @export
#'
#' @examples
#' \dontrun{
#'  # Not exported, prefix with AbMSE:::
#'  txtline <- "Depleted, TRUE"
#'  AbMSE:::getLogical(txtline,nb=1)
#'  txtline2 <- "calcthis, TRUE, FALSE"
#'  AbMSE:::getLogical(txtline2,nb=2)
#' }
getLogical <- function(inline,nb) {  #inline <- txtline; nb=2
  tmp <- unlist(strsplit(inline,","))
  tmp <- removeEmpty(tmp)
  outtmp <- as.logical(as.character(tmp[2:(nb+1)]))
  return(outtmp)
}

#' @title getnas gets the numbers-at-size for all spatial scales
#'
#' @description getnas extracts numbers-at-size for all populations,
#'     SMUs, and the region. There will therefore be numpop + nSMU + 1
#'     columns of numbers-at-size fr the particular year given.
#'
#' @param regD the dynamic portion of the region
#' @param yr which yr from the range available should be summarized
#' @param glob the global variables object
#' @param region the region constants
#'
#' @return a matrix of numpop + nSMU + 1 columns of numbers-at-size
#' @export
#'
#' @examples
#' data(region1)
#' glb <- region1$globals
#' data(testregD)
#' nas <- getnas(testregD,yr=1,glob=glb,region=region1)
#' round(nas[1:30,],1)
getnas <- function(regD,yr,glob,region) {
  numpop <- glob$numpop
  nSMU <- glob$nSMU
  columns <- c(paste0("p",1:numpop),region$SMUnames,"region")
  nas <- matrix(0,nrow=glob$Nclass,ncol=length(columns),
                dimnames=list(glob$midpts,columns))
  for (pop in 1:numpop) nas[,pop] <- regD$Nt[,yr,pop]
  nas[,"region"] <- rowSums(nas[,1:numpop],na.rm=TRUE)
  for (mu in 1:nSMU) {
    pickcol <- which(regD$SMU == mu)
    nas[,(numpop+mu)] <- rowSums(nas[,pickcol],na.rm=TRUE)
  }
  return(nas)
} #end of getnas

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

#' @title getsingleNum extracts a single number from an input line of text
#'
#' @description getsingleNum obtains a number from an input line. The
#'     variable context is the a variable declared in each function within
#'     which getsingleNum is to be called. This cvan then be used in a
#'     the stop function.
#'
#' @param varname the name of the variable to get from intxt
#' @param intxt text to be parsed, usually obtained using readLines
#' @param context the function name within which getsingleNum is
#'     being used, defaults to empty string
#'
#' @return a single number
#' @export
#'
#' @examples
#' \dontrun{
#'  # Not exported, prefix with aMSE:::
#'  context = "Function Example"
#'  txtline <- "replicates, 100"
#'  aMSE:::getsingleNum("replicates",txtline,context=context)
#'  aMSE:::getsingleNum("eeplicates",txtline,context=context)
#' }
getsingleNum <- function(varname,intxt,context="") {
  begin <- grep(varname,intxt)
  if (length(begin) > 0) {
    return(as.numeric(getConst(intxt[begin],1)))
  } else {
    stop(paste0("No data for ",varname," within ",context))
  }
}

#' @title getStr obtains a string from an input text line
#'
#' @description  getStr obtains a string from an input text line in
#'     which any parts are separated by ','. Then, after ignoring the
#'     first component, assumed to be a label, it returns the first
#'     nb parts.
#'
#' @param inline input text line with components separated by ','
#' @param nb number of parts to return
#'
#' @return a vector of character string(s)
#' @export
#'
#' @examples
#'   txt <- "runlabel, development_run, label for this particular run"
#'   getStr(txt,1)
getStr <- function(inline,nb) {
  tmp <- unlist(strsplit(inline,","))
  tmp <- removeEmpty(tmp)
  outconst <- as.character(tmp[2:(nb+1)])
  return(outconst)
} # end of getStr



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
  B0 <- sapply(regC,"[[","B0")
  ExB0 <- sapply(regC,"[[","ExB0")
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
  deplet <- matB/B0
  depletEx <- ExB/ExB0
  legaldepl <- legalmatB/B0
  propprot <- (matB - legalmatB)/matB
  label <- c("B0","matureB","legalmatB","propprot","MSY",
             "exB0","exploitB","SpBDepl","ExBDepl",
             "legalDepl","MSYDepl","LML","bLML","harvestR","catch")
  ans <- rbind(B0,matB,legalmatB,propprot,msy,ExB0,ExB,deplet,
               depletEx,legaldepl,msydepl,lml,blml,harvestR,catch)
  rownames(ans) <- label
  tot <- numeric(length(rownames(ans))); names(tot) <- label
  tot[c(1:3,5:7)] <- rowSums(ans)[c(1:3,5:7)]
  tot["propprot"] <- (tot["matureB"] - tot["legalmatB"])/tot["matureB"]
  wgts <- ans["B0",]/tot["B0"]
  tot["SpBDepl"] <- sum(wgts * ans["SpBDepl",])
  tot["ExBDepl"] <- sum((ans["exB0",]/tot["exB0"]) * ans["ExBDepl",])
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

#' @title getvar a replacement for sapply to obtain scalar constants
#'
#' @description getvar is a replacement for sapply to obtain scalar
#'     constants from regionC and is significantly faster. It should
#'     be used to obtain things like B0, R0, MSY, popq, etc. Still
#'     need to use sapply to pull out vectors.
#'
#' @param regC the constants object for the region
#' @param invar a character variable eg. "B0" or "R0"
#'
#' @return a numpop vector of the invar constants from regionC
#' @export
#'
#' @examples
#' data(testregC)
#' getvar(testregC,"MSY")
#' getvar(testregC,"B0")
getvar <- function(regC,invar) {
  npop <- length(regC)
  recs <- numeric(npop)
  for (i in 1:npop) recs[i] <- regC[[i]][[invar]]
  return(recs)
} # end of getvar

#' @title getvect extracts invar from the popdef vector in regionC
#'
#' @description getvect extracts a numpop vector of invar from the
#'     popdef vector in regionC. Still need to use sapply to pull out
#'     complete vectors such as popdef or maturity etc.
#'
#' @param regC the constants object for the region
#' @param invar a character variable eg. "steeph", "DLMax"
#'
#' @return a numpop vector of invar from the numpop popdefs in regionC
#' @export
#'
#' @examples
#' data(testregC)
#' getvect(testregC,"steeph")
getvect <- function(regC,invar) {
  npop <- length(regC)
  ans <- numeric(npop)
  for (i in 1:npop) ans[i] <- regC[[i]]$popdef[invar]
  return(ans)
} # end of getvect
