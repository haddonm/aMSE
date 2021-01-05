

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

#' @title getlistvar extracts a vector or matrix from zoneC
#'
#' @description getlistvar extracts a vector or matrix from zoneC.
#'    If a vector of scalars, the names relate to populations, if a
#'    matrix the columns relate to populations. Only Me, R0, B0, effB0,
#'    ExB0, effExB0, MSY, MSYDepl, bLML, popq, SaM, SAU, popdef,
#'    LML, Maturity, WtL, Emergent, and MatWt are currently valid
#'    choices. The indexvar = popdef would generate a listing of all
#'    the constants. If you only want a single constant from popdefs
#'    then use indexvar2.
#'
#' @param zoneC the constants components of the simulated zone
#' @param indexvar the name of the variable to be extracted; character
#' @param indexvar2 the name of the variable within popdef to extract
#'
#' @return either a vector or matrix of values depending on variable
#' @export
#'
#' @examples
#'  data(zone)
#'  zoneC <- zone$zoneC
#'  getlistvar(zoneC,"MSY")
#'  getlistvar(zoneC,"B0")
#'  getlistvar(zoneC,"popdef","AvRec")
getlistvar <- function(zoneC,indexvar,indexvar2="") {
  if (is.character(indexvar)) {
    fields <- c("Me","R0","B0","effB0","ExB0","effExB0","MSY",
                "MSYDepl","bLML","popq","SaM","SAU",
                "popdef","LML","Maturity","WtL","Emergent","MatWt")
    vects<- c("popdef","LML","Maturity","WtL","Emergent","MatWt")
    pick <- match(indexvar,fields)
    numpop <- length(zoneC)
    if (is.na(pick)) {
      warning(paste0(fields,"  "))
      stop("Invalid variable name attempted in getlistvar")
    }
    x <- sapply(zoneC,"[[",indexvar)
    if (fields[pick] %in% vects) {
      colnames(x) <- paste0("p",1:numpop)
    } else {
      names(x) <- paste0("p",1:numpop)
    }
    if ((indexvar == "popdef") & (nchar(indexvar2) > 0)) {
      x <- sapply(zoneC,"[[",indexvar)[indexvar2,]
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
#'  txtline <- "Depleted, TRUE"
#'  aMSE:::getLogical(txtline,nb=1)
#'  txtline2 <- "calcthis, TRUE, FALSE"
#'  aMSE:::getLogical(txtline2,nb=2)
getLogical <- function(inline,nb) {  #inline <- txtline; nb=2
  tmp <- unlist(strsplit(inline,","))
  tmp <- removeEmpty(tmp)
  outtmp <- as.logical(as.character(tmp[2:(nb+1)]))
  return(outtmp)
}

#' @title getnas gets the numbers-at-size for all spatial scales
#'
#' @description getnas extracts numbers-at-size for all populations,
#'     SAUs, and the zone. There will therefore be numpop + nSAU + 1
#'     columns of numbers-at-size fr the particular year given.
#'
#' @param zoneD the dynamic portion of the zone
#' @param yr which yr from the range available should be summarized
#' @param glob the global variables object
#' @param zone the zone constants
#'
#' @return a matrix of numpop + 1 columns of numbers-at-size for each population
#'     and for the zone
#' @export
#'
#' @examples
#' data(zone)
#' nas <- getnas(zone$zoneD,yr=1,glob=zone$glb,zone=zone$zone1)
#' round(nas[1:30,],1)
getnas <- function(zoneD,yr,glob,zone) {# zoneD=zoneD;yr=1;glob=glb;zone=zone1
  numpop <- glob$numpop
  nSAU <- glob$nSAU
  columns <- c(paste0("p",1:numpop),"zone") # ,zone$SAUnames
  nas <- matrix(0,nrow=glob$Nclass,ncol=length(columns),
                dimnames=list(glob$midpts,columns))
  for (pop in 1:numpop) nas[,pop] <- zoneD$Nt[,yr,pop]
  nas[,"zone"] <- rowSums(nas[,1:numpop],na.rm=TRUE)
  # for (mu in 1:nSAU) {
  #   pickcol <- which(zoneD$SAU == mu)
  #   nas[,(numpop+mu)] <- rowSums(nas[,pickcol],na.rm=TRUE)
  # }
  return(nas)
} #end of getnas

##  inzone <- zone1; indexVar <- "Nt"
## gets all LF data from a zone across pops and years
#' @title getzoneLF extracts all LF data from a zone across pops and years
#'
#' @description getzoneLF extracts all LF data from a zone across
#'     all populations for each year. Thus an Nclass x Nyrs X numpop
#'     matrix is compressed into a Nclass x Nyrs matrix
#'
#' @param zoneD the dynamic part of the zone
#' @param glb the global constants
#'
#' @return an Nclass x Nyrs matrix containing LF data across all
#'     populations by year
#' @export
#'
#' @examples
#' print("An example needs to be written")
getzoneLF <- function(zoneD,glb) { # need to define years
  numpop <- glb$numpop
  Nyrs <- glb$Nyrs
  storeLF <- matrix(0,nrow=glb$Nclass,ncol=Nyrs,
                    dimnames=list(glb$midpts,1:Nyrs))
  for (yr in 1:Nyrs)
    storeLF[,yr] <- rowSums(zoneD$Nt[,yr,])
  return(storeLF)
} # end of getzoneLF



#' @title getsingleNum find a line of text and extracts a single number
#'
#' @description getsingleNum uses grep to fine an input line. The
#'     variable 'context' is used to identify within which function
#'     getsingleNum is being called. This can then be used in a
#'     the stop function if grep fails to find the 'varname' inside the
#'     'intxt', which is usually obtained using readLines.
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
#'  txtlines <- c("replicates, 100","Some_other_text, 52")
#'  getsingleNum("replicates",txtlines,context=context)
#'  getsingleNum("eeplicates",txtlines,context=context)
#'  getsingleNum("other",txtlines,context=context)
#' }
getsingleNum <- function(varname,intxt,context="") {
  begin <- grep(varname,intxt)
  if (length(begin) > 0) {
    return(as.numeric(getConst(intxt[begin],1)))
  } else {
    stop(paste0("No data for ",varname," within ",context))
  }
}

#' @title getsauzone summarizes zoneD into MSU and zone
#'
#' @description getsauzone rowsums the matrices of matureB, exploitB,
#'     catch, and recruit from the input zoneD into the separate SAUs
#'     and the total into the zone. The harvestR is simply the
#'     respective catch divided by the exploitB, and the cpue are the
#'     individual population cpue weighted relative to the catch taken
#'     from each population in either the SAU or the complete zone.
#'
#' @param zoneD the zoneD after nyrs of dynamics
#'
#' @return a list of six matrices of nSAU columns of SAU summaries,
#'     and one column for the zone
#' @export
#'
#' @examples
#' print("wait on an example")
getsauzone <- function(zoneD) { # zoneC=zoneC; zoneD=zoneD; glb=glb
  iSAU <- zoneD$SAU
  SAU <- unique(iSAU)
  nSAU <- length(SAU)
  matB <- getsum(zoneD$matureB,iSAU)
  expB <- getsum(zoneD$exploitB,iSAU)
  catch <- getsum(zoneD$catch,iSAU)
  recruit <- getsum(zoneD$recruit,iSAU)
  harvestR <- catch/expB
  cpue <- catch # just ot have a labelled matrix ready
  wtzone <- zoneD$catch/catch[,(nSAU+1)]
  wtsau <- zoneD$catch
  for (mu in 1:nSAU) {
    wtsau[,(iSAU==mu)] <- zoneD$catch[,(iSAU==mu)]/catch[,mu]
    cpue[,mu] <- rowSums(zoneD$cpue[,(iSAU==mu)] * wtsau[,(iSAU==mu)])
  }
  cpue[,(nSAU+1)] <- rowSums(zoneD$cpue * wtzone)
  ans <- list(matB=matB,expB=expB,catch=catch,recruit=recruit,
              harvestR=harvestR,cpue=cpue)
  return(ans)
}  # end of getsauzone



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

#' @title getsum sums each of the main dynamics within zoneD
#'
#' @description getsum is only used by getsauzone to sum each of the main
#'     dynamic variables within zoneD, and hence is not exported.
#'
#' @param inmat what matrix within zoneD to sum into SAU and zone
#' @param index a vector containing an index of populations within SAU
#'
#' @return an nSAU+1 column matrix summarizing each SAU and the zone
#'
#' @examples
#' print("wait on an example")
getsum <- function(inmat,index) {
  nSAU <- length(unique(index))
  nyr <- nrow(inmat)
  matO <- matrix(0,nrow=nyr,ncol=(nSAU+1),
                 dimnames=list(1:nyr,c(paste0("SAU",1:nSAU),"zone")))
  for (mu in 1:nSAU) matO[,mu] <- rowSums(inmat[,(index==mu)])
  matO[,(nSAU+1)] <- rowSums(inmat)
  return(matO)
}

#' @title getunFished - extracts all data relating to year 1; unfished.
#'
#' @description getunFished - extracts all data relating to year 1;
#'    unfished. Of course, this will only work as long as movezoneYear
#'    hasn't been called, which would disrupt the contents of the first
#'    year. But if applied to zone rather than zone1 etc, this should be
#'    fine. The outputs include: R0, B0, B0crypt, popdef, MSY,
#'    MSYDepl, LML, ExB0, and Nt
#'
#' @param zoneC the constant part of the zone
#' @param zoneD the dynamic part of the zone
#' @param glb contains the global variables used everywhere
#'
#' @return a list of multiple components relating to the unfished stock this
#'     includes a list of each of the equilibrium populations and the zone as
#'     the last element of the list
#' @export
#' @examples
#' \dontrun{
#'   data(zone) # would normally use zone <- makeequilzone(resdir,"control.csv")
#'   unfish <- getunFished(zone$zoneC,zone$zoneD,zone$glb)
#'   str(unfish,max.level=1)
#' }
getunFished <- function(zoneC,zoneD,glb) {  # inzone=zone
  if((zoneC[[1]]$B0 - zoneD$matureB[1,1]) > 0.001) {
    warning("Unfished characters taken from a depleted zone in getunFished \n")
  }
  numpop <- glb$numpop
  nofished <- vector("list",(numpop+1))
  label <- c("R0","B0","B0crypt","popdef","MSY","MSYDepl","LML","ExB0",
             "Nt")
  for (pop in 1:numpop) {
    Nt <- zoneD$Nt[,1,pop]
    R0 <- zoneC[[pop]]$R0
    B0 <- zoneC[[pop]]$B0
    B0crypt <- B0 -
            sum(zoneC[[pop]]$Emergent * Nt * zoneC[[pop]]$MatWt)/1e06
    popdef <- zoneC[[pop]]$popdef
    MSY <- zoneC[[pop]]$MSY
    MSYDepl <- zoneC[[pop]]$MSYDepl
    LML <- zoneC[[pop]]$LML[1]
    ExB0 <- zoneC[[pop]]$ExB0
    ans <- list(R0,B0,B0crypt,popdef,MSY,MSYDepl,LML,ExB0,Nt)
    names(ans) <- label
    nofished[[pop]] <- ans
  }
  R0 <- 0
  B0 <- 0
  B0crypt <- 0
  MSY <- 0
  ExB0 <- 0
  Nt <- numeric(glb$Nclass)
  for (pop in 1:numpop) { #sums across all pops to go into place numpop+1
    R0 <- R0 + nofished[[pop]]$R0
    B0 <- B0 + nofished[[pop]]$B0  ## or use: sum(sapply(zone,"[[","B0"))
    B0crypt <- B0crypt + nofished[[pop]]$B0crypt
    MSY <- MSY + nofished[[pop]]$MSY
    ExB0 <- ExB0 + nofished[[pop]]$ExB0
    Nt <- Nt + nofished[[pop]]$Nt
  }
  zoneT <- list(R0,B0,B0crypt,MSY,ExB0,Nt)
  names(zoneT) <- c("R0","B0","B0crypt","MSY","ExB0","Nt")
  nofished[[numpop+1]] <- zoneT
  return(nofished)
}  # End of getUnfished


#' @title getzoneprod zone scale summary of product matrix from doproduction
#'
#' @description getzoneprod takes in the product matrix from doproduction and
#'     sums the mature and exploitable biomassand the total catches. Then it
#'     recalculates, for the zone, the approximate annual harvest rate, the
#'     depletion level, and the relative catch rate. The annual harvest rate and
#'     relative catch rate are simply the arithmetic mean of the values for all
#'     populations, which the depletion is the column of mature biomass divided
#'     by the first value = B0
#'
#' @param product The productivity array from doproduction containing the
#'     range of imposed harvest rates, and the resulting outputs for each
#'     population
#'
#' @return a matrix containing the approximate productivity matrix for the zone
#' @export
#'
#' @examples
#' data(zone)
#' zoneprod <- getzoneprod(zone$product)
#' head(zoneprod,20)
getzoneprod <- function(product) {
  numrow <- dim(product)[1]
  rows <- rownames(product[,,1])
  columns <- colnames(product[,,1])
  numpop <- dim(product)[3]
  zoneprod <- matrix(0,nrow=numrow,ncol=length(columns),
                     dimnames=list(rows,columns))
  for (i in 1:numpop) zoneprod <- zoneprod + product[,,i]
  #head(zoneprod,20)
  zoneprod[,"AnnH"] <- apply(product[,"AnnH",],1,mean,na.rm=TRUE)
  zoneprod[,"Deplet"] <- zoneprod[,"MatB"]/zoneprod[1,"MatB"]
  zoneprod[,"RelCE"] <- apply(product[,"RelCE",],1,mean,na.rm=TRUE)
  return(zoneprod)
} # end of getzoneprod

#' @title getzoneprops extracts the depletion level for a given year
#'
#' @description getzoneprops extracts a set of properties for each
#'     population and summarizes them for the zone as well. These
#'     properties include effB0, matureB, legalmatB, propprot, MSY,
#'     effexB0, exploitB, SpBDepl, ExBDepl, legalDepl, MSYDepl, LML,
#'     and bLML. See the model documentation for the meaning of each
#'     of these.
#'
#' @param zoneC the constants for the zone being simulated
#' @param zoneD the dynamic part of the zone being simulated
#' @param glb the global constants
#' @param year which years information is wanted, default=1
#'
#' @return a matrix of population properties with a column of totals
#' @export
#'
#' @examples
#' data(zone)
#' round(getzoneprops(zone$zoneC,zone$zoneD,zone$glb),4)
getzoneprops <- function(zoneC,zoneD,glb,year=1) { #zoneC=zoneC; zoneD=zoneDD;glb=zone$glb; year=1
  numpop <- glb$numpop
  Nclass <- glb$Nclass
  sau <- getvar(zoneC,"SAU")
  B0 <- getvar(zoneC,"B0")
  ExB0 <- getvar(zoneC,"ExB0")
  blml <- getvar(zoneC,"bLML")
  msy <- getvar(zoneC,"MSY")
  msydepl <- getvar(zoneC,"MSYDepl")
  matB <- zoneD$matureB[year,]
  ExB <- zoneD$exploitB[year,]
  harvestR <- zoneD$harvestR[year,]
  catch <- zoneD$catch[year,]
  legalmatB <- numeric(numpop); names(legalmatB) <- 1:numpop
  lml <- sapply(zoneC,"[[","LML")[year,]
  for (pop in 1:numpop) { #   pop=1
    MatWt <- zoneC[[pop]]$MatWt
    first <- which(glb$midpts >= lml[pop])[1]
    legalmatB[pop] <- sum(zoneC[[pop]]$MatWt[first:Nclass] *
                            zoneD$Nt[first:Nclass,year,pop])/1000000.0
  }
  deplet <- matB/B0
  depletEx <- ExB/ExB0
  legaldepl <- legalmatB/B0
  propprot <- (matB - legalmatB)/matB
  label <- c("SAU","B0","matureB","legalmatB","propprot","MSY",
             "exB0","exploitB","SpBDepl","ExBDepl",
             "legalDepl","MSYDepl","LML","bLML","harvestR","catch")
  ans <- rbind(sau,B0,matB,legalmatB,propprot,msy,ExB0,ExB,deplet,
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
  colnames(ans) <- c(paste0("p",1:numpop),"zone")
  return(ans)
}  # end of getzoneprops

#' @title getvar a replacement for sapply to obtain scalar constants
#'
#' @description getvar is a replacement for sapply to obtain scalar
#'     constants from zoneC and is significantly faster. It should
#'     be used to obtain things like B0, R0, MSY, popq, etc. Still
#'     need to use sapply to pull out vectors.
#'
#' @param zoneC the constants object for the zone
#' @param invar a character variable eg. "B0" or "R0"
#'
#' @return a numpop vector of the invar constants from zoneC
#' @export
#'
#' @seealso getlistvar
#'
#' @examples
#' data(zone)
#' zoneC <- zone$zoneC
#' getvar(zoneC,"MSY")
#' getvar(zoneC,"B0")
getvar <- function(zoneC,invar) {
  npop <- length(zoneC)
  recs <- numeric(npop)
  for (i in 1:npop) recs[i] <- zoneC[[i]][[invar]]
  return(recs)
} # end of getvar

#' @title getvect extracts invar from the popdef vector in zoneC
#'
#' @description getvect extracts a numpop vector of invar from the
#'     popdef vector in zoneC. Still need to use sapply to pull out
#'     complete vectors such as popdef or maturity etc.
#'
#' @param zoneC the constants object for the zone
#' @param invar a character variable eg. "steeph", "DLMax"
#'
#' @return a numpop vector of invar from the numpop popdefs in zoneC
#' @export
#'
#' @seealso getlistvar
#'
#' @examples
#' data(zone)
#' zoneC <- zone$zoneC
#' getvect(zoneC,"steeph")
getvect <- function(zoneC,invar) {
  npop <- length(zoneC)
  ans <- numeric(npop)
  for (i in 1:npop) ans[i] <- zoneC[[i]]$popdef[invar]
  return(ans)
} # end of getvect
