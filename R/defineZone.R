#  hutils::listFunctions("C:/Users/User/Dropbox/A_Code/aMSE/R/defineZone.R")

#' @title defineBlock subdivides defined number of populations into blocks
#'
#' @description subdivides the defined number of populations into blocks
#'     in the numbers defined in the constant blkpop read into 'constants'.
#'     It allocates the populations sequentially to numblk blocks; if
#'     numblk does not divide into numpop exactly, any remainder is put
#'     into the final block. For block can read SAU = spatial assessment
#'     unit.
#'
#' @param numblk number of blocks in the simulated zone, from 'constants'
#' @param blknum number of populations in each block, from 'constants'
#' @param numpop the total number of populations in the zone
#'
#' @return a vector of length numpop, indexing which block each pop is in.
#' @export
#'
#' @examples
#'   nblock <- 4
#'   blkpop <- c(3,4,8,2)
#'   blockI <- defineBlock(nblock,blkpop,numpop=17)
#'   print(blockI)
#'   cat(length(blockI),sum(blkpop),"\n")
defineBlock <- function(numblk,blknum,numpop) {
  index <- numeric(numpop)
  if (sum(blknum) != numpop) warning("Wrong numpop or block definitions")
  if (numblk == 1) {
    index <- rep(1,numpop)
  } else {
    begin <- 1
    finish <- blknum[1]
    for (i in 1:(numblk-1)) {
      index[begin:finish] <- rep(i,blknum[i])
      begin <- finish + 1
      finish <- finish + blknum[i+1]
    }
    index[begin:numpop] <- rep(numblk,blknum[i+1])
  }
  return(index)
} # end of defineBlock

#' @title definepops makes a parameter vector, popdef, for each population
#'
#' @description definepops makes a parameter vector, popdef, for each
#'     population. Each population will have slightly different
#'     biological properties sampled from probability density
#'     distributions defined in const=condDat$constants. The input can
#'     be either the whole of condDat or just condDat$constants, which
#'     contains the definitions of the probability density functions
#'     used to define each population.
#'
#' @param inSAU the number of SAUs defined in the zone, a scalar
#' @param inSAUindex a vector of the block index for each population
#' @param const a matrix with rows containing the PDF parameters for
#'     each of the biological parameters. The columns of the matrix
#'     reflect values for each populations from readdatafile
#' @param glob the globals object from readzonefile
#'
#' @return a matrix with a row for each population and whose columns
#'     are parameters defining the biological parameters for each population.
#' @export
#'
#' @examples
#' \dontrun{
#'   data(zone1)
#'   glb <- zone1$globals
#'   data(constants)
#'   ans <- makezoneC(zone1,constants)
#'   zoneC <- ans$zoneC
#'   popdefs <- ans$popdefs
#'   print(popdefs)
#'  }
definepops <- function(inSAU,inSAUindex,const,glob) {
  #  inSAU=nSAU; inSAUindex=SAUindex; const=saudat; glob=glb
  numpop <- glob$numpop
  hyrs <- glob$hyrs
  columns <- c("DLMax","L50","L95","SigMax","SaMa","SaMb","Wta","Wtb","Me",
               "L50C","deltaC","AvRec","SelP1","SelP2","Nyrs","steeph",
               "MaxCE","L50mat","SAU","lambda","qest","scalece")
  popdefs <- matrix(0,nrow=numpop,ncol=length(columns),
                    dimnames=list(1:numpop,columns))
  popdefs[,"Nyrs"] <- rep(hyrs,numpop) # Num Years - for selectivity?
  for (pop in 1:numpop) {  # pop=1
    popdefs[pop,"Me"] <- rnormz(1,mean=const["Me",pop],
                               sd=const["sMe",pop])
    popdefs[pop,"DLMax"] <- rnormz(1,mean=const["DLMax",pop],
                                  sd=const["sDLMax",pop])   # DLMax
    popdefs[pop,"L50"] <- rnormz(1,mean=const["L50",pop],
                                const["sL50",pop])
    popdefs[pop,"L95"] <- popdefs[pop,"L50"] +
      rnormz(1,mean=const["L50inc",pop],const["sL50inc",pop]) #L95
    popdefs[pop,"SigMax"] <- rnormz(1,mean=const["SigMax",pop],
                                   sd=const["sSigMax",pop])   # SigMax
    popdefs[pop,"L50mat"] <- rnormz(1,mean=const["L50Mat",pop],
                                   const["sL50Mat",pop])   # L50Mat
    popdefs[pop,"SaMa"] <- const["SaMa",pop]
    popdefs[pop,"SaMb"] <- -popdefs[pop,"SaMa"]/popdefs[pop,"L50mat"]
    popdefs[pop,"SelP1"] <- const["selL50p",pop]
    popdefs[pop,"SelP2"] <-  const["selL95p",pop]
    popdefs[pop,"L50C"] <- rnormz(1,mean=const["L50C",pop],
                                 const["sL50C",pop])
    popdefs[pop,"deltaC"] <-rnormz(1,mean=const["deltaC",pop],
                                  const["sdeltaC",pop])
    popdefs[pop,"Wtb"] <- rnormz(1,mean=const["Wtb",pop],const["sWtb",pop])
    # popdefs[pop,"Wta"] <- const["Wtbtoa",pop] *
    #   popdefs[pop,"Wtb"]^const["sWtbtoa",pop]
    popdefs[pop,"Wta"] <- rnormz(1,mean=const["Wta",pop],const["sWta",pop])
    popdefs[pop,"steeph"] <- rnormz(1,mean=const["defsteep",pop],
                                   const["sdefsteep",pop])
    popdefs[pop,"AvRec"] <- rlnormz(1,meanlog=const["AvRec",pop],
                                   const["sAvRec",pop])
    popdefs[pop,"MaxCE"] <- 0 # gets value in makezone
    popdefs[pop,"SAU"] <- const["SAU",pop]
    popdefs[pop,"lambda"] <- const["lambda",pop]
    popdefs[pop,"qest"] <-const["qest",pop]
    popdefs[pop,"scalece"] <- 0 # gets value in makezone
  }
  test <- popdefs[,"L95"] - popdefs[,"L50"]
  if (any(test < 0)) stop("L95 < L50 - adjust the input data file  \n")
  return(popdefs)
} # end of definepops

#' @title doproduction estimates a production curve for each population
#'
#' @description doproduction estimates a production curve for each
#'     population in the simulated zone. It does this by sequentially
#'     applying a series of increasing harvest rates and running the
#'     populations to equilibrium to discover, empirically, the yield
#'     and other details, such as exploitable biomass, mature biomass,
#'     the actual harvest rate, the catch or yield, the mature
#'     depletion, and the relative cpue. Keep in mind that the time
#'     taken to run this function depends on the number of populations
#'     but mainly on the length of the harvest rate sequence, which is
#'     determined by lowlim, uplim, and inc. The lonjger the sequence
#'     the slower the run. However, the estimates of MSY and related
#'     statistics are expected to have better accuracy (resolution)
#'     the finer the harvest rate increments. So when initiating the
#'     zone it is best to have a long sequence of finely incremented
#'     harvest rates that take a long time.
#'
#' @param zoneC the constants components of the simulated zone
#' @param zoneD the dynamic components of the simulated zone
#' @param glob the general global variables
#' @param lowlim the lower limit of harvest rate applied, default=0.0
#' @param uplim the upper limit of harvest rate applied, default=0.35
#' @param inc the harvest rate increment at each step, default=0.005
#'
#' @return an array of the six productivity variables by numpop by
#'     the number of harvest rates applied
#' @export
#'
#' @examples
#' print("wait")
#' #  zoneC=zoneC; zoneD=zoneD; glob=glb; lowlim=0.0;uplim=0.4;inc=0.01
doproduction <- function(zoneC,zoneD,glob,lowlim=0.0,uplim=0.35,inc=0.005) {
  numpop <- glob$numpop
  Nclass <- glob$Nclass
 # hyrs <- glob$hyrs
  larvdisp <- glob$larvdisp
  initH <- seq(lowlim,uplim,inc)
  nH <- length(initH)
  columns <- c("ExB","MatB","AnnH","Catch","Deplet","RelCE")
  results <- array(0,dim=c(nH,6,numpop),dimnames=list(initH,columns,1:numpop))
  for (aH in 1:nH) { # aH=1 ; yr=2
    zoneP <- runthreeH(zoneC=zoneC,zoneD,inHarv=rep(initH[aH],numpop),glob)
    results[aH,"ExB",] <- zoneP$exploitB[1,]
    results[aH,"MatB",] <- zoneP$matureB[1,]
    results[aH,"AnnH",] <- zoneP$harvestR[1,]
    results[aH,"Catch",] <- zoneP$catch[1,]
    results[aH,"Deplet",] <- zoneP$deplsB[1,]
    results[aH,"RelCE",] <- zoneP$cpue[1,]
  } # end of yr loop
  return(results)
} # end of doproduction

#' @title driftrec adjusts the recruitment allowing for larval drift
#'
#' @description driftrec adjusts the recruitment levels in the stock
#'     dynamics allowing for larval drift. The input drift levels can
#'     be a constant for each population, or a vector, with a value
#'     for each population if their rates are assumed to different. If
#'     a constant valuer is input it is extended to a vector of the
#'     same length as the number of populations. Each population is
#'     assumed to connect only to the populations either side of it.
#'     The populations on the end only connect into the simulated
#'     populations and so only lose half the recp value for that
#'     population, which is sent to the population it is beside in the
#'     conditioning. If the zone is not conditioned using yield or
#'     productivity from spatial data then productivity will be
#'     allocated randomly and the recp (recruitment proportion lost)
#'     may as well be set to zero.
#'
#' @param recs the vector of recruits calculated internally to oneyear
#' @param recp the recruitment proporiton lost. Either a constant (eg
#'     0.02 for 2 percent) or a vector, a value for each population.
#'
#' @return a vector of revised recruitment levels accounting for
#'     larval drift.
#'
#' @export
#'
#' @examples
#'   recs <- c(243699,49747,1285137,1492742,923282,273427)
#'   newrec <- driftrec(recs,0.04)
#'   newrec0 <- driftrec(recs,0.0)
#'   recs
#'   newrec
#'   newrec0
#'   sum(recs)
#'   sum(newrec)
#'   sum(newrec0)
driftrec <- function(recs,recp) {
  nrec <- length(recs)
  if (length(recp) == 1) recp <- rep(recp,nrec)
  addr <- recp/2.0
  subr <- 1 - addr
  subrec <- 1 - recp
  newrec <- numeric(nrec)
  newrec[1] <- (recs[1] * subr[1]) + (recs[2] * addr[2])
  newrec[nrec] <- (recs[nrec] * subr[nrec]) +
    (recs[nrec-1] * addr[nrec-1])
  for (i in 2:(nrec-1))
    newrec[i] <- (recs[i-1] * addr[i-1]) + (recs[i] * subrec[i]) +
    (recs[i+1] * addr[i+1])
  return(newrec)
} # end of driftrec an internal function


#' @title fillzoneDef characterizes simulated Zone; holds zone properties
#'
#' @description fillzoneDef Characterizes the simulated Zone; holds the
#'     zone properties. Lists the date and time the zone was made, its size
#'     and basic biological properties: MSY, B0, Average bLML, Average SaM,
#'     and more details
#'
#' @param zoneC the constnats for the zone to be characterized.
#' @param zoneD the dynamic parts of the rgion being characterized
#' @param prod the production object from doproduction
#'
#' @return A list of 10 objects containing the properties of the zone;
#'     has the class zoneDefinition - for S3 methods
#'
#' \itemize{
#'   \item Production the productivity of each populations separately
#'   \item defpop each population's defining parameters
#'   \item Struct the zone's hierarchical structure
#'   \item nBlock the number of blocks
#'   \item numpop the number of populations
#'   \item SummaryPop biological properties of populations
#'   \item SummaryZone summary of the zone's properties
#'   \item BlockProp summary of the block properties
#'   \item Date date/time of fillzoneDef activation; if the function is run
#'      soon after a zone's creation this can be a proxy for the zone
#'      creation.
#' }
#' @export
#'
#' @examples
#' # needs summaryPop and summaryZone
#' txt2 <- 'always use examples rather than example'
fillzoneDef <- function(zoneC,zoneD,prod) {  # inzone=zone; prod=production
   defNames <- c("production","defpop","struct","nBlock","numpop",
                 "summaryPop","summaryZone","blockProp","date")
   zDef <- vector("list",length(defNames))
   names(zDef) <- defNames
   numpop <- length(zoneC)
   pops <- seq(1,numpop,1)
   popdefs <- t(sapply(zoneC,"[[","popdef"))
   blockI <- popdefs[,"SAU"]
   nblock <- length(unique(blockI))
   struct <- matrix(0,nrow=numpop,ncol=3,
                    dimnames=list(pops,c("Blocks","Population","Hexagon")))
   struct[,1] <- blockI
   struct[,2] <- pops
   ans <- 5 #zoneProperty(inzone)
   #ndef <- length(zone[[1]]$popdef)
   #defname <- names(zone[[1]]$popdef)

   zDef[["production"]] <- prod
   zDef[["defpop"]] <- popdefs  # already defined
   zDef[["struct"]] <- struct
   zDef[["nBlock"]] <- nblock
   zDef[["numpop"]] <- numpop
   zDef[["summaryPop"]] <- 5 #ans$summaryMatrix
   zDef[["summaryZone"]] <- 5 #ans$ZoneSummary
   zDef[["blockProp"]] <- 5 #summaryBlock(inzone)
   zDef[["date"]] <- date()
   class(zDef) <- "zoneDefinition"
   return(zDef)
}  # end of fillzoneDef

#' @title findF1 approximates the F0.1 harvest rate for each population
#'
#' @description findF1 generates estimates of the annual harvest rate
#'     that would be equivalent to the F0.1 for each population. It
#'     does this by calculating the gradient of the production curve
#'     for the selected pop number and selecting the best approximation
#'     by searching for the gradient that is cloest to 0.1.
#'
#' @param product the 3D array from doproduction
#'
#' @return either the index within the production array or the vector
#'     of gradients
#' @export
#'
#' @examples
#' \dontrun{
#' # data files changed, no longer works
#' data(product) #for pop=1, a 28percent drop from Hmsy leads to a
#' findF1(product=product) # loss of 3 tonnes, 4 percent of MSY
#' findmsy(product)  # compare the AnnH, Deplet, and RelCE levels.
#' }
findF1 <- function(product) {
  npop <- dim(product)[3]
  label <- c(colnames(product),"index")
  xval <- matrix(0,nrow=npop,ncol=length(label),
                 dimnames=list(1:npop,label))
  for (pop in 1:npop) {
    harv <- product[,"AnnH",pop]
    nH <- length(harv)
    catch <- product[,"Catch",pop]
    grad <- numeric(nH-1)
    for (i in 1:(nH-1)) {
      divisor <- harv[i+1] - harv[i]
      numerator <- catch[i+1] - catch[i]
      grad[i] <- numerator/divisor
    }
    pickF1 <-  which.closest(0.1,grad/grad[1])
    xval[pop,] <- c(product[pickF1,,pop],pickF1)
  }
  return(xval)
} # end of findF1

#' @title findmsy identifies the closest productivity value to MSY
#'
#' @description findmsy for each population in the zone, identifies
#'     the closest productivity value to MSY and returns the vector
#'     of productivity values for the selected harvest rate. The Catch
#'     in each populaiton = MSY and the other variables relate to the
#'     value at MSY
#'
#' @param product array of productivity values output from doproduction
#'
#' @return a matrix of numpop rows containing the approximate MSY
#'     related productivity values and the index of the approx MSY
#' @export
#'
#' @examples
#' \dontrun{  # takes far too long
#'   data(zone1)
#'   glb <- zone1$globals
#'   data(constants)
#'   ans <- makezoneC(zone1,constants)
#'   zoneC <- ans$zoneC
#'   glb <- ans$glb
#'   ans <- makezone(glb,zoneC)
#'   zoneC <- ans$zoneC
#'   zoneD <- ans$zoneD
#'   ans2 <- modzoneC(zoneC,zoneD,glb)
#'   zoneC <- ans2$zoneC  # zone constants
#'   product <- ans2$product
#'   approxMSY <- findmsy(product)
#'   print(approxMSY)
#' }
findmsy <- function(product) {  # product=production
  catch <- product[,"Catch",]
  numpop <- ncol(catch)
  label <- c(colnames(product),"index")
  xval <- matrix(0,nrow=numpop,ncol=length(label),
                 dimnames=list(1:numpop,label))
  for (pop in 1:numpop) { # pop=1
    pick <- which.max(catch[,pop])
    xval[pop,] <- c(product[pick,,pop],pick)
  }
  return(xval)
} # end of findmsy


#' @title logistic a Logistic selectivity function
#'
#' @description logistic a logistic selectivity function. This uses the
#'     logistic function 1/(1+exp(-log(19.0)*(lens-inL50)/(delta))),
#'     where delta = inL95 - inL50. This explicitly defines the SM50
#'     but uses delta (which is SM95-SM50) as the second. This ensures
#'     that when adding variation to parameters, to vary between
#'     populations, when SM95 and SM50 are close together it is not
#'     possible for SM50 to become larger than SM95. Be careful using the
#'     knifeedge option. Strictly knifeedge selectivity would entail the
#'     selectivity values being zero up to the knife-edge and then being 1.0.
#'     This is not what happens here. Instead the knifeedge option literally
#'     sets all values to zero at and below the value of knifeedge but leaves
#'     any curve above that value as it is. Hence this is not strict knife-edge
#'     selectivity. However, it does provide a selectivity curve that reflects
#'     the selectivity of a diver led fleet working on abalone.
#'
#' @param inL50 is the length at 50 percent selection
#' @param delta is the difference between the 95 percent selection and
#'     the 50 percent selection
#' @param lens a vector of lengths for which the logistic maturity value
#'     will be calculated
#' @param knifeedge defaults to 0. If knifeedge is set to a particular
#'     length then the selectivity less than the value of knifeedge is set
#'     to zero.
#' @param maxLML default = 0. The parameter is to allow for a slot selectivity
#'
#' @return A vector of length(lens) containing the predicted selectivity at
#'    length values
#' @export
#'
#' @examples
#' \dontrun{
#' inL50 <- 100.0
#' delta <- 8.0
#' lens <- seq(2,210,2)
#' select <- logistic(inL50,delta,lens)
#' select <- logistic(inL50,delta,lens,knifeedge=105)
#' select <- logistic(inL50,delta,lens,maxLML=185)
#' }
logistic <- function(inL50,delta,lens,knifeedge=0,maxLML=0) {
   ans <- 1/(1+exp(-log(19.0)*(lens-inL50)/(delta)))
   if (knifeedge > 0) {
      pick <- which(lens < knifeedge)
      if (length(pick) > 0) ans[pick] <- 0.0
   }
   if (maxLML > 0) {
     pick <- which(lens > maxLML)
     if (length(pick) > 0) ans[pick] <- 0.0
   }
   return(ans)
} # end of logistic

#' @title makeabpop generates full population structure in an unfished state
#'
#' @description makeabpop generates the full population structure in an
#'    unfished state. The structure is pre-determined: Me R0 B0 effB0
#'    ExB0, effExB0 MSY MSYDepl bLML scalece SaM popdef (vector of constants)
#'    LML G Maturity WtL Emergent Select SelWt MatWt SAUname. See the
#'    AbMSE documentation to see the full definition of a zone, made
#'    up of a zoneC and a zoneD. Notice the presence of effB0 and
#'    effExB0, these relate to the influence of larval dispersal on each
#'    populations productivity. The effective B0 relates to the unfished
#'    mature biomass after larval dispersal occurs and the population
#'    achieves equilibrium.
#'
#' @param popparam the vector of biological parameters values that define
#'     the specific properties of this population. Obtained from
#'     popdefs, which is produced by definepops
#' @param midpts the center values of the size classes dervied from
#'     the zone data file
#' @param initLML the LML expected in the projection years 2 - Nyrs;
#'     a vector obtained from the controlfile
#'
#' @return a list of numpop lists of 19 objects as detailed above.
#' @export
#'
#' @examples
#' print("See code for makezoneC to see usage of makeabpop")
makeabpop <- function(popparam,midpts,initLML) {
  #  popparam=popdef;midpts=midpts;initLML=initialLML
  #(DLMax,L50,L95,SigMax,SaMa,SaMb,Wta,Wtb,Me,L50C,L95C,R0 SelP[1],SelP[2],Nyrs,steeph,MaxCE,L50Mat block
  #   1    2   3    4      5    6   7   8   9 10   11   12  13     14      15   16     17    18     19
  numYr <- popparam["Nyrs"]
  N <- length(midpts)
  G <- STM(popparam[1:4],midpts)  # c(DLMax, L50, L95, SigMax)
  mature <- maturity(popparam["SaMa"],popparam["SaMb"],midpts)
  WtL <- WtatLen(popparam["Wta"],popparam["Wtb"],midpts)
  emergent <- logistic(popparam["L50C"],popparam["deltaC"],midpts)
  MatWt <- mature * WtL
  zSelect <- matrix(0,nrow=N,ncol=numYr)#,dimnames=list(midpts,1:numYr))
  blk <- popparam["block"]
  selL50 <- popparam["SelP1"]
  selL95 <- popparam["SelP2"]
  if (ncol(initLML) > 2) { # ensures correct sau LML selected if sauLML
    pickLML <- grep(popparam["SAU"],colnames(initLML))
  } else {
    pickLML <- 2
  }
  histLML <- initLML[,pickLML]
  verLML <- unique(histLML) # find unique LML in the vetor of LML
  for (LML in verLML) {
    Sel <- logistic((LML+selL50),selL95,midpts)
    Sel <- Sel * emergent # include emergence to determine availability
    pick <- which(histLML == LML)
    zSelect[,pick] <- rep(Sel,length(pick))
  }
  zSelWt <- zSelect * WtL
  zLML <- histLML
  Me <- as.numeric(popparam["Me"])
  AvRec <- as.numeric(popparam["AvRec"])  # R0
  qest <- as.numeric(popparam["qest"])
  lambda <- as.numeric(popparam["lambda"])
  scalece <- 0.0  # instead of popq
  B0 <- 0.0
  ExB0 <- 0.0
  SaM <- as.numeric(-popparam["SaMa"]/popparam["SaMb"]) # -a/b
  MSY <- 0
  MSYDepl <- 0
  bLML <- 0
  SAU <- 0
  ans <- list(Me,AvRec,B0,ExB0,MSY,MSYDepl,bLML,scalece,qest,lambda,SaM,
              popparam,zLML,G,mature,WtL,emergent,zSelect,zSelWt,MatWt,SAU)
  names(ans) <- c("Me","R0","B0","ExB0","MSY","MSYDepl","bLML","scalece","qest",
                  "lambda","SaM","popdef","LML","G","Maturity","WtL","Emergent",
                  "Select","SelWt","MatWt","SAU")
  class(ans) <- "abpop"
  return(ans)
} # End of makeabpop

#' @title makeequilzone high level function generates equilibrium zone
#'
#' @description makeequilzone is a high level function that merely hides the
#'     details of generating the original unfished zone after reading in the
#'     data files, estimates the productivity, and sets up the results
#'     directory, rundir, ready to receive files.
#'
#' @param rundir the directory containing the control csv files. It
#'     can/will also act to store results in a manner that will allow them
#'     to be displayed using makehtml.
#' @param ctrlfile the main file that controls the particular run. It contains
#'     the name of the data file that is used to biologically condition the
#'     numpop populations
#' @param doproduct boolean, should the productivity calculations be made
#'     during the conditioning. Set to FALSE conditionOM
#' @param uplimH defines the upper limit of harvest used when estimating the
#'     productivity (also important when initial depletion is not 1.0). The
#'     default = 0.4
#' @param incH defines the interval between H steps when estimating productivity
#'     default = 0.005
#' @param verbose Should progress comments be printed to console, default=TRUE
#'
#' @return a list of zoneC, zoneD, glb, constants, saudat,product, ctrl, and zone1
#' @export
#'
#' @examples
#' print("wait on datafiles")
makeequilzone <- function(rundir,ctrlfile="control.csv",doproduct=TRUE,
                          uplimH=0.4,incH=0.005,verbose=TRUE) {
 #  rundir=rundir;ctrlfile=controlfile;doproduct=FALSE; verbose=TRUE;uplimH=0.35;incH=0.005
  zone1 <- readctrlfile(rundir,infile=ctrlfile,verbose=verbose)
  ctrl <- zone1$ctrl
  glb <- zone1$globals     # glb without the movement matrix
  bysau <- zone1$ctrl$bysau
  opar <- NULL
  parsin <- zone1$condC$parsin
  if (parsin) opar <- as.matrix(zone1$condC$optpars)
  if (is.null(bysau)) bysau <- 0
  if (bysau) {
    saudata <- readsaudatafile(rundir,ctrl$datafile,optpar=opar,verbose=verbose)
    constants <- saudata$constants
    saudat <- saudata$saudat
    zone1$condC$poprec <- saudata$poprec
  } else {
    constants <- readpopdatafile(rundir,ctrl$datafile)
    saudat <- constants
  }
  if (verbose) cat("Files read, now making zone \n")
  out <- setupzone(constants,zone1,doproduct,uplim=uplimH,inc=incH,
                   verbose=verbose) # make operating model
  zoneC <- out$zoneC
  zoneD <- out$zoneD
  glb <- out$glb             # glb now has the movement matrix
  product <- out$product     # important bits usually saved in rundir
  zone1$globals <- glb
  if (parsin) {
    rewritecontrolfile(rundir,zone1,ctrlfile)
    rewritedatafile(rundir,zone1,saudat)
    if (verbose)
      cat("Control and Data files rewritten after amending values from sizemod \n")
  }
  # did the larval dispersal level disturb the equilibrium?
  # zoneD <- testequil(zoneC,zoneD,glb,verbose=verbose)
  # ans <- resetexB0(zoneC,zoneD) # rescale exploitB to avexplB after dynamics
  # zoneC <- ans$zoneC
  # zoneD <- ans$zoneD
  equilzone <- list(zoneC=zoneC,zoneD=zoneD,glb=glb,constants=constants,
                    saudat=saudat,product=product,ctrl=ctrl,zone1=zone1)
  return(equilzone)
} # end of makeequilzone

#' @title makemove produces the movement matrix
#'
#' @description makemove produces the movement matrix, which assumes
#'     that only adjacent populations contribute to each other. Thus,
#'     a central population would lose half its larvae to one side and
#'     half to the other side, but would receive half of their larval
#'     dispersal each. As larval dispersal is modelled as a proportion
#'     those populations that produce more larvae will lose more in
#'     absolute terms. Those populations at the edges of the zone
#'     only lose half the larval dispersal into the zone and it is
#'     assumed that what they lose out of the zone will be matched
#'     by what will come into the zone.
#'
#' @param npop the number of populations
#' @param recs the original unfished recruitment levels estimated for
#'     each population
#' @param ld the larval dispersal rate for all populations
#' @param sigmove currently no variation is assumed to occur. This is
#'     here in case it needs to get implemented, but as there is
#'     little or no data on actual larval dispersal levels currently
#'     this is set to zero and not used.
#'
#' @return a npop x npop matrix describing movement among populations
#' @export
#'
#' @examples
#' recruits <- c(636146,263878,819189,1112513,285025,671573.9)
#' makemove(npop=6,recs=recruits,ld=0.04)
makemove <- function(npop,recs,ld,sigmove=0.0) {
  ldh <- ld/2.0
  if (npop < 3) {
    if (npop == 1) {
      move <- matrix(0,nrow=1,ncol=1,dimnames=list(1:1,1:1))
      move[1,1] <- 1
    }
    if (npop == 2) {
      move <- matrix(c((1-ldh),ldh,ldh,(1-ldh)),nrow=npop,ncol=npop,
                     dimnames=list(1:npop,1:npop))
    }
  } else {
    move <- matrix(0,nrow=npop,ncol=npop,dimnames=list(1:npop,1:npop))
    move[1,1:2] <- c((1-ldh),ldh)
    recseq <- c(ldh,(1-ld),ldh)
    for (r in 2:(npop-1)) {
      move[r,] <- c(rep(0,(r-2)),recseq,rep(0,(npop-3-(r-2))))
    }
    move[npop,(npop-1):npop] <- c(ldh,(1-ldh))
  }
  return(move)
} # end of makemove

#' @title makezoneC makes the constant part of the simulation
#'
#' @description makezoneC makes the constant part of the simulated
#'     zone. Once defined this does not change throughout the
#'     simulation. Once made it still requires makezone to be run
#'     to fill in the B0, ExBo, MSY, MSYDepl, and the scalece values, and
#'     to produce zoneD, the dynamic part of the new zone
#'
#' @param zone the object derived from the readzonefile function
#' @param const the object derived from the readdatafile function
#'
#' @return a list containing the constant part of the simulated zone
#' @export
#'
#' @examples
#' \dontrun{
#' data(zone1)
#' data(constants)
#' ans <- makezoneC(zone=zone1,const=constants)
#' zoneC <- ans$zoneC
#' popdefs <- ans$popdefs
#' str(zoneC,max.level=1)
#' str(zoneC[[1]])  # not complete at this stage
#' print(popdefs)
#' }
makezoneC <- function(zone,const) { # zone=zone1; const=constants
  glb <- zone$globals
  nSAU <- glb$nSAU
  numpop <- glb$numpop
  midpts <- glb$midpts
  blkdef <- zone$SAUpop
  SAUindex <- defineBlock(nSAU,blkdef,numpop)
  glb$sauindex <- SAUindex
  saunames <- glb$saunames
  SAU <- saunames[SAUindex]     #as.numeric(const["SAU",])
  if (zone$randomseed > 0) {
      set.seed(zone$randomseed)
    } else {
      set.seed()
  }
  initialLML <- cbind(year=glb$hyrnames,histLML=rep(zone$initLML,glb$hyrs))
  rownames(initialLML) <- glb$hyrnames
  if (zone$catches > 0) initialLML <- zone$condC$histyr
  if (zone$ctrl$bysau) {
    popdefs <- definepops(nSAU,SAUindex,const,glob=glb) # define pops
  } else {
    popdefs <- popsdefine(const,glob=glb)
  }
  zoneC <- vector("list",numpop)
  for (pop in 1:numpop) {      # pop <- 1
    popdef <- popdefs[pop,]
    zoneC[[pop]] <- makeabpop(popdef,midpts,initialLML)
    tmpL <- oneyrgrowth(zoneC[[pop]],zoneC[[pop]]$SaM) #SaM defined in makeabpop
    zoneC[[pop]]$bLML <- oneyrgrowth(zoneC[[pop]],tmpL) #SaM + 2 year's growth
    zoneC[[pop]]$SAU <- SAU[pop]
  }
  recs <- getvar(zoneC,"R0") #sapply(zoneC,"[[","R0")
  larvdisp <- glb$larvdisp
  move <- makemove(numpop,recs,larvdisp)
  glb$move <- as.matrix(move)
  glb$SAUnum <- as.vector(const["SAU",])
  class(zoneC) <- "zoneC"
  ans <- list(zoneC=zoneC,glb=glb,popdefs=popdefs)
  return(ans)
}  # End of makezoneC

#' @title makezone generates the dynamic parts of the simulated zone
#'
#' @description makezone generates the dynamics components of the
#'     simulated zone and completes the constant components of the
#'     simulated zone. The term 'zone' refers to the upper level
#'     of geographical detail that is used. Thus, for example, in
#'     Tasmania we might simulate a number of statistical blocks with
#'     multiple populations. In combination, we refer to the total as
#'     a zone. The matrices are Nyrs x numpop and the two
#'     arrays are N x Nyrs x numpop.
#'
#' @param glob the global constants defined for the current simulation.
#'     These include numpop, nblock, midptsd, Nclass, and Nyrs
#' @param zoneC the zoneC object from makezoneC
#'
#' @return a list of the dynamics and the constant components of the
#'     simulation
#' @export
#'
#' @examples
#' \dontrun{
#'   data(zone1)
#'   glb <- zone1$globals
#'   data(constants)
#'   ans <- makezoneC(zone1,constants)
#'   zoneC <- ans$zoneC
#'   glb <- ans$glb
#'   ans2 <- makezone(glb,zoneC)
#'   str(ans2,max.level=2)
#'  }
makezone <- function(glob,zoneC) { # glob=glb; zoneC=zoneC;
  hyrs <- glob$hyrs
  numpop <- glob$numpop
  N <- glob$Nclass
  ExplB <- matrix(0,nrow=hyrs,ncol=numpop)
  midyexpB <- matrix(0,nrow=hyrs,ncol=numpop)
  MatB <- matrix(0,nrow=hyrs,ncol=numpop)
  Catch <- matrix(0,nrow=hyrs,ncol=numpop)
  Harvest <- matrix(0,nrow=hyrs,ncol=numpop)
  cpue <- matrix(0,nrow=hyrs,ncol=numpop)
  deplExB <- matrix(0,nrow=hyrs,ncol=numpop)
  deplSpB <- matrix(0,nrow=hyrs,ncol=numpop)
  Recruit <- matrix(0,nrow=hyrs,ncol=numpop)
  CatchN <- array(data=0,dim=c(N,hyrs,numpop))
  Nt <- array(data=0,dim=c(N,hyrs,numpop))
  NumNe <- Nt
  SAU <- getvar(zoneC,"SAU") #as.numeric(sapply(zoneC,"[[","SAU"))
  recs <- getvar(zoneC,"R0") #sapply(zoneC,"[[","R0")
  move <- glob$move
  newrecs <- move %*% recs
  for (pop in 1:numpop) {  # pop=1
    SurvE <- exp(-zoneC[[pop]]$Me)
    hSurv <- exp(-zoneC[[pop]]$Me/2.0)
    recr <- rep(0,N)
    recr[1] <- newrecs[pop]  # change this once I have imported move
    UnitM <- matrix(0,nrow=N,ncol=N)
    diag(UnitM) <- 1.0
    G <- zoneC[[pop]]$G
    Minv <- solve(UnitM - (SurvE * G))
    Nt[,1,pop] <- Minv %*% recr # initial unfished numbers-at-size
    MatB[1,pop] <- sum(zoneC[[pop]]$MatWt*Nt[,1,pop])/1e06
    zoneC[[pop]]$B0 <- MatB[1,pop] # mature biomass at start of year
    newNt1 <- (hSurv * (G %*% Nt[,1,pop]))
    newNt2 <- (hSurv * newNt1)
    popSel <- zoneC[[pop]]$SelWt[,1]
    preexB <- sum(popSel*newNt1)/1e06
    postexB <- sum(popSel*newNt2)/1e06
    ExplB[1,pop] <- (preexB + postexB)/2 # mean exploitB before/after halfyear
    zoneC[[pop]]$ExB0 <- ExplB[1,pop]
    deplExB[1,pop] <- 1.0  # no depletion when first generating zones
    deplSpB[1,pop] <- 1.0
    Recruit[1,pop] <- recr[1]
    #qcalc <- as.numeric(zoneC[[pop]]$popdef["MaxCE"])
    zoneC[[pop]]$scalece<- 0 # used to scale all pops to maxce
  }
  nsau <- glob$nSAU
  sauindex <- glob$sauindex
  lambda <- zoneC[[1]]$lambda  # constant across the zone
  popexb0 <- getlistvar(zoneC,"ExB0")     # pop level expB
  sauexb0 <- sumpop2sau(popexb0,sauindex) # sau level expB
  for (sau in 1:nsau) { # sau=1
    pickpop <- which(sauindex == sau)
    npop <- length(pickpop)
    exb0 <- sauexb0[sau]
    maxcpue <- zoneC[[pickpop[1]]]$qest * (exb0 ^ lambda)
    for (popi in 1:npop) { # popi = 1 # populations in each sau
      zoneC[[pickpop[popi]]]$popdef["MaxCE"] <- maxcpue
      scalece <- exb0/popexb0[pickpop[popi]]
      zoneC[[pickpop[popi]]]$popdef["scalece"] <- scalece
      zoneC[[pickpop[popi]]]$scalece <- scalece
    }
  }
  for (pop in 1:numpop) {  # pop=1
    cpue[1,pop] <- zoneC[[pop]]$qest *
                   ((zoneC[[pop]]$scalece * ExplB[1,pop]) ^ lambda)
  }
  ans <- list(SAU=SAU,matureB=MatB,exploitB=ExplB,midyexpB=midyexpB,
              catch=Catch,harvestR=Harvest,cpue=cpue,recruit=Recruit,
              deplsB=deplSpB,depleB=deplExB,catchN=CatchN,Nt=Nt,NumNe=NumNe)
  return(list(zoneD=ans,zoneC=zoneC))
} # end of makezone

#' @title maturity Logistic maturity curve
#'
#' @description maturity this uses the logistic function:
#'   exp(a+bL)/(1+exp(a+bL)), which has the property that the SM50 = -a/b
#'   and the interquartile distance is 2.Ln(3)/b.
#' @param ina is the intercept of the exponential function
#' @param inb is the gradient of the exponential function
#' @param lens a vector of lengths for which the logistic maturity value
#'     will be calculated
#' @return A vector of length(lens) containing the predicted maturity at
#'     length values
#' @export
#'
#' @examples
#' \dontrun{
#' a <- -14.383
#' b <- 0.146017
#' lens <- seq(2,210,2)
#' Maturity <- maturity(a,b,lens)
#' }
maturity <- function(ina,inb,lens) {
   ans <- exp(ina+inb*lens)/(1+exp(ina+inb*lens))
   return(ans)
} # end of maturity



#' @title modzoneC runs the zone 3 x hyrs to equilibrium
#'
#' @description modzoneC runs the doproduction function so it can then
#'     appropriately fill in each population's MSY and MSYDepl values.
#'     It outputs both the finalized zoneC and the prodution array.
#'     Hence we need to push copies of both zoneC and zoneD into
#'     modzoneC
#'
#' @param zoneC the constant components of the simulated zone
#' @param zoneD the dynamic components of the simulated zone
#' @param glob the general global variables
#' @param lowlim the lower limit of harvest rate applied, default=0.0
#' @param uplim the upper limit of harvest rate applied, default=0.4
#' @param inc the harvest rate increment at each step, default=0.005
#'
#' @return a list containing the updated zoneC and product
#' @export
#'
#' @examples
#' \dontrun{  # the doproduction part takes too long to run
#'   data(zone1)
#'   glb <- zone1$globals
#'   data(constants)
#'   ans <- makezoneC(zone1,constants)
#'   zoneC <- ans$zoneC
#'   ans <- makezone(glb,zoneC)
#'   zoneC <- ans$zoneC
#'   zoneD <- ans$zoneD
#'   ans2 <- modzoneC(zoneC,zoneD,glb)
#'   str(ans2,max.level=2)
#' }   # zoneC=zoneC; zoneD=zoneD; glob=glb; lowlim=0.0;uplim=0.4;inc=0.01
modzoneC <- function(zoneC,zoneD,glob,lowlim=0.0,uplim=0.4,inc=0.005) {
  numpop <- glob$numpop
  production <- doproduction(zoneC=zoneC,zoneD=zoneD,glob=glob,
                             lowlim=lowlim,uplim=uplim,inc=inc)
  xval <- findmsy(production)
  for (pop in 1:numpop) {
    zoneC[[pop]]$MSY <- xval[pop,"Catch"]
    zoneC[[pop]]$MSYDepl <- xval[pop,"Deplet"]
  }
  return(list(zoneC=zoneC,product=production))
} # end of modzoneC

#' @title popsdefine makes a parameter vector, popdef, for each population
#'
#' @description popsdefine makes a parameter vector, popdef, for each
#'     population. However, instead of imputing random variation the input
#'     matrix of constants is already defined by population. This contrasts with
#'     the function definepops which is used if saudata is read instead of
#'     population level data.
#'
#' @param const a matrix with rows containing the constants defined for
#'     each of the biological parameters. The columns of the matrix
#'     reflect values for each populations from readpopdatafile.
#' @param glob the globals object from readctrlfile
#'
#' @seealso{
#'   \link{definepops}, \link{readpopdatafile}, \link{readsaudatafile}
#' }
#'
#' @return a matrix with a row for each population and whose columns
#'     are parameters defining the biological parameters for each population.
#' @export
#'
#' @examples
#' \dontrun{
#'   data(zone1)
#'   glb <- zone1$globals
#'   # this needs fixing
#' }
popsdefine <- function(const,glob) {
  #  inSAU=nSAU; inSAUindex=SAUindex; const=const; glob=glb
  numpop <- glob$numpop
  hyrs <- glob$hyrs
  columns <- c("DLMax","L50","L95","SigMax","SaMa","SaMb","Wta","Wtb","Me",
               "L50C","deltaC","AvRec","SelP1","SelP2","Nyrs","steeph",
               "MaxCE","L50mat","SAU","lambda","qest","scalece")
  popdefs <- matrix(0,nrow=numpop,ncol=length(columns),
                    dimnames=list(1:numpop,columns))
  popdefs[,"Nyrs"] <- rep(hyrs,numpop) # Num Years - why is this here?
  popdefs[,"Me"] <- const["Me",]
  popdefs[,"DLMax"] <- const["DLMax",]
  popdefs[,"L50"] <- const["L50",]
  popdefs[,"L95"] <- const["L95",]
  popdefs[,"SigMax"] <- const["SigMax",]
  popdefs[,"L50mat"] <- const["L50Mat",]
  popdefs[,"SaMa"] <- const["SaMa",]
  popdefs[,"SaMb"] <- const["SaMb",]
  popdefs[,"SelP1"] <- const["selP1",]
  popdefs[,"SelP2"] <- const["selP2",]
  popdefs[,"L50C"] <- const["L50C",]
  popdefs[,"deltaC"] <- const["deltaC",]
  popdefs[,"Wtb"] <- const["Wtb",]
  popdefs[,"Wta"] <- const["Wta",]
  popdefs[,"steeph"] <- const["steeph",]
  popdefs[,"AvRec"] <- const["AvRec",]
  popdefs[,"MaxCE"] <- 0 # gets value in makezone
  popdefs[,"SAU"] <- const["SAU",]
  popdefs[,"lambda"] <- const["lambda",]
  popdefs[,"qest"] <- const["qest",]
  popdefs[,"scalece"] <- 0 # gets value in makezone
  test <- popdefs[,"L95"] - popdefs[,"L50"]
  if (any(test < 0)) stop("L95 < L50 - adjust the input data file  \n")
  return(popdefs)
} # end of popsdefinepops


#' @title print.zoneDefinition S3 method for printing zonedef summary
#'
#' @description print.zoneDefinition S3 method for printing a summary
#'     of the zoneDef.
#' @param x a zonedefinition object
#' @param ... in case there are extra parameters
#' @return nothing, but does print out a summary of the zeonDef
#'
#' @export
#'
#' @examples
#' \dontrun{
#' txt1 <- 'all example code should be able to run'
#' }
print.zoneDefinition <- function(x, ...) {
  cat("\n Zone_Definition  \n")
  cat("\n Date made : ",x$date)
  cat("\n numpop    : ",x$numpop)
  cat("\n numBlock  : ",x$nBlock)
  cat("\n MSY       : ",round(x$summaryZone["MSY"],3))
  cat("\n B0        : ",round(x$summaryZone["B0"],3))
  cat("\n Av bLML   : ",round(x$summaryZone["AvbLML"],1))
  cat("\n Av SaM    : ",round(x$summaryZone["AvSaM"],1))
  cat("\n %prot SpB : ",round(x$summaryZone["%ProtB"],2))
  cat("\n \n ")
  cat("\n  1: $production  - productivity of each population")
  cat("\n  2: $defpop      - each population's definition - see defpops")
  cat("\n  3: $struct      - zones hierarchical structure")
  cat("\n  4: $nBlock      - number of blocks")
  cat("\n  5: $numpop      - number of populations")
  cat("\n  6: $summaryPop  - biological properties of populations")
  cat("\n  7: $summaryZone - summary zone properties")
  cat("\n  8: $blockProp   - summary block properties")
  cat("\n  9: $date        - date/time of zone creation")
  cat("\n \n ")
  NextMethod("print")
}

#' @title print.zone S3 method for printing a summary of a given zone
#'
#' @description print.zone S3 method for printing a summary of a given zone
#'
#' @param x a zone object
#' @param ... in case there are extra parameters
#' @return nothing, but does print out a summary of the zone
#'
#' @export
#'
#' @examples
#' \dontrun{
#' txt1 <- 'all example code should be able to run'
#' }
print.zone <- function(x, ...) {
  cat("Number of populations: ",length(x),"\n")
  cat("Number of years      : ",length(x[[1]]$ExploitB),"\n")
  cat("B0                   : ",round(sum(getlistvar(x,"B0")),3), "\n")
  cat("MSY                  : ",round(sum(getlistvar(x,"MSY")),3), "\n")
  cat("Size of zone (bytes) : ",object.size(x),"\n\n")
  NextMethod("print")
}


#' @title resetLML sets the LML for a given zone from a given year
#'
#' @description  resetLML sets the LML for a given zone from a given year
#'     out to the final year. It changes the selectivity, the SelWt values,
#'     and the LML value per block.
#' @param inzone a zone object
#' @param inLML a vector of LML for each block to be imposed from inyear on.
#' @param inyear the start year for the new LML
#' @param glob the global variables object
#' @return A zone object with the LML changed from inyear on
#' @export
#'
#' @examples
#' \dontrun{
#' txt1 <- 'all example code should be able to run'
#' }
resetLML <- function(inzone,inLML,inyear,glob) {
  Npop <- glob$numpop
  hyrs <- glob$hyrs
  midpts <- glob$midpts
  if (length(inLML) != glob$nblock)
    stop("the number of LML does not match the number of blocks")
  for (pop in 1:Npop) {  # inzone <- zone; pop <- 1
    popdef <- inzone[[pop]]$popdef
    blk <- popdef["block"]
    Sel <- logistic((inLML[blk]+popdef["SelP1"]),
                    (inLML[blk]+popdef["SelP2"]),midpts)
    inzone[[pop]]$Select[,inyear:hyrs] <- Sel
    inzone[[pop]]$SelWt[,inyear:hyrs] <- Sel * inzone[[pop]]$WtL
    inzone[[pop]]$LML[inyear:hyrs] <- inLML[blk]
  }
  return(inzone)
}  # end of resetLML

#' @title resetexB0 resets the unfished exploitable biomass at time zero
#'
#' @description resetexB0 sets the unfished exploitable biomass at time
#'     zero to the average of the exploitable biomass levels before and
#'     after any fishing (when there are zero catches). Hence growth and
#'     natural mortality are included.
#'
#' @param zoneC the constant components of the simulated zone
#' @param zoneD the dynamic components of the simulated zone
#'
#' @return a refreshed zoneC with updated ExB0 values and zoneD$depleB=1
#' @export
#'
#' @examples
#' print("wait on data")
resetexB0 <- function(zoneC,zoneD) {
  numpop <- length(zoneC)
  for (pop in 1:numpop) {
    zoneC[[pop]]$ExB0 <- zoneD$exploitB[1,pop]
    zoneD$depleB[1,pop] <- 1.0
  }
  return(list(zoneC=zoneC,zoneD=zoneD))
} # end of resetexB0

# nsau <- glob$nSAU
# sauindex <- glob$sauindex
# lambda <- zoneC[[1]]$lambda  # constant across the zone
# popexb0 <- getlistvar(zoneC,"ExB0")     # pop level expB
# sauexb0 <- sumpop2sau(popexb0,sauindex) # sau level expB
# for (sau in 1:nsau) { # sau=1
#   pickpop <- which(sauindex == sau)
#   npop <- length(pickpop)
#   exb0 <- sauexb0[sau]
#   maxcpue <- zoneC[[pickpop[1]]]$qest * (exb0 ^ lambda)
#   for (popi in 1:npop) { # popi = 1 # populations in each sau
#     zoneC[[pickpop[popi]]]$popdef["MaxCE"] <- maxcpue
#     zoneC[[pickpop[popi]]]$qest <-
#       as.numeric(maxcpue/(popexb0[pickpop[popi]] ^ lambda))
#   }
# }
# for (pop in 1:numpop) {  # pop=1
#   cpue[1,pop] <- zoneC[[pop]]$qest * (ExplB[1,pop] ^ lambda)
# }

#' @title setupzone makes zone's constant, dynamic, and productivity parts
#'
#' @description setupzone makes the zone's constant, dynamic, and
#'     productivity parts returning them, along with glb, in a list. The
#'     objective of this function is to generate the orignal unfished
#'     equilibrium zone of nSAU SAU, and numpop populations
#'
#' @param constants the population constants derived from readdatafile
#' @param zone1 the zonal object driving the construction
#' @param doproduct boolean, should the productivity calculations be made
#'     during the conditioning. defined in do_MSE and makeequilzone
#' @param uplim the upper limit of harvest rate applied, default=0.4
#' @param inc the harvest rate increment at each step, default=0.005
#' @param verbose Should progress comments be printed to console, default=TRUE
#'
#' @return a list of zoneC, zoneD, product, and glb the main
#'     components of the zone
#' @export
#'
#' @examples
#' \dontrun{
#' data(constants)
#' data(zone1)
#' out <- setupzone(constants,zone1)
#' zoneC <- out$zoneC
#' glb <- out$glb
#' str(zoneC[[1]])
#' str(glb)
#' }
setupzone <- function(constants,zone1,doproduct,uplim=0.4,inc=0.005,verbose=TRUE) {
  # constants=constants; zone1=zone1; doproduct=TRUE; uplim=0.4; inc=0.001; verbose=TRUE
  ans <- makezoneC(zone1,constants) # initiates zoneC
  zoneC <- ans$zoneC
  glb <- ans$glb
  ans <- makezone(glob=glb,zoneC=zoneC) # make zoneD, add cpue, qest to zoneC
  zoneC <- ans$zoneC  # zone constants
  zoneD <- ans$zoneD  # zone dynamics
  product <- NULL
  if (doproduct) {
    if (verbose)
      cat("Now estimating population productivity between H ",inc," and ",
          uplim, "\n")
    # adds productivity, and MSY, MSYdepl to zoneC if doproduct=TRUE
    ans <- modzoneC(zoneC=zoneC,zoneD=zoneD,glob=glb,uplim=uplim,inc=inc)
    zoneC <- ans$zoneC  # zone constants
    product <- ans$product  # productivity by population
  }
  out <- list(zoneC=zoneC, zoneD=zoneD, product=product,glb=glb)
  return(out)
} # end of setupzone


#' @title STM Generates the Size Transition Matrix for Inverse Logistic
#'
#' @description STM With the input of the four parameters inside a vector,
#'     and a vector of initial lengths or mid-points of size classes STM
#'     generates a square transition matrix with the probabilities of
#'     growing from each initial size into the same or larger sizes.
#'     Negative growth is disallowed. All columns in the matrix sum to one.
#' @param p a vector of four parameters in the following order
#'     MaxDL the maximum growth increment of the inverse logistic,
#'     L50 the initial length at which the inflexion of the growth
#'     increment curve occurs, L95 - the initial length that defines the
#'     95th percentile of the growth increments, SigMax - the maximum
#'     standard deviation of the normal distribution used to describe the
#'     spread of each distribution of growth increments
#' @param mids a vector of initial lengths which also define the width of
#'     each size class thus, mids from 2 - 210 would imply 2mm size
#'     classes 1 - 3 = 2, 3 - 5 = 4, etc
#' @return A square matrix with dimension =the length of the mids vector
#' @export
#'
#' @references Haddon, M., Mundy, C., and D. Tarbath (2008) Using an
#'     inverse-logistic model to describe growth increments of blackip
#'     abalone (Haliotis rubra) in Tasmania. Fisheries Bulletin 106: 58-71
#' @examples
#'  param <- c(25.0,120.0,170.0,4.0)
#'  midpts <- seq(2,210,2)
#'  G <- STM(param,midpts)
#'  print(round(G[1:30,1:8],4))
STM <- function(p,mids) { #    # p <- popparam[1:4]; mids <- midpts
   n <- length(mids)
   G <- matrix(0,nrow=n,ncol=n)
   cw <- mids[2]-mids[1]
   SigL  <- p[4]/((1+exp(log(19.0)*(mids-p[3])/(mids[n]-p[3]))))
   MeanL <- mids + (p[1]/((1+exp(log(19.0)*(mids-p[2])/(p[3]-p[2])))))
   for (j in 1:n) {
      for (i in 1:n) {
         Prob <- (1-pnorm(mids[i]+cw/2.0,MeanL[j],SigL[j],FALSE))
         if (i < j)  { G[i,j] <- 0.0 }
         if (i == j) { G[i,j] <- Prob }
         if (i > j)  { G[i,j] <- Prob - (1-pnorm(mids[i-1]+cw/2.0,MeanL[j],
                                                 SigL[j],FALSE)) }
      }
   }
   G[n,] <- G[n,]+ (1-colSums(G)) # plus group
   class(G) <- "STM"
   return(G)
} # end of STM

#' @title testequil runs a zone for hyrs and determines stability
#'
#' @description testequil runs a given zone for hyrs at the given
#'     harvest rate, and then tests that the last values of matureB,
#'     exploitB, recruitment, and spawning biomass depletion are the
#'     same as the first (to three decimal places). It reports this
#'     to the console if verbose=TRUE. This is used with a harvest rate
#'     of zero and no variation in recruitment when defining the equilibrium
#'     zone under the application of the movement matrix.
#'
#' @param zoneC the constants components of the simulated zone
#' @param zoneD the dynamic components of the simulated zone
#' @param glb the global variables
#' @param inH a vector of numpop harvest rates
#' @param verbose should results go to the console, default=TRUE
#'
#' @return the dynamics component with hyrs of dynamics
#' @export
#'
#' @examples
#' \dontrun{  # modzoneC takes too long to run because of doproduction
#'  data(zone)
#'  zoneDe <- testequil(zoneC=zone$zoneC,zoneD=zone$zoneD,glb=zone$glb)
#' }    #zoneC=zoneC; zoneD=zoneD; glb=glb; inH=0.0; verbose=TRUE
testequil <- function(zoneC,zoneD,glb,inH=0.0,verbose=TRUE) {
  hyrs <- glb$hyrs
  Nclass <- glb$Nclass
  npop <- glb$numpop
  inHarv <- rep(inH,npop)
  for (yr in 2:hyrs)
    zoneD <- oneyearD(zoneC=zoneC,zoneD=zoneD,
                      inHt=inHarv,year=yr,sigmar=1e-08,
                      Ncl=Nclass,npop=npop,movem=glb$move)
  zoneD <- restart(oldzoneD=zoneD,hyrs=hyrs,npop=npop,N=Nclass,zero=FALSE)
  for (yr in 2:hyrs)  # repeat to be sure
    zoneD <- oneyearD(zoneC=zoneC,zoneD=zoneD,
                      inHt=inHarv,year=yr,sigmar=1e-08,
                      Ncl=Nclass,npop=npop,movem=glb$move)
  if (verbose) {
    if (all(trunc(zoneD$matureB[1,],2) == trunc(zoneD$matureB[hyrs,],2))) {
      cat("matureB Stable \n")
    } else {
      cat("matureB varies \n")
    }
    expldiff <- trunc(zoneD$exploitB[1,],2) == trunc(zoneD$exploitB[hyrs,],2)
    if (all(expldiff)) {
      cat("exploitB Stable \n")
    } else {
      diffe <- abs(trunc(zoneD$exploitB[1,],2) - trunc(zoneD$exploitB[hyrs,],2))
      label <- paste0("exploitB varies ",max(diffe,na.rm=TRUE)," \n")
      cat(label)
    }
    if (all(trunc(zoneD$recruit[1,],2) == trunc(zoneD$recruit[hyrs,],2))) {
      cat("recruitment Stable \n")
    } else {
      cat("recruitment varies \n")
    }
    if (all(round(zoneD$deplsB[1,],2) == round(zoneD$deplsB[hyrs,],2))) {
      cat("spawning depletion Stable \n")
    } else {
      cat("spawning depletion varies \n")
    }
  }
  zoneD <- restart(oldzoneD=zoneD,hyrs=hyrs,npop=npop,N=Nclass,zero=TRUE)
  return(zoneD)
} # end of testequil


#' @title WtatLen Power function to describe weight at length relationship
#'
#' @description WtatLen This can be used to generate any power function.
#' @param ina is the intercept of the power function
#' @param inb is the gradient (or the explonent) of the power function
#' @param lens a vector of lengths for which the weight will be calculated
#' @return A vector of length(lens) containing the predicted values
#'
#' @export
#'
#' @examples
#' \dontrun{
#' a <- 0.0000542856
#' b <- 3.161415
#' lens <- seq(2,210,2)
#' wtL <- WtatLen(a,b,lens)
#' }
WtatLen <- function(ina,inb,lens) {
   ans <- ina*lens^inb
   return(ans)
}  #end of WtatLen

