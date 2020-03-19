
#' @title defineBlock subdivides defined number of populations into blocks
#'
#' @description subdivides the defined number of populations into blocks
#'     in the numbers defined in the constant blkpop read into 'constants'.
#'     It allocates the populations sequentially to numblk blocks; if
#'     numblk does not divide into numpop exactly, any remainder is put
#'     into the final block. For block can read SMU = spatial management
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
#'     distributions defined in bConst=condDat$constants. The input can
#'     be either the whole of condDat or just condDat$constants, which
#'     contains the definitions of the probability density functions
#'     used to define each population.
#'
#' @param inSMU the number of SMUs defined in the region, a scalar
#' @param inSMUindex a vector of the block index for each population
#' @param bConst a matrix with rows containing the PDF parameters for
#'     each of the biological parameters. The columns of the matrix
#'     reflect values for each populations from readdatafile
#' @param glob the globals object from readregionfile
#'
#' @return a matrix with a row for each population and whose columns
#'     are parameters defining the biological parameters for each population.
#' @export
#'
#' @examples
#' print("need to wait on built in data sets")
definepops <- function(inSMU,inSMUindex,bConst,glob) {
  #  inSMU=nSMU; inSMUindex=SMUindex; bConst=const; glob=glb
  numpop <- glob$numpop
  Nyrs <- glob$Nyrs
  columns <- c("DLMax","L50","L95","SigMax","SaMa","SaMb","Wta","Wtb","Me",
               "L50C","deltaC","AvRec","SelP1","SelP2","Nyrs","steeph",
               "MaxCE","L50mat","SMU")
  popdefs <- matrix(0,nrow=numpop,ncol=length(columns),
                    dimnames=list(1:numpop,columns))
  popdefs[,"Nyrs"] <- rep(Nyrs,numpop) # Num Years - why is this here?
  for (pop in 1:numpop) {  # blk=2
    popdefs[pop,"Me"] <- rnorm(1,mean=bConst["Me",pop],
                               sd=bConst["sMe",pop])
    popdefs[pop,"DLMax"] <- rnorm(1,mean=bConst["DLMax",pop],
                                  sd=bConst["sMaxDL",pop])   # DLMax
    popdefs[pop,"L50"] <- rnorm(1,mean=bConst["L50",pop],
                                bConst["sL50",pop])
    popdefs[pop,"L95"] <- popdefs[pop,"L50"] +
      rnorm(1,mean=bConst["L50inc",pop],bConst["sL50inc",pop]) #L95
    popdefs[pop,"SigMax"] <- rnorm(1,mean=bConst["SigMax",pop],
                                   sd=bConst["sSigMax",pop])   # SigMax
    popdefs[pop,"L50mat"] <- rnorm(1,mean=bConst["L50Mat",pop],
                                   bConst["sL50Mat",pop])   # L50Mat
    popdefs[pop,"SaMa"] <- bConst["SaMa",pop]
    popdefs[pop,"SaMb"] <- -popdefs[pop,"SaMa"]/popdefs[pop,"L50mat"]
    popdefs[pop,"SelP1"] <- bConst["selL50p",pop]
    popdefs[pop,"SelP2"] <-  bConst["selL95p",pop]
    popdefs[pop,"L50C"] <- rnorm(1,mean=bConst["L50C",pop],
                                 bConst["sL50C",pop])
    popdefs[pop,"deltaC"] <-rnorm(1,mean=bConst["deltaC",pop],
                                  bConst["sdeltaC",pop])
    popdefs[pop,"Wtb"] <- rnorm(1,mean=bConst["Wtb",pop],
                                bConst["sWtb",pop])  # wt parameters
    popdefs[pop,"Wta"] <- bConst["Wtbtoa",pop] *
      popdefs[pop,"Wtb"]^bConst["sWtbtoa",pop]
    popdefs[pop,"steeph"] <- rnorm(1,mean=bConst["defsteep",pop],
                                   bConst["sdefsteep",pop])
    popdefs[pop,"AvRec"] <- rlnorm(1,meanlog=bConst["AvRec",pop],
                                   bConst["sAvRec",pop])
    popdefs[pop,"MaxCE"] <- rnorm(1,mean=bConst["MaxCEpars",pop],
                                  bConst["sMaxCEpars",pop])
    popdefs[pop,"SMU"] <- bConst["SMU",pop]
  }
  test <- popdefs[,"L95"] - popdefs[,"L50"]
  if (any(test < 0)) stop("L95 < L50 - adjust the input data file  \n")
  return(popdefs)
} # end of definepops


#' @title doproduction - estimates production curve for all populations
#'
#' @description doproduction estimates a surplus production curve for all
#'    populations. This estimates a zone's production by estimate MSY for
#'    each population separately and then combining. do.zoneProd is more
#'    accurate for a zone being fished as a whole. 'doproduction' is
#'    usually applied sequentially to all populations in a zone when using
#'    the makeZone function.
#'
#' @param inpop the population to be assessed; an 'abpop'
#' @param uplim the upper limit on teh harvest rate to test, defaults to
#'    0.4; a scalar
#' @return a list made up of the MSY, the depletion level at MSY, and a
#'    matrix with six columns: ExB, MatB, AnnH, Catch, Deplet, and RelCE.
#' @export
#'
#' @examples
#' \dontrun{
#' txt1 <- 'generate a zone and apply do.zoneProd and print the results'
#' str(txt1)
#' }
doproduction <- function(inpop,uplim=0.4) { # inpop=zone[[1]]; uplim=0.35
   imph <- seq(0.01,uplim,0.01)
   numrow <- length(imph)              # hyr <- 1; year <- 2
   Nyrs <- length(inpop$MatureB)
   columns <- c("ExB","MatB","AnnH","Catch","Deplet","RelCE")
   results <- matrix(0,nrow=numrow,ncol=6,dimnames=list(imph,columns))
   for (hyr in 1:numrow) {    # do the dynamics for numrow diff H values
      for (year in 2:Nyrs) { # always leaving yr1 the same hyr=1;year=2
         ExpB <- inpop$ExploitB[year-1]
         catch <- imph[hyr] * ExpB
         out <- oneyear(inpop,catch,year)
         inpop$ExploitB[year] <- out$ExploitB
         inpop$MatureB[year] <- out$MatureB
         inpop$Catch[year] <- out$Catch
         inpop$HarvestR[year] <- out$Harvest
         inpop$Nt[,year] <- out$Nt
      }

      results[hyr,"ExB"] <- inpop$ExploitB[Nyrs]
      results[hyr,"MatB"] <- inpop$MatureB[Nyrs]
      results[hyr,"AnnH"] <- inpop$HarvestR[Nyrs]
      results[hyr,"Catch"] <- inpop$Catch[Nyrs]
      results[hyr,"Deplet"] <- results[hyr,"MatB"]/inpop$B0
      results[hyr,"RelCE"] <- inpop$popq * results[hyr,"ExB"]
   }
   msy <- max(results[,"Catch"])
   pick <- which(results[,"Catch"] == msy)
   msyDepl <- results[pick,"Deplet"]
   ans <- list(msy,msyDepl,results)
   names(ans) <- c("MSY","MSYDepl","Productivity")
   return(ans)
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
#' @param inzone the zone to be characterized.
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
fillzoneDef <- function(inzone,prod) {  # inzone=zone; prod=production
   defNames <- c("production","defpop","struct","nBlock","numpop",
                 "summaryPop","summaryZone","blockProp","date")
   zDef <- vector("list",length(defNames))
   names(zDef) <- defNames
   numpop <- length(inzone)
   pops <- seq(1,numpop,1)
   popdefs <- t(sapply(inzone,"[[","popdef"))
   blockI <- popdefs[,"block"]
   nblock <- length(unique(blockI))
   struct <- matrix(0,nrow=numpop,ncol=3,
                    dimnames=list(pops,c("Blocks","Population","Hexagon")))
   struct[,1] <- blockI
   struct[,2] <- pops
   ans <- zoneProperty(inzone)
   #ndef <- length(zone[[1]]$popdef)
   #defname <- names(zone[[1]]$popdef)

   zDef[["production"]] <- production
   zDef[["defpop"]] <- popdefs  # already defined
   zDef[["struct"]] <- struct
   zDef[["nBlock"]] <- nblock
   zDef[["numpop"]] <- numpop
   zDef[["summaryPop"]] <- ans$summaryMatrix
   zDef[["summaryZone"]] <- ans$ZoneSummary
   zDef[["blockProp"]] <- summaryBlock(inzone)
   zDef[["date"]] <- date()
   class(zDef) <- "zoneDefinition"
   return(zDef)
}  # end of fillzoneDef

#' @title logistic a Logistic selectivity function
#'
#' @description logistic a logistic selectivity function. This uses the
#'     logistic function 1/(1+exp(-log(19.0)*(lens-inL50)/(delta))),
#'     where delta = inL95 - inL50. This explicitly defines the SM50
#'     but uses delta (which is SM95-SM50) as the second. This ensures
#'     that when adding variation to parameters, to vary between
#'     populations, when SM95 and SM50 are close together it is not
#'     possible for SM50 to become larger than SM95.
#'
#' @param inL50 is the length at 50 percent selection
#' @param delta is the difference between the 95 percent selection and
#'     the 50 percent selection
#' @param lens a vector of lengths for which the logistic maturity value
#'     will be calculated
#' @param knifeedge defaults to 0. If knifeedge is set to a particular
#'     length then the selectivity less than the value of knifeedge is set
#'     to zero.
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
#' }
logistic <- function(inL50,delta,lens,knifeedge=0) {
   ans <- 1/(1+exp(-log(19.0)*(lens-inL50)/(delta)))
   if (knifeedge > 0) {
      pick <- which(lens < knifeedge)
      if (length(pick) > 0) ans[1:pick[1]] <- 0.0
   }
   return(ans)
} # end of logistic

#' @title makeabpop generates full population structure in an unfished state
#'
#' @description makeabpop generates the full population structure in an
#'    unfished state. The structure is pre-determined: MaxDL L50 L95 MaxSig
#'    Me Mc R0 A0 B0 steeph ExploitB MatureB MatBCypt HarvestR Catch popdef
#'    MSY MSYDepl LML bLML popq SaM cpue CatchN deplExB Recruit ExB0 G
#'    Maturity WtL Emergent Select Nemerg SelWt MatWt. See the AbMSE
#'    documentation to see the full definition of an abpop and the structure
#'    of a zone
#'
#' @param popparam the vector of biological parameters values that define
#'     the specific properties of this population. Obtained from 'popdefs'
#' @param midpts the center values of the size classes from
#'     condDat{parseFile}
#' @param projLML the LML expected in the projection years 2 - Nyrs;
#'     matrix read in from datafile
#'
#' @return a list of objects as shown in the description.
#' @export
#'
#' @examples
#' \dontrun{
#'  data(condDat)
#'  bConst <- condDat$constants
#'  blkpop <- condDat$blkpop
#'  numpop <- condDat$numpop
#'  blockI <- defineBlock(nblock,blkpop,numpop)
#'  popdefs <- definepops(nblock,blockI,bConst)
#'  pop1 <- makeabpop(popdefs[1,],condDat$midpts,
#'                    condDat$ProjLML[popdefs[1,20]])
#'  str(pop1,max.level=1)
#' }
makeabpop <- function(popparam,midpts,projLML) {  #popparam=popdef;midpts=midpts;projLML=projLML
  #(DLMax,L50,L95,SigMax,SaMa,SaMb,Wta,Wtb,Me,L50C,L95C,R0 SelP[1],SelP[2],Nyrs,steeph,MaxCE,L50Mat block
  #   1    2   3    4      5    6   7   8   9 10   11   12  13     14      15   16     17    18     19
  numYr <- popparam["Nyrs"]
  N <- length(midpts)
  G <- STM(popparam[1:4],midpts)  # c(DLMax, L50, L95, SigMax)
  mature <- maturity(popparam["SaMa"],popparam["SaMb"],midpts)
  WtL <- WtatLen(popparam["Wta"],popparam["Wtb"],midpts)
  emergent <- logistic(popparam["L50C"],popparam["deltaC"],midpts)
  MatWt <- mature * WtL
  zSelect <- matrix(0,nrow=N,ncol=numYr,dimnames=list(midpts,1:numYr))
  blk <- popparam["block"]
  selL50 <- popparam["SelP1"]
  selL95 <- popparam["SelP2"]
  verLML <- unique(projLML) # find different LML in projections
  for (LML in verLML) {
    Sel <- logistic((LML+selL50),selL95,midpts)
    Sel <- Sel * emergent # include emergence to determine availability
    pick <- which(projLML == LML)
    zSelect[,pick] <- rep(Sel,length(pick))
  }
  zSelWt <- zSelect * WtL
  zLML <- projLML
  Me <- popparam["Me"]
  AvRec <- popparam["AvRec"]  # R0
  catq <- 0.0
  B0 <- 0.0
  ExB0 <- 0.0
  SaM <- -popparam["SaMa"]/popparam["SaMb"] # -a/b
  MSY <- 0
  MSYDepl <- 0
  bLML <- 0
  SMUname <- ""
  ans <- list(Me,AvRec,B0,ExB0,MSY,MSYDepl,bLML,catq,SaM,popparam,zLML,
              G,mature,WtL,emergent,zSelect,zSelWt,MatWt,SMUname)
  names(ans) <- c("Me","R0","B0","ExB0","MSY","MSYDepl","bLML","popq",
                  "SaM","popdef","LML","G","Maturity","WtL","Emergent",
                  "Select","SelWt","MatWt","SMUname")
  class(ans) <- "abpop"
  return(ans)
} # End of makeabpop

#' @title makeregionC makes the constant part of the simulation
#'
#' @description makeregionC makes the constant part of the simulated
#'     region. Once defined this does not change throughout the
#'     simulation. Once made it still requires makeregion to be run
#'     to fill in the B0, ExBo, MSY, MSYDepl, and the popq values, and
#'     to produce regionD, the dynamic part of the new region
#'
#' @param region the object derived from the readregionfile function
#' @param const the object derived from teh readdatafile function
#'
#' @return a list containing the constant part of the simulated region
#' @export
#'
#' @examples
#' print("Wait on built in data")
makeregionC <- function(region,const) { #  region=reg1; const=constants
  glb <- region$globals
  nSMU <- glb$nSMU
  numpop <- glb$numpop
  midpts <- glb$midpts
  blkdef <- region$SMUpop
  SMUindex <- defineBlock(nSMU,blkdef,numpop)
  SMUnames <- region$SMUnames
  if (region$randomseed > 0) set.seed(region$randomseed)
  projectionLML <- region$projLML
  historicalLML <- region$histLML
  if (region$condition)
    projLML <- historicalLML else projLML <- projectionLML
  # imph <- seq(0.01,uplim,0.01)   # harvest rate for productivity
  # numrow <- length(imph)
  pops <- trunc(const["popnum",])
  # columns <- c("ExB","MatB","AnnH","Catch","Deplet","RelCE")
  #  production <- array(0,dim=c(numpop,numrow,6),
  #                      dimnames=list(pops,imph,columns))
  popdefs <- definepops(nSMU,SMUindex,const,glob=glb) # define all pops in one go
  regionC <- vector("list",numpop)
  for (pop in 1:numpop) {      # pop <- 1
    popdef <- popdefs[pop,]
    regionC[[pop]] <- makeabpop(popdef,midpts,projLML)
    #  out <- doproduction(zone[[pop]],uplim=uplim)
    #  zone[[pop]]$MSY <- out$MSY
    #  zone[[pop]]$MSYDepl <- out$MSYDepl
    #   production[pop,,] <- out$Productivity
    tmpL <- oneyrgrowth(regionC[[pop]],regionC[[pop]]$SaM)
    regionC[[pop]]$bLML <- oneyrgrowth(regionC[[pop]],tmpL)
    regionC[[pop]]$SMUname <- SMUnames[SMUindex[pop]]
  }
  class(regionC) <- "regionC"
  ans <- list(regionC=regionC,popdefs=popdefs)
  return(ans)
}  # End of makeregionC


#' @title makeregion generates the parts of the simulated region
#'
#' @description makeregion generates the dynamics components of the
#'     simulated region and completes the constant components of the
#'     simulated region. The term 'region' refers to the upper level
#'     of geographical detail that is used. Thus, for example, in
#'     Tasmania we might simulate a number of statistical blocks with
#'     multiple populations. In combination, we refer to the total as
#'     a region. The matrices are Nyrs x numpop and the two
#'     arrays are N x Nyrs x numpop.
#'
#' @param glob the global constants defined for the current simulation.
#'     These include numpop, nblock, midptsd, Nclass, and Nyrs
#' @param regC the regionC object from makeregionC
#' @param uplim default=0.4, the upper bound on harvest rate used when
#'     estimation the production curves and statistic
#'
#' @return a list of the dynamics and the constant components of the
#'     simulation
#' @export
#'
#' @examples
#'  print("wait on internal data sets")
makeregion <- function(glob,regC,uplim=0.4) {
  Nyrs <- glob$Nyrs
  numpop <- glob$numpop
  N <- glob$Nclass
  ExplB <- matrix(0,nrow=Nyrs,ncol=numpop)
  MatB <- matrix(0,nrow=Nyrs,ncol=numpop)
  Catch <- matrix(0,nrow=Nyrs,ncol=numpop)
  Harvest <- matrix(0,nrow=Nyrs,ncol=numpop)
  cpue <- matrix(0,nrow=Nyrs,ncol=numpop)
  deplExB <- matrix(0,nrow=Nyrs,ncol=numpop)
  deplSpB <- matrix(0,nrow=Nyrs,ncol=numpop)
  Recruit <- matrix(0,nrow=Nyrs,ncol=numpop)
  CatchN <- array(data=0,dim=c(N,Nyrs,numpop))
  Nt <- array(data=0,dim=c(N,Nyrs,numpop))
  for (pop in 1:numpop) {  # pop=1
    SurvE <- exp(-regC[[pop]]$Me)
    recr <- rep(0,N)
    recr[1] <-regC[[pop]]$R0
    UnitM <- matrix(0,nrow=N,ncol=N)
    diag(UnitM) <- 1.0
    Minv <- solve(UnitM - (SurvE * regC[[pop]]$G ))  # (I - SG)-1
    Nt[,1,pop] <- Minv %*% recr          # [(I - SG)-1]R
    ExplB[1,pop] <- sum(regC[[pop]]$SelWt[,1]*Nt[,1,pop])/1e06
    regC[[pop]]$ExB0 <- ExplB[1,pop]
    deplExB[1,pop] <- 1.0  # no depletion when first generating regions
    MatB[1,pop] <- sum(regC[[pop]]$MatWt*Nt[,1,pop])/1e06
    regC[[pop]]$B0 <- MatB[1,pop]
    deplSpB[1,pop] <- 1.0
    Recruit[1,pop] <- recr[1]
    regC[[pop]]$popq <- regC[[pop]]$popdef["MaxCE"]/ExplB[1,pop]
    cpue[1,pop] <- 1000.0 * regC[[pop]]$popq * ExplB[1,pop]
  }
  ans <- list(matureB=MatB,exploitB=ExplB,catch=Catch,
              harvestR=Harvest,cpue=cpue,recruit=Recruit,
              deplsB=deplSpB,depleB=deplExB,catchN=CatchN,Nt=Nt)
  return(list(regionD=ans,regionC=regC))
} # end of makeregionD



#' @title maturity Logistic maturity curve
#'
#' @description maturity this uses the logistic function:
#'   exp(a+b*L)/(1+exp(a+b*L)), which has the property that the SM50 = -a/b
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

#' @title oneyear do one year's dynamics for one input population abpop
#'
#' @description oneyear do one year's dynamics for one input population.
#'    Used to step through populations of a region. Its dynamics are
#'    such that it first undergoes growth, then half natural mortality.
#'    This allows an estimate of exploitable biomass before fishing
#'    occurs. The remaining dynamics involve the removal of the catch,
#'    the application of the last half of natural mortality and the
#'    addition of recruits. Which allows the exploitable biomass to be
#'    estimated after fishing. The recruitment occurs in oneyearD so
#'    that larval dispersal can be accounted for.
#'
#' @param inpopC a single population from a regionC, an abpop
#' @param inNt the numbers at size for the year previous to the year
#'     of dynamics. These are projected into the active year.
#' @param Nclass the number of size classes used to describe growth.
#'     used to define vectors
#' @param inCat a literal catch in tonnes to be removed during the
#'     year, a scalar
#' @param yr the year in the dynamics being worked on. The first year
#'     is generated when the zone is defined or when it is initially
#'     depleted. All dynamics are appllied from year 2 - Nyrs; scalar
#' @param calcharv should we calculate the harvest rate from inCat.
#'     Default=TRUE. If FALSE, implies that inCat is a harvest rate
#'
#' @return a list containing ExploitB, MatureB, Catch, Harvest, Nt,
#'     ce, and CatchN used to update the given pop in yr + 1
#' @export
#'
#' @examples
#' print("need to wait on built in data sets")
oneyear <- function(inpopC,inNt,Nclass,inCat,yr,calcharv=TRUE) {  #
  #  yr=2; pop=1; inpopC=regionC[[pop]]; inNt=regionD$Nt[,yr-1,pop]; Nclass=glb$Nclass; inCat=0.0;
  MatWt <- inpopC$MatWt/1e06
  SelectWt <- inpopC$SelWt[,yr]/1e06
  selyr <- inpopC$Select[,yr]
  Ne <- numeric(Nclass)
  Cat <- numeric(Nclass)
  Os <- exp(-inpopC$Me/2)
  MatureB <- sum(MatWt*inNt)
  NumNe <- (Os * (inpopC$G %*% inNt))
  ExploitB <- sum(SelectWt * NumNe) #SelectWt=Select*WtL
  oldExpB <- ExploitB   # ExploitB after growth and 0.5NatM
  if (calcharv) {
    Ht <- inCat/ExploitB
    if (Ht > 0.9) {  # can I eliminate this step? Perhaps use logits.
      Ht <- 0.9
      warning(paste0("Harvest rates > 0.9 in year ",yr," in pop ",
                    (inpopC$popdef["block"])))
    }
  } else {
    Ht = inCat
  }
  Fish <- 1-(Ht*selyr)
  newNt <- (Os * (Fish * NumNe)) #+ Rec # Nt - catch - 0.5M, and + Rec
  Cat <- (Ht*selyr) * NumNe  #numbers at size in the catch
  ExploitB <- sum(SelectWt * newNt)
  MatureB <- sum(MatWt*newNt) #+ MatBC
  Catch <- sum(inpopC$WtL*Cat)/1e06
  Harvest <- (2.0 * Catch)/(oldExpB + ExploitB)  # average of the start and end
  ce <- inpopC$popq * ((oldExpB + ExploitB)/2) * 1000.0  #ExploitB
  ans <- list(ExploitB,MatureB,Catch,Harvest,newNt,ce,Cat)
  names(ans) <- c("ExploitB","MatureB","Catch","Harvest","Nt","ce",
                  "CatchN")
  return(ans)
} # End of oneyear


#' @title oneyearD conducts one year's dynamics in the simulation
#'
#' @description onyearD conducts one year's dynamics in the simulation
#'     returning the revised regionD, which will have had a single year
#'     of activity included in each of its components.
#'
#' @param regC the constant portion of the region with a list of
#'     properties for each population
#' @param regD the dynamics portion of the region, with matrices and
#'     arrays for the dynamic variables of the dynamics of the
#'     operating model
#' @param Ncl the number of size classes used to describe size
#' @param incatch a vector of catches taken in the year from each
#'     population
#' @param year the year of the dynamics, would start in year 2 as year
#'     1 is the year of initiation.
#' @param sigmar the variation in recruitment dynamics, set to 1e-08
#'     when searching for equilibria.
#' @param npop the number of populations, the global numpop
#' @param deltarec the rate of larval dispersal
#' @param calcharv should we calculate the harvest rate from inCat.
#'     Default=TRUE. If FALSE, implies that inCat is a harvest rate
#'
#' @return a list containing a revised dynamics list
#' @export
#'
#' @examples
#' print("wait for built in data sets")
oneyearD <- function(regC,regD,Ncl,incatch,year,sigmar,npop,deltarec,
                     calcharv=TRUE) {
  # npop=6; regC=regionC; regD=regionD; Ncl=glb$Nclass; incatch=rep(0.0,npop); year=2; sigmar=1e-08; deltarec=reg1$larvdisp
  matb <- numeric(npop)
  for (popn in 1:npop) {  #  popn=1
    out <- oneyear(inpopC=regC[[popn]],inNt=regD$Nt[,year-1,popn],
                   Nclass=Ncl,inCat=incatch[popn],yr=year,
                   calcharv=calcharv)
    regD$exploitB[year,popn] <- out$ExploitB
    regD$matureB[year,popn] <- out$MatureB
    regD$catch[year,popn] <- out$Catch
    regD$harvestR[year,popn] <- out$Harvest
    regD$cpue[year,popn] <- out$ce
    regD$Nt[,year,popn] <- out$Nt
    matb[popn] <- out$MatureB
  }
  steep <- sapply(regC,"[[","popdef")["steeph",]
  r0 <- sapply(regC,"[[","R0")
  b0 <- sapply(regC,"[[","B0")
  recs <- oneyearrec(steep,r0,b0,matb,sigR=sigmar)
  newrec <- driftrec(recs,deltarec)
  regD$recruit[year,] <- newrec
  regD$Nt[1,year,] <- newrec
  regD$deplsB[year,] <- regD$matureB[year,]/b0
  regD$depleB[year,] <- regD$exploitB[year,]/sapply(regC,"[[","ExB0")
  return(regD)
} # end of oneyearD

#' @title oneyearrec calculates the Beverton-Holt recruitment
#'
#' @description oneyearrec calculates the Beverton-Holt recruitment for a
#'    single population in a single year; parameterized with steepness,
#'    R0, and B0. Includes the logical parameter 'recvar' which determines
#'    whether recruitment variation is expressed or not.
#' @param steep the steepness for the population; scalar
#' @param R0 the unfished recruitment levels for the population; scalar
#' @param B0 the unfished spawning biomass; scalar
#' @param Bsp the current spawning biomass; scalar
#' @param sigR standard deviation of the recruitment residuals;
#'     scalar. set this to 1e-08 to avoid recruitment variability
#' @return an absolute number of recruits from a given spawning biomass
#' @export
#'
#' @examples
#' \dontrun{
#' steep <- 0.707
#' R0 <- 319971
#' B0 <- 313
#' Bsp <- 147
#' insigmar <- 0.3
#' oneyearrec(steep,R0,B0,Bsp,insigmar)
#' oneyearrec(steep,R0,B0,Bsp,insigmar)
#' }
oneyearrec <- function(steep,R0,B0,Bsp,sigR) {
   epsilon <- exp(rnorm(1,mean=0,sd=sigR) - (sigR * sigR)/2)
   rec <- ((4*steep*R0*Bsp)/((1-steep)*B0+(5*steep-1)*Bsp)) * epsilon
   return(rec)
} # end of oneyearrec

#' @title oneyrgrowth one years' growth for a population and initial size
#'
#' @description oneyrgrowth one years' growth for a given population and
#'    initial size. Used to determine the size after two year's growth as
#'    a measure of the expected bLML - biological Legal Minimum Length. A
#'    reflection of the two year rule in Tasmania. To get this the function
#'    should obviously be run twice.
#' @param inpop the abpop to be grown forward
#' @param startsize default = 2, but often set to the size at 50 percent
#'    maturity or the vector of midpts to determine the growth of Nt
#' @return the expected mean size after one year's growth
#' @export
#'
#' @examples
#' \dontrun{
#' nblock <- 2
#' filename <- datafiletemplate(numblock=nblock,filename="block3.csv")
#' condDat <- readdatafile(filename)
#' bConst <- condDat$constants
#' blkpop <- condDat$blkpop
#' numpop <- condDat$numpop
#' blockI <- defineBlock(nblock,blkpop,numpop)
#' popdefs <- definepops(nblock,blockI,bConst)
#' pop1 <- makeabpop(popdefs[1,],condDat$midpts,
#'                   condDat$ProjLML[popdefs[1,20]])
#' oneyrgrowth(pop1,startsize=70)
#' }
oneyrgrowth <- function(inpop,startsize=2) {
  pdef <- inpop$popdef
  MaxDL <- pdef[1]
  L50 <- pdef[2]
  L95 <- pdef[3]
  delta <- MaxDL/(1+exp(log(19)*(startsize - L50)/(L95-L50)))
  outsize <- startsize + delta
  names(outsize) <- NULL
  return(outsize)
} # end of onyrgrowth


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
  cat("B0                   : ",round(sum(getlistVar(x,"B0")),3), "\n")
  cat("MSY                  : ",round(sum(getlistVar(x,"MSY")),3), "\n")
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
  Nyrs <- glob$Nyrs
  midpts <- glob$midpts
  if (length(inLML) != glob$nblock)
    stop("the number of LML does not match the number of blocks")
  for (pop in 1:Npop) {  # inzone <- zone; pop <- 1
    popdef <- inzone[[pop]]$popdef
    blk <- popdef["block"]
    Sel <- logistic((inLML[blk]+popdef["SelP1"]),
                    (inLML[blk]+popdef["SelP2"]),midpts)
    inzone[[pop]]$Select[,inyear:Nyrs] <- Sel
    inzone[[pop]]$SelWt[,inyear:Nyrs] <- Sel * inzone[[pop]]$WtL
    inzone[[pop]]$LML[inyear:Nyrs] <- inLML[blk]
  }
  return(inzone)
}  # end of resetLML

#' @title restart transfers final year values of regD into the first year
#'
#' @description restart transfers the final year values from the
#'     dynamics part of the region (regionD), into the first year.
#'     This is used, for example, when searching for an equilibrium
#'     state if there is larval dispersal > 0.0. Of if one sets the
#'     initial depletion to anything other than 1.0. Contains the
#'     option of setting every other cell to zero, which is the
#'     default.
#'
#' @param oldregD the old regionD contianing the dynamics as run for
#'     Nyrs.
#' @param nyrs The number of years of dynamics, the global Nyrs
#' @param npop The number of populations, the global numpop
#' @param N the number of size classes, teh global Nclass
#' @param zero should the arrays be otherwise filled with zeros? The
#'     default = TRUE
#'
#' @return a list containing a revised dynamics list
#' @export
#'
#' @examples
#' print("wait for built in data sets")
restart <- function(oldregD,nyrs,npop,N,zero=TRUE) { # oldregD=regionD; nyrs=Nyrs; npop=npop; N=Nclass
  ExplB <- matrix(0,nrow=nyrs,ncol=npop)
  MatB <- matrix(0,nrow=nyrs,ncol=npop)
  Catch <- matrix(0,nrow=nyrs,ncol=npop)
  Harvest <- matrix(0,nrow=nyrs,ncol=npop)
  cpue <- matrix(0,nrow=nyrs,ncol=npop)
  deplExB <- matrix(0,nrow=nyrs,ncol=npop)
  deplSpB <- matrix(0,nrow=nyrs,ncol=npop)
  Recruit <- matrix(0,nrow=nyrs,ncol=npop)
  CatchN <- array(data=0,dim=c(N,nyrs,npop))
  Nt <- array(data=0,dim=c(N,nyrs,npop))
  regD <- list(matureB=MatB,exploitB=ExplB,catch=Catch,
               harvestR=Harvest,cpue=cpue,recruit=Recruit,
               deplsB=deplSpB,depleB=deplExB,catchN=CatchN,Nt=Nt)
  if (!zero) regD <- oldregD
  regD$matureB[1,] <- oldregD$matureB[nyrs,]
  regD$exploitB[1,] <- oldregD$exploitB[nyrs,]
  regD$catch[1,] <- oldregD$catch[nyrs,]
  regD$harvestR[1,] <- oldregD$harvestR[nyrs,]
  regD$cpue[1,] <- oldregD$cpue[nyrs,]
  regD$recruit[1,] <- oldregD$recruit[nyrs,]
  regD$deplsB[1,] <- oldregD$deplsB[nyrs,]
  regD$depleB[1,] <- oldregD$depleB[nyrs,]
  regD$catchN[,1,] <- oldregD$catchN[,nyrs,]
  regD$Nt[,1,] <- oldregD$Nt[,nyrs,]
  return(regD)
} # end of restart

#' @title runthree conducts the dynamics with constant catch 3 times
#'
#' @description runthree is used when searching numerically for an
#'     equilibrium and it conducts the Nyrs dynamics three times, each
#'     time through it replaces year 1 with year Nyrs. Thus if Nyrs is
#'     40 it conducts 3 * 39 years of dynamics (117 years).
#'
#' @param regC the constant portion of the region with a list of
#'     properties for each population
#' @param regD the dynamics portion of the region, with matrices and
#'     arrays for the dynamic variables of the dynamics of the
#'     operating model
#' @param glob the globals variable from readregionfile
#' @param catch a vector of numpop catches to be held constant
#'     across all years.
#'
#' @return a list containing a revised dynamics list, regionD
#' @export
#'
#' @examples
#' print("wait on built in data sets")
runthree <- function(regC,regD,glob,catch) {
  npop <- glob$numpop
  Nclass <- glob$Nclass
  Nyrs <- glob$Nyrs
  larvdisp <- glob$larvdisp
  for (iter in 1:3) {
    for (yr in 2:Nyrs)
      regD <- oneyearD(regC=regC,regD=regD,Ncl=Nclass,
                       incatch=catch,year=yr,sigmar=1e-08,npop=npop,
                       deltarec=larvdisp)
    regD <- restart(oldregD=regD,nyrs=Nyrs,npop=npop,N=Nclass,zero=TRUE)
  }
  return(regD)
} # end of runthree

#' @title STM Generates the Size Transition Matrix for Inverse Logistic
#'
#' @description STM With the input of the four parameters inside a vector,
#'     and a vector of initial lengths or mid-points of size classes STM
#'     generates a square transition matrix with the probabilities of
#'     growing from each initial size into the same or larger sizes.
#'     Negative growth is disallowed. All columns in the matrix sunm to one.
#' @param p a vector of four parameters in the following order
#'     MaxDL the maximum growth increment of the inverse logistic,
#'     L50 the initial length at which the inflexion of the growth
#'     increment curve occurs, L95 - the initial length that defines the
#'     95th percentile of the growth increments, SigMax - the maximum
#'     standard deviaiton of the normal distribution used to describe the
#'     spread of each distribution of growth increments
#' @param mids a vector of initial lengths which also define the width of
#'     each size class thus, mids from 2 - 210 woul dimply 2mm size
#'     classes 1 - 3 = 2, 3 - 5 = 4, etc
#' @return A square matrix with dimension =the length of the mids vector
#' @export
#'
#' @references Haddon, M., Mundy, C., and D. Tarbath (2008) Using an
#'     inverse-logistic model to describe growth increments of blackip
#'     abalone (Haliotis rubra) in Tasmania. Fisheries Bulletin 106: 58-71
#' @examples
#' \dontrun{
#'  param <- c(25.0,120.0,170.0,4.0)
#'  midpts <- seq(2,210,2)
#'  G <- STM(param,midpts)
#'  print(round(G[1:30,1:8],4))
#' }
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
   rownames(G) <- mids
   colnames(G) <- mids
   class(G) <- "STM"
   return(G)
} # end of STM

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

