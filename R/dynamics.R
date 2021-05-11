

#' @title depleteSAU resets zoneD, approximately to an input depletion level
#'
#' @description depleteSAU resets the depletion level of the whole
#'     zone and does this by searching for the harvest rate that
#'     leads to the mature biomass in each population is as close as
#'     possible to the desired depletion level. This means
#'     the individual populations will likely vary around the target
#'     depletion, so the depletion across the zone will only be
#'     approximately at the target depletion. The
#'     depletion is measured relative to the effective B0 as that
#'     takes account of any larval dispersal. This function uses the
#'     production curve array to search for harvest rates that bound
#'     the target depletion and then re-searches across those bounds
#'     using len intervals.
#'
#' @param zoneC the constants components of the simulated zone
#' @param zoneD the dynamic components of the simulated zone
#' @param glob the general global variables
#' @param initdepl a vector of target depletion levels for each SAU. This
#'     is found in zone1$condC, but could obviously be modified in each run.
#' @param product the production curve matrix from doproduction
#' @param len the number of intervals in the trial harvest rates
#'
#' @return a revised zoneD object
#' @export
#'
#' @examples
#' \dontrun{
#'   data(zone)
#'   glb <- zone$glb
#'   depl <- rep(0.3,glb$nSAU)
#'   zoneDD <- depleteSAU(zone$zoneC,zone$zoneD,glb,initdepl=depl,zone$product)
#'   sum((zoneDD$matureB[1,]/sum(zoneDD$matureB[1,]))*zoneDD$deplsB[1,])
#'   mean(zoneDD$deplsB[1,])
#' }  # zoneC=zone$zoneC;zoneD=zone$zoneD;glob=zone$glb;initdepl=origdepl;product=zone$product;len=12
depleteSAU <- function(zoneC,zoneD,glob,initdepl,product,len=15) {
  #  zoneC=zoneC; zoneD=zoneD; glob=glb;  product=product; depl=initdepl;len=15
  # use product to find bounds on H around the desired depletion level
  if (length(initdepl) != glob$nSAU) stop("Need a depletion for each population \n")
  npop <- glob$numpop  # need to deplete by pop
  depl <- initdepl[glob$sauindex] # each pop should have depl of its SAU
  harvests <- matrix(0,nrow=len,ncol=npop,dimnames=list(1:len,1:npop))
  popDepl <- harvests # need a matrix to hold the depletion levels
  initH <- as.numeric(rownames(product))
  # bound the harvest rates nearest to the desired depletion
  for (pop in 1:npop) {
    pick <- which.closest(depl[pop],product[,"Deplet",pop])
    if (pick == nrow(product)) {
      mssg <- paste0("Initial maximum harvest rate when estimating ",
                     "productivity was too low to generate desired ",
                     "depletion level")
      stop(mssg)
    }
    if (pick == 1)
      stop("Do not use depleteSAU with a depletion level of 1.0 \n")
    lowl <- initH[pick-1]
    upl <- initH[pick+1]
    harvests[,pop] <- seq(lowl,upl,length=len)
  }
  for (harv in 1:len) { # harv=1
    zoneDD <- runthreeH(zoneC=zoneC,zoneD,inHarv=harvests[harv,],glob)
    popDepl[harv,] <- zoneDD$deplsB[1,]
  }
  pickharv <- numeric(npop)
  for (pop in 1:npop) {
    pickdepl <- which.closest(depl[pop],popDepl[,pop],index=TRUE)
    pickharv[pop] <- harvests[pickdepl,pop]
  }
  zoneDD <- runthreeH(zoneC=zoneC,zoneD,inHarv=pickharv,glob)
  return(zoneDD)
} # end of depleteSAU

#' @title dohistoricC imposes the historical catches on an unfished zone
#'
#' @description dohistoricC is used during the conditioning of the zone/region
#'     and imposes the historical catches, held in the zone data object
#'     obtained using readzonefile, onto an unfished initial zone, and it does
#'     this by imposing the time-series of catches to each SAU/block. The
#'     operation is through the use of oneyearC, which imposes one year's
#'     catches, which are a vector of SAU catches for each year of the series.
#'
#' @param zoneDD The input unfished dynamic zone, zoneD, object
#' @param zoneC the zone constants object, zoneC
#' @param glob the globals object
#' @param condC from the zone1 object, contains historical fisheries data plus
#' @param calcpopC a function that takes the output from hcrfun and gernerates
#'     the actual catch per population expected in the current year.
#' @param sigR the recruitment variation included
#' @param sigB the variability introduced to the catches by population by
#'     fishers not knowing the distribution of exploitable biomass exactly. What
#'     this value should be is unknown, the default=1e-08, is arbitrary but
#'     avoids any effective fisher allocation error between populations.
#'
#' @return a zoneD object
#' @export
#'
#' @examples
#' print("wait on some data sets")
dohistoricC <- function(zoneDD,zoneC,glob,condC,calcpopC,sigR=1e-08,sigB=1e-08) {
  # zoneC=zone$zoneC; zoneDD=zone$zoneD;glob=zone$glb;condC=zone$zone1$condC
  # sigR=1e-08; sigB=1e-08; Ncl=glob$Nclass;sauindex=glob$sauindex;movem=glob$move;
  histC <- condC$histCatch
  yrs <- condC$histyr[,"year"]
  nyrs <- length(yrs)
  sauindex <- glob$sauindex
  r0 <- getvar(zoneC,"R0") #sapply(zoneC,"[[","R0")
  b0 <- getvar(zoneC,"B0") #sapply(zoneC,"[[","B0")
  exb0 <- getvar(zoneC,"ExB0")
  for (year in 2:nyrs) {  # year=2 # ignores the initial unfished year
    catchsau <- histC[year,]
    hcrout <- list(acatch=catchsau)
    popC <- calcpopC(hcrout,exb=zoneDD$exploitB[year-1,],sauindex,sigmab=sigB)
    inN <- zoneDD$Nt[,year-1,]
    out <- oneyearsauC(zoneCC=zoneC,inN=inN,popC=popC,year=year,
                       Ncl=glob$Nclass,sauindex=sauindex,movem=glob$move,
                       sigmar=sigR,r0=r0,b0=b0,exb0=exb0)
    dyn <- out$dyn
    zoneDD$exploitB[year,] <- dyn["exploitb",]
    zoneDD$midyexpB[year,] <- dyn["midyexpB",]
    zoneDD$matureB[year,] <- dyn["matureb",]
    zoneDD$catch[year,] <- dyn["catch",]
    zoneDD$harvestR[year,] <- dyn["catch",]/out$dyn["exploitb",]
    zoneDD$cpue[year,] <- dyn["cpue",]
    zoneDD$recruit[year,] <- dyn["recruits",]
    zoneDD$deplsB[year,] <- dyn["deplsB",]
    zoneDD$depleB[year,] <- dyn["depleB",]
    zoneDD$Nt[,year,] <- out$NaL
    zoneDD$catchN[,year,] <- out$catchN
  }
  return(zoneDD)
} # end of dohistoricC

#' @title imperr calculates population catches from sau catches with error
#'
#' @description imperr converts aspirational sau catches into population level
#'     catches while introducing management implementation error. Here this
#'     error is implemented as Log-Normal errors on diver intuitions
#'     concerning the relative abundance in each population. The error is
#'     imposed separately on the populations in each SAU.
#'
#' @param catchsau the predicted or aspirational catch per SAU from the Harvest
#'     control rule
#' @param exb the exploitable biomass at the end of the previous year. In the
#'     first year of the projections this would be the last year of the
#'     conditioning.
#' @param sauindex the SAU index for each population
#' @param sigmab the Log-Normal standard deviation of implementation error. The
#'     default value = 1e-08, which effectively means no errors.
#'
#' @return a vector of population catches for the year to be imposed after the
#'     estimation of exploitable biomass
#' @export
#'
#' @examples
#' print("wait on suitable internal data sets")
imperr <- function(catchsau,exb,sauindex,sigmab=1e-08) {
  # exb=zoneDD$exploitB[endyr,]; suaindex=glb$sauindex; sigmaB=0.1; catchsau=acatch
  TAC <- sum(catchsau,na.rm=TRUE)
  totexb <- sum(exb,na.rm=TRUE)
  npop <- length(exb)
  nSAU <- length(catchsau)
  sauexb <- tapply(exb,sauindex,sum,na.rm=TRUE) * exp(rnorm(nSAU,mean=0,sd=sigmab))
  sauexb <- sauexb * (totexb/sum(sauexb,na.rm=TRUE))
  popC <- catchsau[sauindex] * (exb/sauexb[sauindex])
  popC <- popC * TAC/sum(popC)
  return(popC)
} # end of imperr

#' @title oneyear one year's harvest dynamics for one abpop
#'
#' @description oneyear do one year's dynamics for one input population.
#'    Used to step through populations of a zone. Its dynamics are
#'    such that it first undergoes growth, then half natural mortality.
#'    This allows an estimate of exploitable biomass before fishing
#'    occurs. The remaining dynamics involve the removal of the catch,
#'    the application of the last half of natural mortality and the
#'    addition of recruits. Which allows the exploitable biomass to be
#'    estimated after fishing. The recruitment occurs in oneyearD so
#'    that larval dispersal can be accounted for. Thus, oneyear requires
#'    oneyearD.
#'
#' @param inpopC a single population from a zoneC, an abpop
#' @param inNt the numbers at size for the year previous to the year
#'     of dynamics. These are projected into the active year.
#' @param Nclass the number of size classes used to describe growth.
#'     used to define vectors
#' @param inH a literal annual harvest rate as a proportion to be
#'     removed as catch during the year, a scalar
#' @param yr the year in the dynamics being worked on. The first year
#'     is generated when the zone is defined or when it is initially
#'     depleted. All dynamics are applied from year 2 - Nyrs; scalar
#'
#' @return a list containing ExploitB, MatureB, Catch, Harvest, Nt,
#'     ce, CatchN, and midyexpB used to update the given pop in yr + 1
#' @export
#'
#' @examples
#' print("need to wait on built in data sets")
oneyear <- function(inpopC,inNt,Nclass,inH,yr) {  #
  # yr=2; pop=2; inpopC=zoneC[[pop]]; inNt=zoneD$Nt[,yr-1,pop];
  # Nclass=glb$Nclass; inH=0.2;
  MatWt <- inpopC$MatWt/1e06
  SelectWt <- inpopC$SelWt[,yr]/1e06
  selyr <- inpopC$Select[,yr]
  Ne <- numeric(Nclass)
  Cat <- numeric(Nclass)
  Os <- exp(-inpopC$Me/2)
  #  MatureB <- sum(MatWt*inNt)
  NumNe <- (Os * (inpopC$G %*% inNt))
  ExploitB <- sum(SelectWt * NumNe) #SelectWt=Select*WtL
  oldExpB <- ExploitB   # oldExpB = ExploitB after growth and 0.5NatM
  Fish <- 1-(inH*selyr)
  newNt <- (Os * (Fish * NumNe)) #+ Rec # Nt - catch - 0.5M,
  Cat <- (inH*selyr) * NumNe  #numbers at size in the catch
  newExpB <- sum(SelectWt * newNt)
  avExpB <- (newExpB + oldExpB)/2.0 #av start and end
  MatureB <- sum(MatWt*newNt) #+ MatBC
  Catch <- sum(inpopC$WtL*Cat)/1e06
 # Harvest <- Catch/avExpB  # uses average of the start and end
  ce <- inpopC$popq * avExpB * 1000.0  # Need to add CPUE variation
  ans <- list(newExpB,MatureB,Catch,inH,newNt,ce,Cat,oldExpB)
  names(ans) <- c("ExploitB","MatureB","Catch","Harvest","Nt","ce",
                  "CatchN","midyexpB")
  return(ans)
} # End of oneyear

#' @title oneyearcat one year's catch dynamics for one abpop
#'
#' @description oneyear does one year's dynamics for one input
#'     population. Where the fishing is based on a given catch per
#'     population, which would be determined by any harvest control
#'     rule. It is used to step through populations of a zone.
#'     Its dynamics are such that each population first undergoes
#'     growth, then half natural mortality. This allows an estimate of
#'     exploitable biomass before fishing occurs. The remaining
#'     dynamics involve the removal of the catch, after estimating the
#'     harvest rate from catch/exploitb, the application of the last
#'     half of natural mortality and the addition of recruits. This
#'     allows the exploitable biomass to be estimated after fishing.
#'     The cpue uses the average of the pre-fishing and post-fishing
#'     exploitB to smooth over any changes brought about by fishing.
#'     The recruitment occurs in oneyearC so that larval dispersal can
#'     be accounted for. oneyearcat is not used independently of
#'     oneyearC.
#'
#' @param inpopC a single population from a zoneC, an abpop
#' @param inNt the numbers at size for the year previous to the year
#'     of dynamics. These are projected into the active year.
#' @param Nclass the number of size classes used to describe growth.
#'     used to define vectors
#' @param incat the literal annual catch from the population. Derived
#'     from teh harvest control rule for each SAU, then allocated to
#'     each population with respect to its relative catching
#'     performance.
#' @param yr the year in the dynamics being worked on. The first year
#'     is generated when the zone is defined or when it is initially
#'     depleted. All dynamics are appllied from year 2 - Nyrs; scalar
#'
#' @return a list containing a vector of ExploitB, MatureB, Catch, and ce, and
#'     a matrix NaL containing Nt and CatchN used to update a pop in yr + 1
#' @export
#'
#' @examples
#' print("need to wait on built in data sets")
oneyearcat <- function(inpopC,inNt,Nclass,incat,yr) {  #
  # yr=2; pop=2; inpopC=zoneC[[pop]]; inNt=zoneD$Nt[,yr-1,pop];
  # Nclass=glb$Nclass; inH=0.05;
  MatWt <- inpopC$MatWt/1e06
  SelectWt <- inpopC$SelWt[,yr]/1e06
  selyr <- inpopC$Select[,yr]
  Ne <- numeric(Nclass)
  Cat <- numeric(Nclass)
  Os <- exp(-inpopC$Me/2)
  NumNe <- (Os * (inpopC$G %*% inNt))
  midyexpB <- sum(SelectWt * NumNe) #SelectWt=Select*WtL =midyrexB
  estH <- min(incat/midyexpB,0.8) # no more than 0.8 harvest rate
  Fish <- 1-(estH*selyr)
  newNt <- (Os * (Fish * NumNe))
  Cat <- (estH*selyr) * NumNe  #numbers at size in the catch
  ExploitB <- sum(SelectWt * newNt) # end of year exploitable biomass
  avExpB <- (midyexpB + ExploitB)/2.0 #av start and end
  MatureB <- sum(MatWt*newNt)
  Catch <- sum(inpopC$WtL*Cat)/1e06
  ce <- inpopC$popq * avExpB * 1000.0  #ExploitB
  vect <- c(exploitb=ExploitB,midyexpB=midyexpB,matureb=MatureB,
            catch=Catch,cpue=ce)
  ans <- list(vect=vect,NaL=newNt,catchN=Cat,NumNe=NumNe)
  return(ans)
} # End of oneyearcat


#' @title oneyearsauC conducts one year's dynamics using catch not harvest
#'
#' @description oneyearsauC conducts one year's dynamics in the simulation
#'     using historical SAU catches rather than harvest rates. The harvest rates
#'     are estimated after first estimating the exploitable biomass.
#'     returning the revised zoneD, which will have had a single year
#'     of activity included in each of its components. uses the function
#'     imperr to introduce implementation error to the actual catches.
#'
#' @param zoneCC the constant portion of the zone with a list of
#'     properties for each population
#' ##param exb the vector of exploitable biomass for each population for the
#' ##     previous year
#' @param inN the numbers-at-length for each population for the previous year
#' @param popC a vector of actual catches to be taken in the year from each SAU
#' @param year the year of the dynamics the exb and inN would be from year-1
#' @param Ncl the number of size classes used to describe size, global Nclass
#' @param sauindex the sau index for each population from glb
#' @param movem the larval dispersal movement matrix, global move
#' @param sigmar the variation in recruitment dynamics, set to 1e-08
#'     when searching for an equilibrium.
#' @param r0 the unfished R0 level used by oneyearrec
#' @param b0 the unfished mature biomass B0, used in recruitment and depletion
#' @param exb0 the unfished exploitable biomass used in depletion
#'
#' @return a list containing a revised dynamics list
#' @export
#'
#' @examples
#' print("Wait on new data")
oneyearsauC <- function(zoneCC,inN,popC,year,Ncl,sauindex,
                        movem,sigmar=1e-08,r0,b0,exb0) {
 # zoneCC=zoneCP;exb=zoneDP$exploitB[year-1,,iter]
#  inN=zoneDP$Nt[,year-1,,iter];catchsau=acatch;year=year
#  Ncl=Nclass;sauindex=sauindex;movem=movem
#  sigmar=sigmar;sigmab=sigmab
#  popC <- imperr(catchsau,exb,sauindex,sigmab)
  npop <- length(popC)
  matb <- numeric(npop)
  ans <- vector("list",npop)
  for (popn in 1:npop) {  # popn=14
    ans[[popn]] <- oneyearcat(inpopC=zoneCC[[popn]],inNt=inN[,popn],
                              Nclass=Ncl,incat=popC[popn],yr=year)
  }
  dyn <- sapply(ans,"[[","vect")
  steep <- getvect(zoneCC,"steeph") #sapply(zoneC,"[[","popdef")["steeph",]
  recs <- oneyearrec(steep,r0,b0,dyn["matureb",],sigR=sigmar)
  recruits <- as.numeric(movem %*% recs)
  deplsB <- dyn["matureb",]/b0
  depleB <- dyn["exploitb",]/exb0 # end of year exploitB
  dyn <- rbind(dyn,recruits,deplsB,depleB)
  NaL <-  sapply(ans,"[[","NaL")
  NaL[1,] <- recruits
  catchN <- sapply(ans,"[[","catchN")
  NumNe <- sapply(ans,"[[","NumNe")
  return(list(dyn=dyn,NaL=NaL,catchN=catchN,NumNe=NumNe))
} # end of oneyearsauC

#' @title oneyearD conducts one year's dynamics on zoneD in the MSE
#'
#' @description oneyearD conducts one year's dynamics on zoneD in the MSE
#'     returning the revised zoneD, which will have had a single year
#'     of activity included in each of its components. This uses zoneC
#'     but always within the environment of another function in which
#'     zoneC (as zoneC) can be found. Used in runthreeH, (and hence
#'     dodepletion and doproduction) and in testequil.
#'
#' @param zoneC the constant portion of the zone with a list of
#'     properties for each population
#' @param zoneD the dynamics portion of the zone, with matrices and
#'     arrays for the dynamic variables of the dynamics of the
#'     operating model
#' @param inHt a vector of harvest rates taken in the year from each
#'     population
#' @param year the year of the dynamics, would start in year 2 as year
#'     1 is the year of initiation.
#' @param sigmar the variation in recruitment dynamics, set to 1e-08
#'     when searching for equilibria.
#' @param Ncl the number of size classes used to describe size
#' @param npop the number of populations, the global numpop
#' @param movem the larval dispersal movement matrix
#'
#' @return a list containing a revised dynamics list
#' @export
#'
#' @examples
#'  print("wait oin revised data")
oneyearD <- function(zoneC,zoneD,inHt,year,sigmar,Ncl,npop,movem) {
#  zoneC=zoneC;zoneD=zoneD;Ncl=Nclass;inHt=inHarv;year=yr;sigmar=1e-08;npop=npop;movem=glob$move
  matb <- numeric(npop)
  for (popn in 1:npop) {  # year=2; popn=1
    out <- oneyear(inpopC=zoneC[[popn]],inNt=zoneD$Nt[,year-1,popn],
                   Nclass=Ncl,inH=inHt[popn],yr=year)
    zoneD$exploitB[year,popn] <- out$ExploitB
    zoneD$midyexpB[year,popn] <- out$midyexpB
    zoneD$matureB[year,popn] <- out$MatureB
    zoneD$catch[year,popn] <- out$Catch
    zoneD$harvestR[year,popn] <- out$Harvest
    zoneD$cpue[year,popn] <- out$ce
    zoneD$Nt[,year,popn] <- out$Nt
    zoneD$catchN[,year,popn] <- out$CatchN
    matb[popn] <- out$MatureB
  }
  steep <- getvect(zoneC,"steeph") #sapply(zoneC,"[[","popdef")["steeph",]
  r0 <- getvar(zoneC,"R0") #sapply(zoneC,"[[","R0")
  b0 <- getvar(zoneC,"B0") #sapply(zoneC,"[[","B0")
  recs <- oneyearrec(steep,r0,b0,matb,sigR=sigmar)
  newrecs <- movem %*% recs
  zoneD$recruit[year,] <- newrecs
  zoneD$Nt[1,year,] <- newrecs
  zoneD$deplsB[year,] <- zoneD$matureB[year,]/b0
  zoneD$depleB[year,] <- zoneD$exploitB[year,]/getvar(zoneC,"ExB0")
  return(zoneD)
} # end of oneyearD   round(zoneD$Nt[,year,])


#' @title oneyearrec calculates the Beverton-Holt recruitment
#'
#' @description oneyearrec calculates the Beverton-Holt recruitment for the
#'    input populations in a single year; parameterized with steepness,
#'    R0, and B0. To drop variation to insignificant levels set sigmar-1-08
#'
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
#' }
oneyearrec <- function(steep,R0,B0,Bsp,sigR) {
  epsilon <- exp(rnorm(length(Bsp),mean=0,sd=sigR) - (sigR * sigR)/2)
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
#' const <- condDat$constants
#' blkpop <- condDat$blkpop
#' numpop <- condDat$numpop
#' blockI <- defineBlock(nblock,blkpop,numpop)
#' popdefs <- definepops(nblock,blockI,const)
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

#' @title restart transfers final year values of zoneD into the first year
#'
#' @description restart transfers the final year values from the
#'     dynamics part of the zone (zoneD), into the first year.
#'     This is used, for example, when searching for an equilibrium
#'     state if there is larval dispersal > 0.0. Of if one sets the
#'     initial depletion to anything other than 1.0. Contains the
#'     option of setting every other cell to zero, which is the
#'     default.
#'
#' @param oldzoneD the old zoneD containing the dynamics as run for
#'     Nyrs.
#' @param hyrs The number of years of dynamics, the hyrs from glb
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
restart <- function(oldzoneD,hyrs,npop,N,zero=TRUE) { # oldzoneD=zoneD; hyrs=hyrs; npop=npop; N=Nclass
  ExplB <- matrix(0,nrow=hyrs,ncol=npop)
  midyexpB <- matrix(0,nrow=hyrs,ncol=npop)
  MatB <- matrix(0,nrow=hyrs,ncol=npop)
  Catch <- matrix(0,nrow=hyrs,ncol=npop)
  Harvest <- matrix(0,nrow=hyrs,ncol=npop)
  cpue <- matrix(0,nrow=hyrs,ncol=npop)
  deplExB <- matrix(0,nrow=hyrs,ncol=npop)
  deplSpB <- matrix(0,nrow=hyrs,ncol=npop)
  Recruit <- matrix(0,nrow=hyrs,ncol=npop)
  CatchN <- array(data=0,dim=c(N,hyrs,npop))
  Nt <- array(data=0,dim=c(N,hyrs,npop))
  zoneD <- list(SAU=oldzoneD$SAU,matureB=MatB,exploitB=ExplB,midyexpB=midyexpB,
                catch=Catch,harvestR=Harvest,cpue=cpue,recruit=Recruit,
                deplsB=deplSpB,depleB=deplExB,catchN=CatchN,Nt=Nt)
  if (!zero) zoneD <- oldzoneD
  zoneD$matureB[1,] <- oldzoneD$matureB[hyrs,]
  zoneD$exploitB[1,] <- oldzoneD$exploitB[hyrs,]
  zoneD$midyexpB[1,] <- oldzoneD$midyexpB[hyrs,]
  zoneD$catch[1,] <- oldzoneD$catch[hyrs,]
  zoneD$harvestR[1,] <- oldzoneD$harvestR[hyrs,]
  zoneD$cpue[1,] <- oldzoneD$cpue[hyrs,]
  zoneD$recruit[1,] <- oldzoneD$recruit[hyrs,]
  zoneD$deplsB[1,] <- oldzoneD$deplsB[hyrs,]
  zoneD$depleB[1,] <- oldzoneD$depleB[hyrs,]
  zoneD$catchN[,1,] <- oldzoneD$catchN[,hyrs,]
  zoneD$Nt[,1,] <- oldzoneD$Nt[,hyrs,]
  return(zoneD)
} # end of restart

#' @title runthreeH conducts the dynamics with constant catch 3 times
#'
#' @description runthreeH is used when searching numerically for an
#'     equilibrium and it conducts the hyrs dynamics three times, each
#'     time through it replaces year 1 with year hyrs. Thus if hyrs is
#'     40 it conducts 3 * 39 years of dynamics (117 years). This is
#'     not exported. It uses zoneC but always it does this inside
#'     the environment of another function where zoneC can be found
#'     Used inside dodepletion and doproduction
#'
#' @param zoneC the constants components of the simulated zone
#' @param zoneD the dynamics portion of the zone, with matrices and
#'     arrays for the dynamic variables of the dynamics of the
#'     operating model
#' @param inHarv a vector, length numpop, of annual harvest rates to be held
#'     constant across all years.
#' @param glob the globals variable from readzonefile
#' @param maxiter default=3; the number of runs through the equilibrium loop.
#'
#' @return a list containing a revised dynamics list, zoneD
#' @export
#'
#' @examples
#' print("wait on built in data sets")
#' # zoneC=zoneC; zoneD=zoneD; glob=glb; inHarv=rep(0.2,numpop)
runthreeH <- function(zoneC,zoneD,inHarv,glob,maxiter=3) {
  npop <- glob$numpop
  Nclass <- glob$Nclass
  hyrs <- glob$hyrs
  larvdisp <- glob$larvdisp
  for (iter in 1:maxiter) { # iter=1; yr=2
    for (yr in 2:hyrs)
      zoneD <- oneyearD(zoneC=zoneC,zoneD=zoneD,Ncl=Nclass,
                        inHt=inHarv,year=yr,sigmar=1e-08,npop=npop,
                        movem=glob$move)
    zoneD <- restart(oldzoneD=zoneD,hyrs=hyrs,npop=npop,N=Nclass,zero=TRUE)
  }
  return(zoneD)
} # end of runthreeH

