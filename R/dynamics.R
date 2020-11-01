
#' @title addrepvar adds variation to the start of a generic simulation run
#'
#' @description addrepvar conducts projyrs of simulation at the same constant
#'     harvest rate used to deplete each population to the desired level, only
#'     it does this with recruitment variability turned on, which sets up a
#'     series of replicates with different starting positions.
#'
#' @param zoneC the constant zone object. Its selectivity patterns should be
#'     altered to reflect the projection conditions using modzoneCSel
#' @param zoneDR the dynamic zone object expanded to include all replicates,
#'     using makezoneDR
#' @param inHt the initial harvest strategy that was used to deplete the zone
#' @param glb the globals object
#' @param ctrl the ctrl object
#' @param nyrs defaults to 20, which appears to be enough to produce the
#'     required variation
#'
#' @return a revised zoneDR object
#' @export
#'
#' @examples
#' print("wait on suitable data")
addrepvar <- function(zoneC,zoneDR,inHt,glb,ctrl,nyrs=20) {
  # zoneC=zoneC;zoneDR=zoneDR;inHt=zoneDR$harvestR;glb=glb;ctrl=ctrl
  varzoneDR <- zoneDR
  sigmar <- ctrl$withsigR
  npop <- glb$numpop
  Ncl <- glb$Nclass
  movem <- glb$move
  reps <- ctrl$reps
  matb <- numeric(npop)
  for (iter in 1:reps) {
    for (year in 2:nyrs) {
      for (popn in 1:npop) { # year=2; iter=1; popn=1
        out <- oneyear(inpopC=zoneC[[popn]],inNt=zoneDR$Nt[,year-1,popn,iter],
                       Nclass=Ncl,inH=inHt[1,popn,iter],yr=year)
        zoneDR$exploitB[year,popn,iter] <- out$ExploitB
        zoneDR$matureB[year,popn,iter] <- out$MatureB
        zoneDR$catch[year,popn,iter] <- out$Catch
        zoneDR$harvestR[year,popn,iter] <- out$Harvest
        zoneDR$cpue[year,popn,iter] <- out$ce
        zoneDR$Nt[,year,popn,iter] <- out$Nt
        zoneDR$catchN[,year,popn,iter] <- out$CatchN
        matb[popn] <- out$MatureB
      } # pop
      steep <- getvect(zoneC,"steeph") #sapply(zoneC,"[[","popdef")["steeph",]
      r0 <- getvar(zoneC,"R0") #sapply(zoneC,"[[","R0")
      b0 <- getvar(zoneC,"B0") #sapply(zoneC,"[[","B0")
      recs <- oneyearrec(steep,r0,b0,matb,sigR=sigmar)
      newrecs <- movem %*% recs
      zoneDR$recruit[year,,iter] <- newrecs
      zoneDR$Nt[1,year,,iter] <- newrecs
      zoneDR$deplsB[year,,iter] <- zoneDR$matureB[year,,iter]/b0
      zoneDR$depleB[year,,iter] <- zoneDR$exploitB[year,,iter]/getvar(zoneC,"ExB0")
    }   # year loop        zoneDR$matureB[,,1]
  }     # rep loop
  varzoneDR$matureB[1,,] <- zoneDR$matureB[nyrs,,]
  varzoneDR$exploitB[1,,] <- zoneDR$exploitB[nyrs,,]
  varzoneDR$catch[1,,] <- zoneDR$catch[nyrs,,]
  varzoneDR$recruit[1,,] <- zoneDR$recruit[nyrs,,]
  varzoneDR$cpue[1,,] <- zoneDR$cpue[nyrs,,]
  varzoneDR$catchN[,1,,] <- zoneDR$catchN[,nyrs,,]
  varzoneDR$Nt[,1,,] <- zoneDR$Nt[,nyrs,,]
  return(varzoneDR)
} # end of addrepvar

#' @title depleteSAU resets zoneD to an input depletion level
#'
#' @description depleteSAU resets the depletion level of the whole
#'     zone and does this by searching for the harvest rate that
#'     leads to the sum of the mature biomass, across populations,
#'     divided by the sum of the effective B0 across populations is
#'     as close as possible to the desired depletion level. This means
#'     the individual populations will likely vary around the target
#'     depletion, but across the zone it will be correct. The
#'     depletion is measured relative to the effective B0 as that
#'     takes account of any larval dispersal. This function uses the
#'     production curve matrix to search for harvest rates that bound
#'     the target depletion and then re-searches across those bounds
#'     using len intervals.
#'
#' @param zoneC the constants components of the simulated zone
#' @param zoneD the dynamic components of the simulated zone
#' @param glob the general global variables
#' @param depl a vector of target depletion levels for each SAU
#' @param product the production curve matrix from doproduction
#' @param len the number of intervals in the trial harvest rates
#'
#' @return a revised zoneD object
#' @export
#'
#' @examples
#' \dontrun{
#'   data(zone1)
#'   glb <- zone1$globals
#'   data(constants)
#'   ans <- makezoneC(zone1,constants)
#'   zoneC <- ans$zoneC
#'   glb <- ans$glb  # now contains the move matrix
#'   ans <- makezone(glb,zoneC)
#'   zoneC <- ans$zoneC
#'   zoneD <- ans$zoneD
#'   ans2 <- modzoneC(zoneC,zoneD,glb)
#'   zoneC <- ans2$zoneC  # now has MSY and deplMSY
#'   product <- ans2$product
#'   zoneDD <- dodepletion(zoneC,zoneD,glb,depl=0.3,product)
#'   sum((zoneDD$matureB[1,]/sum(zoneDD$matureB[1,]))*zoneDD$deplsB[1,])
#'   mean(zoneDD$deplsB[1,])
#' }
depleteSAU <- function(zoneC,zoneD,glob,depl,product,len=15) {
  #  zoneC=zoneC; zoneD=zoneD; glob=glb;  product=product; depl=rep(0.8,8);len=15
  # use product to find bounds on H around the desired depletion level
  nSAU <- glob$nSAU
  harvests <- matrix(0,nrow=len,ncol=nSAU,dimnames=list(1:len,zoneD$SAU))
  sauDepl <- harvests # need a matrix to hold the depletion levels
  if (length(depl) != nSAU) stop("Need a depletion for each SAU \n")
  initH <- as.numeric(rownames(product))
  for (sau in 1:nSAU) {
    pick <- which.closest(depl[sau],product[,"Deplet",sau])
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
    harvests[,sau] <- seq(lowl,upl,length=len)
  }
  for (harv in 1:len) { # harv=1
    zoneDD <- runthreeH(zoneC=zoneC,zoneD,inHarv=harvests[harv,],glob)
    sauDepl[harv,] <- zoneDD$deplsB[1,]

  }
  pickharv <- numeric(nSAU)
  for (sau in 1:nSAU)
    pickharv[sau] <- which.closest(depl[sau],harvests[,sau],index=FALSE)
  zoneDD <- runthreeH(zoneC=zoneC,zoneD,inHarv=pickharv,glob)
  return(zoneDD)
} # end of depleteSAU

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
#' @param sel the selectivity vector or matrix to be used
#' @param selwt the selectivity x WtL vector or matrix to be used
#' @param pop the index number of the population
#'
#' @return a list containing ExploitB, MatureB, Catch, Harvest, Nt,
#'     ce, and CatchN used to update the given pop in yr + 1
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
  oldExpB <- ExploitB   # ExploitB after growth and 0.5NatM
  Fish <- 1-(inH*selyr)
  newNt <- (Os * (Fish * NumNe)) #+ Rec # Nt - catch - 0.5M, and + Rec
  Cat <- (inH*selyr) * NumNe  #numbers at size in the catch
  newExpB <- sum(SelectWt * newNt)
  avExpB <- (newExpB + oldExpB)/2.0 #av start and end
  MatureB <- sum(MatWt*newNt) #+ MatBC
  Catch <- sum(inpopC$WtL*Cat)/1e06
 # Harvest <- Catch/avExpB  # uses average of the start and end
  ce <- inpopC$popq * avExpB * 1000.0  #ExploitB
  ans <- list(avExpB,MatureB,Catch,inH,newNt,ce,Cat) # use inH not Harvest
  names(ans) <- c("ExploitB","MatureB","Catch","Harvest","Nt","ce",
                  "CatchN")
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
#' @return a list containing ExploitB, MatureB, Catch, Harvest, Nt,
#'     ce, and CatchN used to update the given pop in yr + 1
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
  #  MatureB <- sum(MatWt*inNt)
  NumNe <- (Os * (inpopC$G %*% inNt))
  ExploitB <- sum(SelectWt * NumNe) #SelectWt=Select*WtL
  oldExpB <- ExploitB   # ExploitB after growth and 0.5NatM
  estH <- min(incat/ExploitB,0.8) # no more than 0.8 harvest rate
  Fish <- 1-(estH*selyr)
  newNt <- (Os * (Fish * NumNe)) #+ Rec # Nt - catch - 0.5M, and + Rec
  Cat <- (estH*selyr) * NumNe  #numbers at size in the catch
  ExploitB <- (sum(SelectWt * newNt) + oldExpB)/2.0 #av start and end
  MatureB <- sum(MatWt*newNt) #+ MatBC
  Catch <- sum(inpopC$WtL*Cat)/1e06
  Harvest <- min(Catch/ExploitB,0.8)  # average of the start and end
  ce <- inpopC$popq * ExploitB * 1000.0  #ExploitB
  ans <- list(ExploitB,MatureB,Catch,Harvest,newNt,ce,Cat)
  names(ans) <- c("ExploitB","MatureB","Catch","Harvest","Nt","ce",
                  "CatchN")
  return(ans)
} # End of oneyearcat

#' @title oneyearC conducts one year's dynamics using catch not harvest
#'
#' @description oneyearD conducts one year's dynamics in the simulation
#'     using catches rather than harvest rates. The harvest rates are
#'     estimated after first estimating the exploitable biomass.
#'     returning the revised zoneD, which will have had a single year
#'     of activity included in each of its components.
#'
#' @param zoneC the constant portion of the zone with a list of
#'     properties for each population
#' @param zoneD the dynamics portion of the zone, with matrices and
#'     arrays for the dynamic variables of the dynamics of the
#'     operating model
#' @param catchp a vector of catches to be taken in the year from each
#'     population
#' @param year the year of the dynamics, would start in year 2 as year
#'     1 is the year of initiation.
#' @param sigmar the variation in recruitment dynamics, set to 1e-08
#'     when searching for equilibria.
#' @param Ncl the number of size classes used to describe size, global Nclass
#' @param npop the number of populations, the global numpop
#' @param movem the larval dispersal movement matrix, global move
#'
#' @return a list containing a revised dynamics list
#' @export
#'
#' @examples
#'  data(constants)
#'  data(zone1)
#'  ans <- makezoneC(zone1,constants) # classical equilibrium
#'  zoneC <- ans$zoneC
#'  glb <- ans$glb
#'  ans <- makezone(glb,zoneC) # now make zoneD
#'  zoneC <- ans$zoneC  # zone constants
#'  zoneD <- ans$zoneD
#'  Nc <- glb$Nclass
#'  nyrs <- glb$Nyrs
#'  catch <- 360.0 # larger than total MSY ~ 310t
#'  B0 <- getvar(zoneC,"B0")
#'  totB0 <- sum(B0)
#'  prop <- B0/totB0
#'  catchpop <- catch * prop
#'  for (yr in 2:nyrs)
#'    zoneD <- oneyearC(zoneC=zoneC,zoneD=zoneD,Ncl=Nc,
#'                    catchp=catchpop,year=yr,sigmar=1e-08,
#'                    npop=glb$numpop,movem=glb$move)
#'  str(zoneD)
#'  round(zoneD$catchN[60:105,1:5,1],1)
oneyearC <- function(zoneC,zoneD,catchp,year,sigmar,Ncl,npop,movem) {
  matb <- numeric(npop)
  for (popn in 1:npop) {  # year=2
    out <- oneyearcat(inpopC=zoneC[[popn]],inNt=zoneD$Nt[,year-1,popn],
                      Nclass=Ncl,incat=catchp[popn],yr=year)
    zoneD$exploitB[year,popn] <- out$ExploitB
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
} # end of oneyearC


#' @title oneyearD conducts one year's dynamics on zoneD in the MSE
#'
#' @description oneyearD conducts one year's dynamics on zoneD in the MSE
#'     returning the revised zoneD, which will have had a single year
#'     of activity included in each of its components. This uses zoneC
#'     but always within the environment of another function in which
#'     zoneC (as zoneC) can be found. Used in runthreeH, (and hence
#'     dodepletion and doproduction) and in testequil
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
#'
#' @examples
#'  data(constants)
#'  data(zone1)
#'  ans <- makezoneC(zone1,constants) # classical equilibrium
#'  zoneC <- ans$zoneC
#'  glb <- ans$glb
#'  ans <- makezone(glb,zoneC) # now make zoneD
#'  zoneC <- ans$zoneC  # zone constants used as zoneC in oneyearD
#'  zoneD <- ans$zoneD
#'  numpop <- glb$numpop
#'  harvest <- rep(0.2,numpop) # not exported so needs aMSE:::
#'  zoneD <- aMSE:::oneyearD(zoneC=zoneC,zoneD=zoneD,
#'                    inHt=harvest,year=2,sigmar=1e-06,
#'                    Ncl=glb$Nclass,npop=numpop,movem=glb$move)
#'  str(zoneD)
#'  round(zoneD$catchN[60:105,1:5,1],1)
#'  zoneC=zoneC;zoneD=zoneD;Ncl=Nclass;inHt=inHarv;year=yr;sigmar=1e-08;npop=npop;movem=glob$move
oneyearD <- function(zoneC,zoneD,inHt,year,sigmar,Ncl,npop,movem) {
  matb <- numeric(npop)
  for (popn in 1:npop) {  # year=2; popn=1
    out <- oneyear(inpopC=zoneC[[popn]],inNt=zoneD$Nt[,year-1,popn],
                   Nclass=Ncl,inH=inHt[popn],yr=year)
    zoneD$exploitB[year,popn] <- out$ExploitB
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
#' oneyearrec(steep,R0,B0,Bsp,insigmar)
#' }   # steep=steep; R0=mover0; B0=b0; Bsp=b0
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
restart <- function(oldzoneD,nyrs,npop,N,zero=TRUE) { # oldzoneD=zoneD; nyrs=Nyrs; npop=npop; N=Nclass
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
  zoneD <- list(SAU=oldzoneD$SAU,matureB=MatB,exploitB=ExplB,
                catch=Catch,harvestR=Harvest,cpue=cpue,recruit=Recruit,
                deplsB=deplSpB,depleB=deplExB,catchN=CatchN,Nt=Nt)
  if (!zero) zoneD <- oldzoneD
  zoneD$matureB[1,] <- oldzoneD$matureB[nyrs,]
  zoneD$exploitB[1,] <- oldzoneD$exploitB[nyrs,]
  zoneD$catch[1,] <- oldzoneD$catch[nyrs,]
  zoneD$harvestR[1,] <- oldzoneD$harvestR[nyrs,]
  zoneD$cpue[1,] <- oldzoneD$cpue[nyrs,]
  zoneD$recruit[1,] <- oldzoneD$recruit[nyrs,]
  zoneD$deplsB[1,] <- oldzoneD$deplsB[nyrs,]
  zoneD$depleB[1,] <- oldzoneD$depleB[nyrs,]
  zoneD$catchN[,1,] <- oldzoneD$catchN[,nyrs,]
  zoneD$Nt[,1,] <- oldzoneD$Nt[,nyrs,]
  return(zoneD)
} # end of restart

#' @title runthree conducts the dynamics with constant catch 3 times
#'
#' @description runthree is used when searching numerically for an
#'     equilibrium and it conducts the Nyrs dynamics three times, each
#'     time through it replaces year 1 with year Nyrs. Thus if Nyrs is
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
#'
#' @examples
#' print("wait on built in data sets")
#' # zoneC=zoneC; zoneD=zoneD; glob=glb; inHarv=rep(0.2,numpop)
runthreeH <- function(zoneC,zoneD,inHarv,glob,maxiter=3) {
  npop <- glob$numpop
  Nclass <- glob$Nclass
  Nyrs <- glob$Nyrs
  larvdisp <- glob$larvdisp
  for (iter in 1:maxiter) { # iter=1; yr=2
    for (yr in 2:Nyrs)
      zoneD <- oneyearD(zoneC=zoneC,zoneD=zoneD,Ncl=Nclass,
                        inHt=inHarv,year=yr,sigmar=1e-08,npop=npop,
                        movem=glob$move)
    zoneD <- restart(oldzoneD=zoneD,nyrs=Nyrs,npop=npop,N=Nclass,zero=TRUE)
  }
  return(zoneD)
} # end of runthree

