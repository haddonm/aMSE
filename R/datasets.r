
#' @title blockE13 is summary fishery data for East coast block 13
#'
#' @description blockE13 is summary fishery data for East coast block
#'     13 blacklip abalone (\emph{Haliotis rubra}), as might be used
#'     to condition the operating model. It can be used for exploring
#'     the performance measure estimators getgrad1 and getscore.
#'
#' @name blockE13
#'
#' @docType data
#'
#' @section contents:
#' \itemize{
#'   \item year the calander year for the fishery data
#'   \item Coeff is the back-transformed standardized coefficient for
#'       the year, appropriately bias-corrected for the Log-Normal
#'       errors.
#'   \item SE is the standard error of the coefficient
#'   \item Scaled is the Coeff scaled to the mean of the time-series
#'   \item cpue is the Scaled values multiplied by the bias-corrected
#'       geometric mean of the entire data-set
#'   \item catch is the total catch from the original data
#'   \item effort is the total effort as hours from the original data
#' }
#'
#' @examples
#'  data(blockE13)
#'  print(blockE13)
NULL

#' @title constants is conditioning data for 16 populations and 8 SAU
#'
#' @description constants is a data.frame of parameters for blacklip
#'     abalone (\emph{Haliotis rubra}) used to condition the
#'     operating model as an example when running the aMSE function
#'     examples. It describes a zone containing eight hypothetical
#'     Tasmanian blocks (the level of SAU) with a total of sixteen
#'     populations (two per SAU). If this were a data-set in a CSV file it
#'     would require the readdatafile function to read the file. The
#'     function datafiletemplate will generate a template of the required
#'     format to be read by readdatafile, which can be edited to suit a
#'     given fishery. once transposed and with random variation included,
#'     this data becomes the contents of popdefs.
#'
#' @name constants
#'
#' @docType data
#'
#' @section contents:
#' \itemize{
#'   \item popnum the index to the population
#'   \item SAU the spatial management unit number = blocks in Tasmania
#'   \item DLMax maximum growth increment used in the inverse logistic
#'   \item sMAxDL the sd variability of DLMax, Normal variation
#'   \item L50  the shell length at 50% of maximum growth increment
#'   \item sL50 the sd variability of L50, Normal variation
#' }
#'
#' @examples
#'  data(constants)
#'  print(constants)
NULL


#' @title ctrl the control file for a particular MSE run
#'
#' @description crtl contains the information required to conduct a
#'     particular MSE run. It identifies directories, filenames, the
#'     number of replicates,m and other variables
#'
#' @name ctrl
#'
#' @docType data
#'
#' @section contents:
#' \itemize{
#'   \item runlabel  the identifying name for the run
#'   \item zonefile  filename containing the zone data, see zone1
#'   \item datafile  filename containing the population defintions, see constants
#'   \item hcrfile  filename continaing the details of the HCR used
#'   \item outdir  the output directory, containing a results subdir
#'   \item reps  how many replicates in this instance, usually 1000
#'   \item initdepl  the initial depletion level for the zone
#'   \item assessinterval  how often should the zone be assessed
#'   \item withsigR  the level of recruitment variability
#'   \item withsigB  the level of noise on biomass estimates
#'   \item withsigce  the level of noise on cpue estimates
#' }
#'
#' @examples
#'  data(ctrl)
#'  print(ctrl)
NULL

#' @title midg is an abalone tagging data-set from the Actaeons
#'
#' @description midg is a tagging data-set for blacklip abalone
#'     (\emph{Haliotis rubra}) from the middle ground in the Actaeons
#'     in Tasmania's Block 13. All individuals were recaptured in 2003,
#'     the site number was 478, at Latitude -43.54 longitude 146.99,
#'     and there are 347 observations. All Dt = 1 year.
#'
#' @name midg
#'
#' @docType data
#'
#' @format A data.frame of abalone tagging data
#' \describe{
#'   \item{RecapL}{the length at recapture}
#'   \item{Lt}{the length at tagging}
#'   \item{Dt}{The time interval between tagging and recapture, in this
#'       instance they are all listed as 1 year}
#'   \item{DL}{the growth increment in mm}
#' }
#'
#' @section Subjects:
#'  \itemize{
#'    \item growth curves
#'    \item inverse logistic, von Bertalanffy, Gompertz
#'    \item Static model fitting
#'  }
#'  @export
#'
#' @source Thanks to the Institute of Marine and Antarctic Science,
#'     which is part of the University of Tasmania, and especially to
#'     Dr Craig Mundy, leader of the Abalone Group, for permission to use
#'     this data collected in 2003.
#'
#' @examples
#'  data(midg)
#'  head(midg,20)
#'  oldpar <- par(no.readonly=TRUE)
#'  plot(midg$Lt,midg$DL,type="p",pch=16,cex=1.0,xlim=c(5,180))
#'  abline(h=0,col=1)
#'  par(oldpar)
NULL

#' @title zone1 the constants common to a zone
#'
#' @description zone1 contains the constants relating to the whole
#'     zone rather than the populations. See the listing below.
#'
#' @name zone1
#'
#' @docType data
#'
#' @format A list of constants that are uniform across all populations
#'     in a zone
#' \describe{
#'   \item{SAUnames}{the names of each SAU}
#'   \item{SAUpop}{a vector of how many populations in each SAU}
#'   \item{minc}{the mid-point of the minimum size class}
#'   \item{cw}{the size-class width in mm}
#'   \item{larvdisp}{the rate of larval dispersal as a proportion}
#'   \item{randomseed}{used if results need repeating}
#'   \item{initLML}{the LML used to generate teh equilibrium population}
#'   \item{condC}{a list containing any historical data used to condition the
#'       initial zone}
#'   \item{projC}{a list containing the details required to prepare the zoneC
#'       and zoneD for the projections}
#'   \item{globals}{a list of global constants, containing numpop,
#'       nSAU, midpts, Nclass, Nyrs, and larvdisp}
#'   \item{ctrl}{the list containing control information for the run, including
#'       the datafile for the constants, the reps, the variation to be included
#'       when projecting}
#'   \item{condition}{if zero then no historical conditioning, if >0, then the
#'       number of years of conditioning catch data, complements condC}
#'   \item{projection}{if zero then no projection conducted, if >0 then the
#'       number of years of projection}
#' }
#'
#' @examples
#'  data(zone1)
#'  str(zone1,max.level=1)
NULL

#' @title tasab is a matrix of abalone maturity-at-length data
#'
#' @description tasab is a 715 x 4 matrix of maturity-at-length data
#'     for blacklip abalone (\emph{Haliotis rubra}) from two sites
#'     along the Tasmanian west coast. All data was collected in
#'     February 1995, but details, such as site name, accurate
#'     location, statistical block, year, month, and other
#'     details have been omitted for brevity.
#'
#' @name tasab
#'
#' @docType data
#'
#' @format A data.frame of maturity-at-length data
#' \describe{
#'   \item{site}{an identifier for the two different sites sampled}
#'   \item{sex}{I = immature, M = male, F = female}
#'   \item{length}{the shell length in mm}
#'   \item{mature}{was the animal mature = 1 or not = 0}
#' }
#'
#' @section Subjects:
#'  \itemize{
#'    \item maturity ogives or logistic curves
#'    \item Binomial likelihoods
#'  }
#'  @export
#'
#' @source Thanks to the Institute of Marine and Antarctic Science,
#'     which is part of the University of Tasmania, and especially to
#'     Dr Craig Mundy, leader of the Abalone Group, for permission to use
#'     this data collected in February 1995.
#'
#' @examples
#'  data(tasab)
#'  head(tasab,20)
#'  table(tasab$site,tasab$sex)
NULL


