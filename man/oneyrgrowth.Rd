% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamics.R
\name{oneyrgrowth}
\alias{oneyrgrowth}
\title{oneyrgrowth one years' growth for a population and initial size}
\usage{
oneyrgrowth(inpop, startsize = 2)
}
\arguments{
\item{inpop}{the abpop to be grown forward}

\item{startsize}{default = 2, but often set to the size at 50 percent
maturity or the vector of midpts to determine the growth of Nt}
}
\value{
the expected mean size after one year's growth
}
\description{
oneyrgrowth one years' growth for a given population and
   initial size. Used to determine the size after two year's growth as
   a measure of the expected bLML - biological Legal Minimum Length. A
   reflection of the two year rule in Tasmania. To get this the function
   should obviously be run twice.
}
\examples{
\dontrun{
nblock <- 2
filename <- datafiletemplate(numblock=nblock,filename="block3.csv")
condDat <- readdatafile(filename)
const <- condDat$constants
blkpop <- condDat$blkpop
numpop <- condDat$numpop
popdefs <- definepops(nblock,blockI,const)
pop1 <- makeabpop(popdefs[1,],condDat$midpts,
                  condDat$ProjLML[popdefs[1,20]])
oneyrgrowth(pop1,startsize=70)
}
}
