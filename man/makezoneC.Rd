% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defineZone.R
\name{makezoneC}
\alias{makezoneC}
\title{makezoneC makes the constant part of the simulation}
\usage{
makezoneC(zone, const)
}
\arguments{
\item{zone}{the object derived from the readzonefile function}

\item{const}{the object derived from the readdatafile function}
}
\value{
a list containing the constant part of the simulated zone
}
\description{
makezoneC makes the constant part of the simulated
    zone. Once defined this does not change throughout the
    simulation. Once made it still requires makezone to be run
    to fill in the B0, ExBo, MSY, MSYDepl, and the scalece values, and
    to produce zoneD, the dynamic part of the new zone
}
\examples{
\dontrun{
data(zone1)
data(constants)
ans <- makezoneC(zone=zone1,const=constants)
zoneC <- ans$zoneC
popdefs <- ans$popdefs
str(zoneC,max.level=1)
str(zoneC[[1]])  # not complete at this stage
print(popdefs)
}
}
