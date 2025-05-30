% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getfunctions.R
\name{getvar}
\alias{getvar}
\title{getvar a replacement for sapply to obtain scalar constants}
\usage{
getvar(zoneC, invar)
}
\arguments{
\item{zoneC}{the constants object for the zone}

\item{invar}{a character variable eg. "B0" or "R0"}
}
\value{
a numpop vector of the invar constants from zoneC
}
\description{
getvar is a replacement for sapply to obtain scalar
    constants from zoneC and is significantly faster. It should
    be used to obtain things like B0, R0, MSY, scalece, etc. Still
    need to use sapply to pull out vectors.
}
\examples{
data(zone)
zoneC <- zone$zoneC
getvar(zoneC,"MSY")
getvar(zoneC,"B0")
}
\seealso{
{
 \link{getlistvar}
}
}
