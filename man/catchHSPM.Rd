% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comparescenarios.R
\name{catchHSPM}
\alias{catchHSPM}
\title{catchHSPM generates boxplots of catch HSPMs and output statistics}
\usage{
catchHSPM(
  rundir,
  hspm,
  glbc,
  scenes,
  filen = "",
  aavcyrs = 10,
  maxvals = c(0.25, 750, 1500),
  wid = 12
)
}
\arguments{
\item{rundir}{the directory in which comparisons are being made. It is best
to be a separate directory form any particular scenario.}

\item{hspm}{a list of the harvest strategy performance measures for actual
catches as calculated by calccatchHSPM}

\item{glbc}{a list of the globals objects from each scenario}

\item{scenes}{the names of the different scenarios being compared}

\item{filen}{the filename for saving the plot, the default="", which sends
the plot to the console.}

\item{aavcyrs}{the number of projection years to use for the aavc calculation
the default = 10}

\item{maxvals}{the maximum y-axis value for the aavc, the sum5, and the sum10,
default=c(0.25,750,1500)}

\item{wid}{the width of hteplot, default = 12, room for 3 scenarios}
}
\value{
a list of the boxplot statistics for the three HSPMs x the number of
    scenarios compared
}
\description{
catchHSPM is used to generate boxplots of the aavc, the sum5,
    and the sum10 for projected catches. The number of years over which the
    aavc is calculated can be modified using the aavcyrs argument. In
    addition to the boxplots, the function outputs the statistics describing
    the statistics from the boxplots.
}
\examples{
print("wait on data sets")
# rundir=rundir;hspm=outcatchHSPM;glbc=glbc;scenes=scenes
# filen="compare_catches_boxplots.png";aavcyrs=10
}
\seealso{
{
   \link{getprojyrs}, \link{getprojyraavc},  \link{getprojyrC}
}
}
