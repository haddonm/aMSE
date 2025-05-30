% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comparescenarios.R
\name{boxbysau}
\alias{boxbysau}
\title{boxbysau generates boxplots of sau x aavc or sum5 or sum10}
\usage{
boxbysau(
  rundir,
  hspm,
  glbc,
  scenes,
  compvar = "aavc",
  filen = "",
  aavcyrs = 10,
  maxval = 0
)
}
\arguments{
\item{rundir}{the directory in which comparisons are being made. It is best
to be a separate directory form any particular scenario.}

\item{hspm}{a list of the harvest strategy performance measures for actual
catches as calculated by calccatchHSPM}

\item{glbc}{a list of the globals objects from each scenario}

\item{scenes}{the names of the different scenarios being compared}

\item{compvar}{which variable should be compared. The options are 'aavc',
'sum5', and 'sum10'. default = 'aavc'}

\item{filen}{the filename for saving the plot, the default="", which sends
the plot to the console.}

\item{aavcyrs}{the number of projection years to use for the aavc calculation
the default = 10}

\item{maxval}{the maximum y-axis value for the aavc, default=0, which means
the maximum for each boxplot is taken from the data for that plot.}
}
\value{
a list of the boxplot statistics for the aavc x number of sau x the
    number of scenarios compared
}
\description{
boxbysau is used to generate boxplots of the aavc, or sum5,
    of sum10 for each sau, in each of the scenarios. The number of years over
    which the aavc is calculated can be modified using the aavcyrs argument.
    If the replicates differ between scenarios then these boxplots will be
    omitted.
}
\examples{
print("wait on data sets")
}
\seealso{
{
   \link{getprojyrs}, \link{getprojyraavc},  \link{catchHSPM}
}
}
