% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/docondition.R
\name{optimizerecdevs}
\alias{optimizerecdevs}
\title{optimizerecdevs improves the fit to the recdevs for a given sau}
\usage{
optimizerecdevs(
  rundir,
  sau,
  controlfile,
  calcpopC,
  wtsc = 0.02,
  maxiter = 100,
  yearrange = 1980:2016,
  verbose = FALSE,
  mincount = 100,
  plottoconsole = FALSE,
  optimmethod = "Nelder-Mead"
)
}
\arguments{
\item{rundir}{the rundir for the scenario}

\item{sau}{which sau to apply this to}

\item{controlfile}{the name of the control file}

\item{calcpopC}{the funciton, from the HS file that allocates catches
across each of the populations within each SAU}

\item{wtsc}{what weight to give to the size-composition data,default=0.1}

\item{maxiter}{The maximum number of iterations to run before stopping.
Default = 100, one needs at least 400 to see a real difference}

\item{yearrange}{which recdev years should be 'conditioned'. Default=
1980:2016}

\item{verbose}{text output to console, default=FALSE}

\item{mincount}{the minimum number of iterations in the solver, default=100}

\item{plottoconsole}{should the final plot be sent to the console = TRUE
of the rundir = FALSE. default=FALSE}

\item{optimmethod}{which optim method to use? default='Nelder-Mead'}
}
\value{
a scalar value which is the total SSQ for the selected years for
    the CPUE and sizecomps. It also alters the recdevs in the controlfile!
}
\description{
optimizerecdevs an important part of optimizing the match
    between the dynamics observed in the fishery and those of the
    operating model is to adjust the annual recruitment deviates so that
    they improve the fit to the observed CPUE and observed numbers-at-size.
    optimizerecdevs attempts to do that.
}
\examples{
print("wait on data sets")
}
\seealso{
{
    \link{adjustavrec}
}
}
