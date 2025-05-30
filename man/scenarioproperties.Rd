% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comparescenarios.R
\name{scenarioproperties}
\alias{scenarioproperties}
\title{scenarioproperties tabulates important properties of each scenario}
\usage{
scenarioproperties(scenes, glbc, ctrlc, condCc)
}
\arguments{
\item{scenes}{the runnames of each scenario}

\item{glbc}{the globals object or list for each scenario}

\item{ctrlc}{the ctrl object for each scenario}

\item{condCc}{the condC object for each scenario}
}
\value{
a list containing a matrix of a set of scenario properties to check
    comparability, and a vector identifying whether all columns are identical
    = 1 or some are different = 0
}
\description{
scenarioproperties is used with the do_comparison function to
    tabulate a set of important properties of each scenario so that the
    comparability of each scenario can be quickly assessed.
}
\examples{
print("wait on data sets")
}
\seealso{
{
  \link{do_comparison}
}
}
