% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inputfiles.R
\name{datafiletemplate}
\alias{datafiletemplate}
\title{datafiletemplate generates a template input datafile akin to M15h75}
\usage{
datafiletemplate(indir, filename = "saudataEG.csv")
}
\arguments{
\item{indir}{the directory into which to place the generated data file}

\item{filename}{the name for the generated datafile, a character
string, defaults to saudata_test.csv, which is the default within
the ctrlfiletemplate function}
}
\value{
a standard definition data file, to be read by readdatafile
    whose name and path is returned invisibly
}
\description{
datafiletemplate generates a standard input datafile
    to use as a template, go in and edit it appropriately to suit
    your own needs. It contains the probability distributions that
    are sampled to provide the necessary biological constants for
    each population. It also contains the proportional distribution of
    recruitment levels determined for Tasmania by examining the GPS
    data-logger data for the zone and allocating recruitment in
    proportion to the relative yield by area over the last 8 years.
}
\examples{
\dontrun{
 yourdir <- tempdir()
 datafiletemplate(yourdir,"saudata_test.csv")
 constants <- readsaudatafile(yourdir,"saudata_test.csv")
 str(constants,max.level=1)
 print(constants[,1:10])
}
}
\seealso{
{
  \link{ctrlfiletemplate}, \link{readsaudatafile}, \link{readctrlfile}
}
}
