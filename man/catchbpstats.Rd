% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comparescenarios.R
\name{catchbpstats}
\alias{catchbpstats}
\title{catchbpstats saves the boxplot statistics and adds to MSE webpage}
\usage{
catchbpstats(rundir, outtab)
}
\arguments{
\item{rundir}{the directory into which all the comparison results are placed}

\item{outtab}{the output from catchHSPM}
}
\value{
invisibly returns the matrix version of the input list. It also
    saves a csv file and adds to the webpage
}
\description{
catchbpstats takes the results from plotting the boxplots of the
    catch HSPMs aavc, sum5, and sum10 and tabulates the statistical
    properties of the boxplots, ie, the lower and upper whiskers, the
    inter-quartile bounds and the median. These are all output into a single
    table. The aavc, sum5, and sum10 results for each scenario are grouped
    together in the final data.frame.
}
\examples{
print("wait on example data sets")
}
\seealso{
\link{catchHSPM}
}
