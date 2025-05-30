% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getfunctions.R
\name{getLogical}
\alias{getLogical}
\title{getLogical extracts nb logicals from an input line of text}
\usage{
getLogical(inline, nb)
}
\arguments{
\item{inline}{text line to be parsed, usually obtained using readLines}

\item{nb}{the number of logicals to extract, if nb is longer than the
number of logicals within inline the vector will contain NAs}
}
\value{
a vector of length nb
}
\description{
getLogical obtains nb logicals from an input line
}
\examples{
 txtline <- "Depleted, TRUE"
 getLogical(txtline,nb=1)
 txtline2 <- "calcthis, TRUE, FALSE"
 getLogical(txtline2,nb=2)
}
