% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getfunctions.R
\name{getmedbysau}
\alias{getmedbysau}
\title{getmedbysau extracts the median projected values for input variables}
\usage{
getmedbysau(invar, glb)
}
\arguments{
\item{invar}{this must be one of the objects from the sauout list.}

\item{glb}{the globals object}
}
\value{
a matrix of the medians across replciates for each sau of the input
    variable
}
\description{
getmedbysau takes the output from the projections summarized
    to an sau scale, found in sauout. This contains matureB, exploitB,
    midyexpB, catch, acatch, harvestR, cpue, recruit, deplsB, depleB, catchN,
    and Nt, each of which is for (hyrs + pyrs) x nsau x reps, or nclass x
    (hyrs + pyrs) x nsau x reps. The input variable can only be one if these.
}
\examples{
print("wait on data sets")
# typical syntax woul dbe getmedbysau(out$sauout$matureB,out$glb)
}
