% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comparescenarios.R
\name{plotzonedevs}
\alias{plotzonedevs}
\title{plotzonedevs plots projection deviates for a variable across scenerios}
\usage{
plotzonedevs(invar, scenes, glb, filen = "")
}
\arguments{
\item{invar}{the sub-object from each scenario as a list, see examples}

\item{scenes}{a vector of names fr each scenario being compared}

\item{glb}{a globals object from any scenario. This is to obtain the
number of projection years, which MUST be the same across
all scenarios (else the matrices will not line up).}

\item{filen}{the filename and path if the plot is to be saved, default = "".
if a filename is given it should be a .png file}
}
\value{
a list of a matrix of the mean and sd devs per scenario, and a matrix
     of the actual deviates, all returned invisibly. It also plots a graph
}
\description{
plotzonedevs provides a plot of a scenario Performance Measure.
    In this case the PM are the deviations from a loess fit to whatever input
    variable has been selected from the zone objects in do_comparison(). One
    needs first generate a list of the selected variable from the the input
    zone objects (one from each scenario). A standard comparison would be on
    the variation of catches through the years of the projection. But this
    could be applied to any of the variables within the zone object, which
    includes matureB, exploitB, midyexpB, catch, acatch, harvestR, cpue,
    recruit, deplsB, depleB. The final catchN and Nt would need separate
    treatment.
}
\examples{
\dontrun{
# example syntax, assumes zone has been extracted from all scenarios
   glb <- glbc[[1]]
   nscen <- length(scenes)
   whichyrs <- (glb$hyrs + 1):(glb$hyrs + glb$pyrs)
   x <- makelist(scenes)
   for (scen in 1:nscen) x[[scen]] <- zone[[scen]]$catch[whichyrs,]
   devout <- plotzonedevs(x,scenes,glb,filen="")
}
}
