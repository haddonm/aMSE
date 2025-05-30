% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotfuns.r
\name{dosau}
\alias{dosau}
\title{dosau plots the conditioning history for the dynamics}
\usage{
dosau(inzone, glb, picksau, histCE, histCatch, yrnames, recdev)
}
\arguments{
\item{inzone}{the conditioned zone (zoneDD) after it has been converted to
sau scale by using getsauzone}

\item{glb}{the global constants object}

\item{picksau}{which sau should be plotted}

\item{histCE}{the historical cpue series from condC}

\item{histCatch}{the historical catch series from condC}

\item{yrnames}{the years of historical catches eg 1973:2019}

\item{recdev}{the recdevs for a single SAU}
}
\value{
nothing but it does generate a plot
}
\description{
dosau plots the deplsB, cpue, matureB, catch, harvestR, and
    recruits for a given SAU during the years of conditioning. This aims to
    assist the conditioning process by illustrating the state of the sau
    during and at the end of conditioning on the fishery. If extra is TRUE
    then it also plots the exploitable biomass and the depletion of the
    exploitable biomass. A loess fit is also plotted onto the recruitments
    plot to give some insight into the recruitment deviates. Unfortunately,
    the 'true' predicted recruitment levels (without variation) cannot be
    easily estimated because the implementation of larval movement disturbs
    how many recruits are present each year in each population and hence
    each SAU.
}
\examples{
print("wait on suitable data-sets")
}
