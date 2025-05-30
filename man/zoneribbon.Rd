% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comparescenarios.R
\name{zoneribbon}
\alias{zoneribbon}
\title{zoneribbon a ribbon plot the quantiles x scenes for the input variable}
\usage{
zoneribbon(
  rundir,
  scenes,
  invar,
  glbc,
  varname,
  category,
  console = TRUE,
  q90 = TRUE,
  intens = 127,
  addleg = "topleft",
  addmedian = 3,
  add1line = TRUE
)
}
\arguments{
\item{rundir}{the directory in which all results are held for a scenario or
comparison of scenarios}

\item{scenes}{the names of the different scenarios being compared}

\item{invar}{a list of the zone-scale replicate runs for the varname from
each scenario. the zone object is generated inside docomparison}

\item{glbc}{the list of the global objects.}

\item{varname}{just the name of the variable being plotted, to make sure the
figures all have the correct labelling}

\item{category}{if a file is saved this argument identifies the webpage tab
into which the plot will be placed.}

\item{console}{should the plot be sent to the console or saved for use in the
web-page output? default=TRUE, set to FALSE to save it}

\item{q90}{should the 90th quantiles be plotted? default=TRUE. If FALSE then
the 95th quantiles are used.}

\item{intens}{if polys=TRUE then intens signifies the intensity of colour on
a scale of 0 - 255. 127 is about 50 percent dense.}

\item{addleg}{add a legend or not. If not use "", if yes use any of topleft,
topright, bottomleft, or bottomright, default = bottomright.}

\item{addmedian}{should each polygon also have its median plotted as a line.}

\item{add1line}{default = TRUE, adds a thicker dashed black line at 1.0,
where the catch/msy ratio = 1.0}
}
\value{
a list containing the quantiles for the varname for each of the
    scenarios - invisibly
}
\description{
zoneribbon is used to produce a ribbon plot for the zone across
    a set of scenarios for a given variable (cpue, catch, harvest rate, etc)
    out of the zone object. After running do_comparison
}
\examples{
\dontrun{
 # a possible run where scenarios 1, 13, and 14 are compared.
   result <- do_comparison(rundir=rundir,postfixdir=postfixdir,outdir=outdir,
                           files=files,pickfiles=c(1,2,3),verbose=TRUE)
}
}
\seealso{
{
  \link{doquantribbon}, \link{do_comp_outputs}
}
}
