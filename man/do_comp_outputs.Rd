% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comparescenarios.R
\name{do_comp_outputs}
\alias{do_comp_outputs}
\title{do_comp_outputs extracts components of the output from do_comparison}
\usage{
do_comp_outputs(result, projonly = TRUE)
}
\arguments{
\item{result}{the output from the do_comparison function}

\item{projonly}{should only the projection years be returned? default=TRUE.}
}
\value{
a list of one vector of names and 7 potentially large lists. Most
    importantly, popout has the dynamics at a population scale, dyn has them
    at an sau scale, and zone has them at the zone scale.
}
\description{
do_comp_outputs is a helper function that extracts a number of
    the more important components from each scenario into lists of their
    own. It does this for 'dyn', which contains the replicate outputs from
    the dynamics summarized at the sau scale (this includes Nt and catchN).
    There is also 'zone', which does the same thing but at a zone scale.
    There is also scenes, which contains the names of each scenario, 'prods',
    which contains the productivity characteristics of each scenario, and
    'scores', which contains the HS scores for each internal performance
    measure.
}
\examples{
print("wait on datasets")
# result=result; projout=TRUE
}
\seealso{
{
   \link{do_comparison}
}
}
