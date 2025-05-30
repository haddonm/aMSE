% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defineZone.R
\name{maturity}
\alias{maturity}
\title{maturity Logistic maturity curve}
\usage{
maturity(ina, inb, lens)
}
\arguments{
\item{ina}{is the intercept of the exponential function}

\item{inb}{is the gradient of the exponential function}

\item{lens}{a vector of lengths for which the logistic maturity value
will be calculated}
}
\value{
A vector of length(lens) containing the predicted maturity at
    length values
}
\description{
maturity this uses the logistic function:
  exp(a+bL)/(1+exp(a+bL)), which has the property that the SM50 = -a/b
  and the interquartile distance is 2.Ln(3)/b.
}
\examples{
\dontrun{
a <- -14.383
b <- 0.146017
lens <- seq(2,210,2)
Maturity <- maturity(a,b,lens)
}
}
