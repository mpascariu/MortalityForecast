% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MEM_Plot_fit.R
\name{fitted2dens}
\alias{fitted2dens}
\title{Prepare data for ggplots in plot.fitMaxEntMortality function}
\usage{
fitted2dens(mat, ny = 7, lx0 = 300)
}
\arguments{
\item{mat}{Matrix containing the observed or fitted value.}

\item{ny}{Number of years to be selected from input data and to be added in the plot.}

\item{lx0}{Adjusting parameter.}
}
\description{
Prepare data for ggplots in plot.fitMaxEntMortality function
}
\examples{
x  <- 0:110
y  <- 1965:2014
dx <- MortalityForecast.data$dx[paste(x), paste(y)]
M  <- fitMaxEntMortality(dx, x, y, n = 5)
fitted2dens(fitted(M))
}
\keyword{internal}