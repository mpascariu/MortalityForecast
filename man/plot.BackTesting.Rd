% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_doBackTesting.R
\name{plot.BackTesting}
\alias{plot.BackTesting}
\title{Plot method for objects of the class \code{BackTesting}}
\usage{
\method{plot}{BackTesting}(x, data.out, facet = c("x", "y"),
  which = NULL, ...)
}
\arguments{
\item{x}{An object of the class \code{do.BackTesting}.}

\item{data.out}{Specify the type of data to be returned in output. 
Various life table indices are accepted: 
\code{"qx", "mx", "dx", "lx", "Lx", "Tx", "ex"}.}

\item{facet}{What facets to include? Options: \code{"x", "y"}.}

\item{which}{Which \code{x} or which \code{y} to be plotted. Numerical vector.}

\item{...}{Arguments to be passed to or from other methods.}
}
\description{
Plot method for objects of the class \code{BackTesting}
}
