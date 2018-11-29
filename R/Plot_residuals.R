# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Thu Nov 29 13:12:51 2018
# --------------------------------------------------- #


#' Plot the deviance residuals
#' 
#' Plots the deviance residuals of a Mortality Model model. 
#' Three types of plots are available: scatter plot of residuals by age, 
#' period and cohort, colour map (heatmap) of the residuals, and a black and 
#' white signplot of the residuals.
#' 
#' @param x An object of class \code{residMF}.
#' @param plotType The type of the plot. The alternatives are 
#' \code{"scatter"}(default), \code{"colourmap"}, and \code{"signplot"}.
#' @param reslim Optional numeric vector of length 2, giving the range of the 
#' residuals.
#' @param plotAge Logical value indicating if the age scatter plot should be 
#' produced. This is only used when \code{plotType = "scatter"}.
#' @param plotYear Logical value indicating if the calendar year scatter plot 
#' should be produced. This is only used when \code{plotType = "scatter"}.
#' @param plotCohort Logical value indicating if the cohort scatter plot 
#' should be produced. This is only used when \code{plotType = "scatter"}.
#' @param pch optional symbol to use for the points in a scatterplot. 
#' This is only used when \code{plotType = "scatter"}. See 
#' \code{\link[graphics]{plot}}.
#' @param col Optional colours to use in plotting. If 
#' \code{plotType = "scatter"} this is a single colour to use in the points
#' in the scatter plots, while if \code{plotType = "colourmap"} this should
#' be a list of colours (see help in \code{\link[fields]{image.plot}} 
#' for details). This argument is ignored if \code{plotType = "signplot"}.
#' @param ... Other plotting parameters to be passed to the plotting 
#' functions. This can be used to control the appearance of the plots.
#' @details
#' When \code{plotType = "scatter"} scatter plots of the residuals against age, 
#' calendar year and cohort (year of birth) are produced. 
#'
#' When \code{plotType = "colourmap"} a two dimensional colour map of the 
#' residuals is plotted. This is produced using function 
#' \code{\link[fields]{image.plot}}. See \code{\link[fields]{image.plot}} 
#' for further parameters that can be passed to this type of plots.
#'
#' When \code{plotType = "signplot"} a two dimensional black and white map of the
#' residuals is plotted with dark grey representing negative residuals and 
#' light grey representing positive residuals. This is produced using 
#' function \code{\link[graphics]{image.default}}. 
#'   
#' @source 
#' The code for producing the residual plots is inspired from the one published 
#' in \href{http://github.com/amvillegas/StMoMo}{\code{StMoMo}} R package. 
#' See \code{\link[StMoMo]{plot.resStMoMo}}. 
#' All the credit goes to it's authors: 
#' \href{https://github.com/amvillegas}{Andres Villegas}, Pietro Millossovich 
#' and Vladimir Kaishev.
#' @keywords internal
#' @export 
plot.residMF <- function(x, 
                         plotType = c("scatter", "colourmap", "signplot"), 
                         reslim = NULL, 
                         plotAge = TRUE, 
                         plotYear = TRUE, 
                         plotCohort = TRUE, 
                         pch = 20, 
                         col = NULL, 
                         ...) {
  
  plotType <- match.arg(plotType)
  oldpar   <- par(no.readonly = TRUE)
  L        <- dim(x)
  res      <- as.data.frame(x[1:L[1], 1:L[2]])
  ages     <- as.numeric(rownames(res))
  years    <- as.numeric(colnames(res))
  
  if (is.null(reslim)) {
    maxRes <- max(abs(res), na.rm = TRUE)
    reslim <- c(-maxRes, maxRes)
  }
  if (is.null(col) & plotType == "colourmap") {
    col <- colorRampPalette(brewer.pal(10, "RdBu"))(64)
  }
  if (is.null(col) & plotType == "scatter") col = "black"
  
  switch(plotType, 
         scatter = scatterplotAPC(res, ages, years, 
                                  plotAge = plotAge, plotYear = plotYear, 
                                  plotCohort = plotCohort, pch = pch, 
                                  ylab = "residuals", ylim = reslim, 
                                  col = col, ...),
         colourmap = image.plot(x = years, y = ages, z = t(res), 
                                zlim = reslim, ylab = "age", 
                                xlab = "calendar year", col = col, ...),
         signplot = image.default(x = years, y = ages, z = t(res), 
                                  zlim = reslim, ylab = "age", 
                                  xlab = "calendar year", 
                                  breaks = c(-10e10, 0, 10e10), 
                                  col = grey.colors(2), ...)
  )
  message(paste('Different types of plots can be obtained by using,\nplotType: ', 
                paste(c("scatter", "colourmap", "signplot"), collapse = ", ")))
  par(oldpar)
}


#'Do a scatter plot of a matrix according to age-period-cohorts
#'
#' @param mat Matrix with the data to plot.
#' @param ages Ages corresponding to the rows in \code{mat}.
#' @param years Years corresponding to the columns in \code{mat}.  
#' @param plotAge Logical value indicating if the age scatter plot should be 
#' produced.
#' @param plotYear Logical value indicating if the calendar year scatter plot 
#' should be produced.
#' @param plotCohort Logical value indicating if the cohort scatter plot 
#' should be produced.
#' @param zeroLine Logical value indicating if a horizontal line at zero
#' should be plotted.
#' @inheritParams plot.residMF
#' @keywords internal
scatterplotAPC <- function(mat, 
                           ages, 
                           years, 
                           plotAge = TRUE, 
                           plotYear = TRUE, 
                           plotCohort = TRUE, 
                           zeroLine  = TRUE, 
                           ...) {
  
  nAges  <- length(ages)
  nYears <- length(years)  
  if (nrow(mat) != nAges ||  ncol(mat) != nYears) {
    stop(paste("Mismatch between the dimensions of the plot in data and the",
               "number of years or ages"), call. = FALSE)
  }
  dimnames(mat) <- list(ages, years)
  x = y = t <- NULL # hack to remove note in CRAN check
  mat    <- cbind(x = ages, mat)
  data   <- tidyr::gather(mat, key = t, value = y, -x)
  data$t <- as.numeric(data$t)
  data   <- transform(data, c = t - x) 
  N      <- plotAge + plotYear + plotCohort
  
  if (N > 0) par(mfrow = c(1, N))
  
  if (plotAge) {
    plot(data$x, data$y, type = "p", xlab = "age", ...)
    if (zeroLine) abline(h = 0) 
  }
  if (plotYear) {
    plot(data$t, data$y, type = "p", xlab = "calendar year", ...)
    if (zeroLine) abline(h = 0) 
  }
  if (plotCohort) {
    plot(data$c, data$y, type = "p", xlab = "year of birth", ...)
    if (zeroLine) abline(h = 0) 
  }
}

