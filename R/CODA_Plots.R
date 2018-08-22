

#' Plot fitted parameters and fitted values of a CoDa mortality model
#' @description Plot fitted parameters and fitted values of a CoDa mortality model. 
#' Two types of plots are available: \code{"coef"} to obtain representations of the 
#' three estimated series of parameters and \code{"data"} for visualising the input 
#' and fitted values.
#' @details 
#' When \code{"plotType = coef"} the figure is produced using the basic 
#' plot R function \code{\link[graphics]{plot}}.
#' 
#' When \code{"plotType = data"} the figure is produced using function 
#' \code{\link[graphics]{matplot}}.
#' 
#' This is important to know in order to use the \code{"..."} argument adequately.
#' @param x An object of class \code{"coda"} with the fitted 
#' parameters of the mortality model.
#' @param plotType The type of the plot. The alternatives are 
#' \code{"coef"}(default) and \code{"data"}.
#' @inheritParams graphics::plot.default
#' @examples 
#' # Fit model
#' D <- MortalityForecast.data$dx
#' M <- coda(D, x = 0:110, y = 1960:2016)
#' 
#' # Plot fitted parameters
#' plot(M, plotType = "coef")
#' plot(M, plotType = "coef", type = "p", pch = 19)
#' 
#' # Plot input data and fitted values
#' plot(M, plotType = "data")
#' @export
#' 
plot.coda <- function(x, plotType = c("coef", "data"), 
                      type = "l", ylim = NULL, ylab = "", ...){
  oldpar <- par(no.readonly = TRUE)
  age    <- x$x
  year   <- x$y
  C      <- coef(x)
  plotType  <- match.arg(plotType)
  
  plot_coef <- function(){
    par(mfrow = c(1, 3))
    plot(age, C$ax, type = type, main = "ax", xlab = "age", ylab = ylab, ...)
    plot(age, C$bx, type = type, main = "bx", xlab = "age", ylab = ylab, ...)
    plot(year, C$kt, type = type, main = "kt", xlab = "year", ylab = ylab, ...)
  }
  plot_data <- function(){
    fdx  <- fitted(x)
    idx  <- x$input$data
    if (is.null(ylim)) ylim <- range(idx, fdx)
    par(mfrow = c(1, 2))
    matplot(idx, type = type, ylim = ylim, main = "Input data", 
            xlab = "age", ylab = ylab,
            col = terrain.colors(ncol(idx), alpha = 1), ...)
    matplot(fdx, type = type, ylim = ylim, main = "Fitted values", 
            xlab = "age", ylab = ylab,
            col = terrain.colors(ncol(idx), alpha = 1), ...)
  }
  switch(plotType,
         coef = plot_coef(),
         data = plot_data())
  par(oldpar)
}


#' Plot the residuals of a CoDa Mortality Model
#' 
#' Plot the deviance residuals of a CoDa Mortality Model which are 
#' of class \code{"residuals.coda"}. Three types of plots
#' are available: scatter plot of residuals by age, period and cohort,
#' colour map (heatmap) of the residuals, and a black and white signplot 
#' of the residuals.
#' 
#' @param x An object of class \code{residuals.coda} with the residuals of a 
#' CoDa Mortality Model.
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
#' @param pch Optional symbol to use for the points in a scatterplot. 
#' This is only used when \code{plotType = "scatter"}. See 
#' \code{\link[graphics]{plot}}.
#' @param col Optional colours to use in plotting. If 
#' \code{plotType = "scatter"} this is a single colour to use in the points
#' in the scatter plots, while if \code{plotType = "colourmap"} this should
#' be a list of colours (see help in \code{\link[fields]{image.plot}} 
#' for details). This argument is ignored if \code{plotType = "signplot"}.
#' @param ... Other plotting parameters to be passed to the plotting 
#' functions. This can be used to control the appearance of the plots.
#'
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
#'  residuals is plotted with dark grey representing negative residuals and 
#'  light grey representing positive residuals. This is produced using 
#'  function \code{\link[graphics]{image.default}}. 
#'   
#' @examples
#' # Fit model
#' D <- MortalityForecast.data$dx
#' M <- coda(D, x = 0:110, y = 1960:2016)
#' 
#' # Plot residuals
#' res <- resid(M)
#' plot(res, plotType = "scatter")
#' plot(res, plotType = "colourmap")
#' plot(res, plotType = "signplot")
#' 
#' @source 
#' The code for producing the residual plots is inspired from the one published 
#' in \href{http://github.com/amvillegas/StMoMo}{\code{StMoMo}} R package. 
#' See \code{\link[StMoMo]{plot.resStMoMo}}. 
#' All the credit goes to it's authors: 
#' \href{https://github.com/amvillegas}{Andres Villegas}, Pietro Millossovich 
#' and Vladimir Kaishev.
#' @export 
#' 
plot.residuals.coda <- function(x, plotType = c("scatter", "colourmap", "signplot"), 
                                reslim = NULL, plotAge = TRUE, plotYear = TRUE, 
                                plotCohort  = TRUE, pch = 20, col = NULL, ...) {
  ages  <- as.numeric(rownames(x))
  years <- as.numeric(colnames(x))
  L     <- dim(x)
  res   <- as.data.frame(x[1:L[1], 1:L[2]])
  plotType <- match.arg(plotType)
  
  if (is.null(reslim)) {
    maxRes <- max(abs(res), na.rm = TRUE)
    reslim <- c(-maxRes, maxRes)
  }
  if (is.null(col) & plotType == "colourmap") {
    col <- colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(64)
  }
  if (is.null(col) & plotType == "scatter") {
    col <- "black"
  }
  
  oldpar <- par(no.readonly = TRUE)
  switch(plotType, 
         scatter = scatterplotAPC(res, ages, years, 
                                  plotAge = plotAge, plotYear = plotYear, 
                                  plotCohort = plotCohort, pch = pch, 
                                  ylab = "residuals", ylim = reslim, col = col, 
                                  ...),
         colourmap = fields::image.plot(years, ages, t(res), zlim = reslim, 
                                        ylab = "age", xlab = "calendar year", 
                                        col = col, ...),
         signplot = image.default(years, ages, t(res), zlim = reslim, 
                                  ylab = "age", xlab = "calendar year", 
                                  breaks = c(-10e10, 0, 10e10), 
                                  col = grey.colors(2), ...)
  )
  message(paste('Different types of plots can be obtained by using,\nplotType: ', 
                paste(c("scatter", "colourmap", "signplot"), collapse = ", ")))
  par(oldpar)
}


#' Do a scatter plot of a matrix according to age-period-cohorts
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
#' @param ... Other arguments to pass to the plot function.
#' @keywords internal
#' 
scatterplotAPC <- function(mat, ages, years, plotAge = TRUE, plotYear = TRUE, 
                           plotCohort  = TRUE, zeroLine  = TRUE, ...) {
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

