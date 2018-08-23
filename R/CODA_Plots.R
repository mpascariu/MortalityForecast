

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

