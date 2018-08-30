# Fri Aug 24 12:35:14 2018 --------- Marius D. Pascariu ---


#' ggplot the observed and fitted values of a CoDa-LC mortality model
#' 
#' @inherit plot.fitMaxEntMortality details
#' @inheritParams plot.fitMaxEntMortality
#' @examples 
#' # For examples go to ?fitOeppen
#' @export
plot.fitOeppen <- function(x, plotType = c("fitted", "observed"), 
                           ny = 7, level = 80, ...){
  plot.fitMaxEntMortality(x, plotType, ny, level, ...)
}


#' ggplot the predicted values of a CoDa-LC mortality model
#' 
#' @param x An object of the class \code{\link{predict.fitOeppen}}.
#' @param plotType The type of the plot. The alternatives are 
#' \code{"mean", "lower", "upper"}. Default: \code{"mean"}.
#' @inheritParams plot.predict.fitMaxEntMortality
#' @examples 
#' # For examples go to ?predict.fitOeppen
#' @export
plot.predict.fitOeppen <- function(x, plotType = c("mean", "lower", "upper"), 
                                   ny = 7, level = 80, ...) 
{
  plotType <- match.arg(plotType)
  if (plotType == "mean") {
    mat = x$predicted.values
    P <- ggplotDistribConvergence(mat, x = x$x, y = x$y, ny, level) + 
      labs(subtitle = "Forecast Values - Best estimate")
    
  } else if (plotType == "lower") {
    mat = x$conf.intervals[[1]]
    P <- ggplotDistribConvergence(mat, x = x$x, y = x$y, ny, level) + 
      labs(subtitle = "Forecast Values - lower bound")
    
  } else if (plotType == "upper") {
    mat = x$conf.intervals[[length(x$conf.intervals)/2 + 1]]
    P <- ggplotDistribConvergence(mat, x = x$x, y = x$y, ny, level) + 
      labs(subtitle = "Forecast Values - upper bound")
    
  } 
  suppressMessages(print(P))
}
