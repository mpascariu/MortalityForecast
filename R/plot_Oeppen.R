# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Fri Nov 30 22:15:23 2018
# --------------------------------------------------- #


#' @rdname plot.MEM
#' @export
plot.Oeppen <- function(x, 
                        plotType = c("fitted", "observed"), 
                        ny = 7, 
                        level = 80, 
                        ...){
  plot.MEM(x, plotType, ny, level, ...)
}


#' Plot the Predicted Age-at-Death Distribution
#' 
#' @param x An object generated using the \code{predict} function;
#' @param plotType The type of the plot. The alternatives are 
#' \code{"mean", "lower", "upper"}. Default: \code{"mean"}.
#' @inheritParams plot.predict.MEM
#' @seealso 
#' \code{\link{model.Oeppen}}
#' \code{\link{model.OeppenC}}
#' @examples # For examples go to ?model.Oeppen or ?model.OeppenC
#' @export
plot.predict.Oeppen <- function(x, 
                                plotType = c("mean", "lower", "upper"), 
                                ny = 7, 
                                level = 80, 
                                ...) 
{
  plotType <- match.arg(plotType)
  if (plotType == "mean") {
    mat <- x$predicted.values
    P <- ggplotDistribConvergence(mat, x = x$x, y = x$y, ny, level) + 
      labs(subtitle = "Forecast Values - Best estimate")
    
  } else if (plotType == "lower") {
    mat <- x$conf.intervals[[1]]
    P <- ggplotDistribConvergence(mat, x = x$x, y = x$y, ny, level) + 
      labs(subtitle = "Forecast Values - lower bound")
    
  } else if (plotType == "upper") {
    mat <- x$conf.intervals[[length(x$conf.intervals)/2 + 1]]
    P <- ggplotDistribConvergence(mat, x = x$x, y = x$y, ny, level) + 
      labs(subtitle = "Forecast Values - upper bound")
    
  } 
  suppressMessages(print(P))
}

# ------------------------------------------

#' @rdname plot.MEM
#' @export
plot.OeppenC <- function(x, 
                         plotType = c("fitted", "observed"), 
                         ny = 7, 
                         level = 80, 
                         ...){
  plot.Oeppen(x, plotType, ny, level, ...)
}


#' @rdname plot.predict.Oeppen
#' @export
plot.predict.OeppenC <- function(x, 
                                 plotType = c("mean", "lower", "upper"), 
                                 ny = 7, 
                                 level = 80, 
                                 ...) {
  plot.predict.Oeppen(x, plotType, ny, level, ...)
}





