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

