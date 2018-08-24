

#' The Lee-Carter Mortality Model
#' @inheritParams doMortalityModels
#' @inherit demography::lca details return
#' @seealso 
#' \code{\link{fitHyndmanUllah}}
#' \code{\link{fitOeppen}}
#' @details \insertNoCite{lee1992}{MortalityForecast}
#' @references \insertAllCited{}
#' @export
fitLeeCarter <- function(data, x, y, ...) {
  demo_data <- demogdata(data = data, ages = x, years = y, 
                         pop = data * 0, label = "demography", 
                         name = "mean", lambda = 0, type = "mortality")
  LCfit <- lca(demo_data, ...)
  
  dimnames(LCfit$fitted$y) <- list(x, y)
  
  return(LCfit)
}


#' The Functional Demographic Model
#' @inheritParams doMortalityModels
#' @inherit demography::fdm details return
#' @seealso 
#' \code{\link{fitLeeCarter}}
#' \code{\link{fitOeppen}}
#' @details \insertNoCite{hyndman2007}{MortalityForecast}
#' @references \insertAllCited{}
#' @export
fitHyndmanUllah <- function(data, x, y, ...) {
  demo_data <- demogdata(data = data, ages = x, years = y, 
                         pop = data * 0, label = "demography", 
                         name = "mean", lambda = 0, type = "mortality")
  FDMfit <- fdm(demo_data, ...)
  
  dimnames(FDMfit$fitted$y) <- list(x, y)
  
  return(FDMfit)
}

