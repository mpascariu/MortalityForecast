

#' The Lee-Carter Mortality Model
#' @inheritParams doMortalityModels
#' @inherit demography::lca details return
#' @seealso 
#' \code{\link{fit_HyndmanUllah}}
#' \code{\link{fit_Oeppen}}
#' @details \insertNoCite{lee1992}{MortalityForecast}
#' @references \insertAllCited{}
#' @export
fit_LeeCarter <- function(data, x, y, ...) {
  demo_data <- demogdata(data = data, ages = x, years = y, 
                         pop = data * 0, label = "demography", 
                         name = "mean", lambda = 0, type = "mortality")
  LCfit <- lca(demo_data, ...)
  
  dimnames(LCfit$fitted$y) <- list(x, y)
  
  return(LCfit)
}

