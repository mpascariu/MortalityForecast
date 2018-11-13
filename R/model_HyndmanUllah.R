
#' The Functional Demographic Model
#' @inheritParams doMortalityModels
#' @inheritParams demography::fdm
#' @inherit demography::fdm details return
#' @seealso 
#' \code{\link{fit_LeeCarter}}
#' \code{\link{fit_Oeppen}}
#' @details \insertNoCite{hyndman2007}{MortalityForecast}
#' @references \insertAllCited{}
#' @export
fit_HyndmanUllah <- function(data, x, y, order = 1, transform = TRUE, ...) {
  demo_data <- demogdata(data = data, ages = x, years = y, 
                         pop = data * 0, label = "demography", 
                         name = "mean", lambda = 0, type = "mortality")
  FDMfit <- fdm(demo_data, order = order, transform = transform, ...)
  
  dimnames(FDMfit$fitted$y) <- list(x, y)
  
  return(FDMfit)
}

