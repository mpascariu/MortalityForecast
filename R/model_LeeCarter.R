# Marius D. Pascariu --- Mon Nov 12 14:00:33 2018 ------------------------------

#' The Lee-Carter Mortality Model
#' @param data A data.frame or a matrix containing mortality data 
#' with ages \code{x} as row and time \code{y} as column.
#' @inheritParams doMortalityModels
#' @inheritParams demography::lca
#' @inherit demography::lca details return
#' @seealso 
#' \code{\link{fit_HyndmanUllah}}
#' \code{\link{fit_Oeppen}}
#' @details \insertNoCite{lee1992}{MortalityForecast}
#' @references \insertAllCited{}
#' @export
fit_LeeCarter <- function(data, x, y, adjust = "none", ...) {
  demo_data <- demogdata(data = data, ages = x, years = y, 
                         pop = data * 0, label = "demography", 
                         name = "mean", lambda = 0, type = "mortality")
  LCfit <- lca(demo_data, restype = "rates", adjust = adjust, ...)
  
  dimnames(LCfit$fitted$y) <- list(x, y)
  
  return(LCfit)
}


#' Residuals of the Lee-Carter Mortality Model
#' @param object An object of class \code{"lca"}
#' @param ... Further arguments passed to or from other methods.
#' @export
residuals.lca <- function(object, ...){
  r <- object$residuals$y
  colnames(r) <- object$year
  structure(class = "residMF", as.matrix(r))
}
