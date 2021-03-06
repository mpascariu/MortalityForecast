# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Fri Nov 23 17:43:11 2018
# --------------------------------------------------- #


#' The Functional Demographic Model
#' @inheritParams do.MortalityModels
#' @inheritParams demography::fdm
#' @inherit demography::fdm details
#' @inherit model.Oeppen return
#' @seealso 
#' \code{\link{predict.HyndmanUllah}}
#' @details \insertNoCite{hyndman2007}{MortalityForecast}
#' @references \insertAllCited{}
#' @examples 
#' # Data
#' x  <- 0:89
#' y  <- 1985:2014
#' mx <- HMD_male$mx$GBRTENW[paste(x), paste(y)]
#' 
#' M <- model.HyndmanUllah(data = mx, x = x, y = y) # fit
#' P <- predict(M, h = 20)  # forecast
#' P
#' @export
model.HyndmanUllah <- function(data, 
                               x, 
                               y, 
                               order = 1, 
                               transform = TRUE, 
                               ...) {
  
  input <- c(as.list(environment()))
  x <- x %||% 1:nrow(data)
  y <- y %||% 1:ncol(data)
  
  # Info
  modelLN <- "Functional Demographic Model -- Hyndman-Ullah"   # long name
  modelSN <- "HyndmanUllah"                   # short name
  modelF  <- "log m[x,t] = a[x] + SUM b[x,k]phi[t,k] + r[t,x]eps[x,t]" # formula
  info <- list(name = modelLN, name.short = modelSN, formula = modelF)
  
  demo_data <- demogdata(data = data, ages = x, years = y, 
                         pop = data * 0, label = "demography", 
                         name = "mean", lambda = 0, type = "mortality")
  M <- fdm(demo_data, order = order, transform = transform, ...)
  
  fv <- exp(fitted(M)$y)
  dimnames(fv) <- list(x, y)
  resid <- data - fv
  cf <- NULL
  
  # Exit
  out <- list(input = input, 
              info = info, 
              call = match.call(), 
              fitted.values = fv, 
              observed.values = data,
              coefficients = cf, 
              residuals = resid, 
              x = x, 
              y = y, 
              demography = M)
  out <- structure(class = "HyndmanUllah", out)
  return(out)
}


#' Forecast age-specific death rates using the Hyndman-Ullah mortality model
#' 
#' @param object An object of class \code{HyndmanUllah}.
#' @inheritParams predict.Oeppen
#' @inherit predict.Oeppen return
#' @seealso 
#' \code{\link{model.HyndmanUllah}}
#' @author Marius D. Pascariu and Marie-Pier Bergeron-Boucher
#' @examples # For examples go to ?model.HyndmanUllah
#' @details \insertNoCite{hyndman2007}{MortalityForecast}
#' @references \insertAllCited{}
#' @export
predict.HyndmanUllah <- function(object, 
                                 h, 
                                 level = 95, 
                                 jumpchoice = c("actual", "fit"), 
                                 verbose = TRUE, 
                                 ...){
  
  x <- object$x
  
  # Timeline
  bop <- max(object$y) + 1
  eop <- bop + h - 1
  fcy <- bop:eop
  
  J <- match.arg(jumpchoice)
  P <- forecast(object$demography, h = h, jumpchoice = J, level = level)$rate
  m <- P[c("mean", "lower", "upper")] # forecast mx
  names(m) <- c('mean', paste0('L', level), paste0('U', level))
  
  # Exit
  out <- list(call = match.call(), 
              predicted.values = m[[1]],
              conf.intervals = m[-1], 
              x = x, 
              y = fcy, 
              info = object$info)
  out <- structure(class = 'predict.HyndmanUllah', out)
  return(out)
}


#' @rdname residuals.Oeppen
#' @export
residuals.HyndmanUllah <- function(object, ...){
  residuals_default(object, ...)
}


#' @rdname print_default
#' @export
print.HyndmanUllah <- function(x, ...) {
  print_default(x, ...)
}


#' @rdname print_default
#' @export
print.predict.HyndmanUllah <- function(x, ...) {
  print_predict_default(x, ...)
}



