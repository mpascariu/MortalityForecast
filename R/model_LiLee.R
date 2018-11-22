# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Mon Nov 19 13:58:32 2018
# --------------------------------------------------- #

#' The Li-Lee Mortality Model
#' 
#' Fit the Li-Lee mortality model
#' @inheritParams doMortalityModels
#' @param benchmark A list containing multiple matrices with mortality data 
#' to be used as benchmark. Format: ages \code{x} as row and time \code{y} as column.
#' @seealso 
#' \code{\link{predict.LiLee}}
#' \code{\link{model_LeeCarter}}
#' @details \insertNoCite{li2005}{MortalityForecast}
#' @references \insertAllCited{}
#' @export
model_LiLee <- function(data, benchmark, x = NULL, y = NULL, verbose = TRUE, ...){
  input <- c(as.list(environment()))
  if (any(data == 0)) {
    stop("The input data contains death rates equal to zero at various ages.")
  }
  x <- x %||% 1:nrow(data)
  y <- y %||% 1:ncol(data)

  # Info
  modelLN <- "Li-Lee Mortality Model"       # long name
  modelSN <- "LL"                           # short name
  modelF  <- "log m[x,t] = a[x] + b'[x]k'[t] + e[x,t]" # formula
  info <- list(name = modelLN, name.short = modelSN, formula = modelF)
  
  # Fit benchmark model
  B <- model_LeeCarter(data = benchmark, x = x, y = y)
  B.cmx <- with(B$coefficients, c(bx) %*% t(kt))
  
  # Estimate model parameters: a[x], b[x], k[t]
  ax  <- apply(log(data), 1, mean)
  A.cmx <- sweep(log(data), 1, ax, FUN = "-")
  cmx <- A.cmx - B.cmx
  S   <- svd(cmx)
  kt  <- S$v[, 1] * sum(S$u[, 1]) * S$d[1]
  bx  <- S$u[, 1] / sum(S$u[, 1])
  cf  <- list(ax = as.numeric(ax), bx = as.numeric(bx), kt = as.numeric(kt))
  
  # Variability
  var <- cumsum(S$d^2 / sum(S$d^2))[1]
  
  # Compute fitted values and devinace residuals based on the estimated model
  fv <- sweep(c(bx) %*% t(kt) + B.cmx, 1, ax, FUN = "+") # Fitted values
  fv <- exp(fv) # fitted mx
  dimnames(fv) <- list(x, y)
  
  resid <- data - fv # residuals
  
  # Exit
  out <- list(input = input, info = info, call = match.call(), 
              fitted.values = fv, observed.values = data,
              coefficients = cf, residuals = resid, x = x, y = y,
              benchmark.LeeCarter = B)
  out <- structure(class = 'LiLee', out)
  return(out)
}


#' Forecast age-specific death rates using the Li-Lee model
#' @param object An object of class \code{LiLee}.
#' @inheritParams predict.Oeppen
#' @inherit predict.Oeppen return
#' @seealso 
#' \code{\link{model_LiLee}}
#' @author Marius D. Pascariu and Marie-Pier Bergeron-Boucher
#' @examples 
#' # Data
#' x  <- 0:89
#' y  <- 1985:2014
#' B.mx <- HMD_male$mx$USA[paste(x), paste(y)]
#' mx <- HMD_male$mx$GBRTENW[paste(x), paste(y)]
#' 
#' M <- model_LiLee(data = mx, x = x, y = y, benchmark = B.mx) # fit
#' P <- predict(M, h = 20)  # forecast
#' P
#' @export
predict.LiLee <- function(object,
                          h,
                          order = c(0, 1, 0),
                          include.drift = TRUE,
                          level = c(80, 95),
                          jumpchoice = c("actual", "fit"),
                          method = "ML",
                          verbose = TRUE, ...){
  
  # Benchmark Lee-Carter forecast
  B <- object$benchmark.LeeCarter
  B.pred <- predict(object = B, h, order, include.drift, level, 
                    jumpchoice, method, verbose = FALSE)
  B.cmx <- coef(B)$bx %*% t(B.pred$kt$mean)
  
  # Data
  mx <- t(object$input$data)
  
  # Timeline
  bop <- max(object$y) + 1
  eop <- bop + h - 1
  fcy <- bop:eop
  
  # Identify the k[t] ARIMA order
  C <- coef(object)
  A <- find_arima(C$kt)
  
  # Estimate/fit k[t] time-series model
  kt.arima <- forecast::Arima(y = C$kt, 
                              order = order %||% A$order, 
                              include.drift = include.drift %||% A$drift,
                              method = method)
  
  # Forecast k[t] using the time-series model
  tsf <- forecast(kt.arima, h = h + 1, level = level)  # time series forecast
  fkt <- data.frame(tsf$mean, tsf$lower, tsf$upper) # forecast kt
  Cnames <- c('mean', paste0('L', level), paste0('U', level))
  dimnames(fkt) <- list(c(0, fcy), Cnames)
  
  # Get forecast m[x] based on k[t] extrapolation 
  # Here we are also adjusting for the jump-off
  fmx <- get_mx_values(object = object,
                       kt = fkt, 
                       jumpchoice = match.arg(jumpchoice), 
                       y = fcy, 
                       adj = t(B.cmx))
  pv <- fmx[[1]]
  CI <- fmx[-1]
  
  # Exit
  out <- list(call = match.call(), predicted.values = pv, 
              kt.arima = kt.arima, kt = fkt, 
              conf.intervals = CI, x = object$x, y = fcy, info = object$info,
              benchmark_LeeCarter = B.pred)
  out <- structure(class = "predict.LiLee", out)
  return(out)
}


# S3 ----------------------------------------------

#' Residuals of the Li-Lee Mortality Model
#' @param object An object of class \code{"LeeCarter2"}
#' @inheritParams residuals_default
#' @export
residuals.LiLee <- function(object, ...){
  residuals_default(object, ...)
}


#' Print Li-Lee model
#' @param x An object of class \code{"LeeCarter2"}
#' @inheritParams print_default
#' @keywords internal
#' @export
print.LiLee <- function(x, ...) {
  print_default(x, ...)
}


#' Summary LiLee
#' @param object An object of class \code{"LeeCarter2"}.
#' @inheritParams print.LiLee
#' @keywords internal
#' @export
summary.LiLee <- function(object, ...) {
  axbx <- data.frame(ax = object$coefficients$ax, 
                     bx = object$coefficients$bx,
                     row.names = object$x)
  kt <- data.frame(kt = object$coefficients$kt)
  out = structure(class = "summary.LiLee", 
                  list(A = axbx, K = kt, call = object$call, info = object$info,
                       y = object$y, x_ = object$x))
  return(out)
}


#' Print summary.LiLee
#' @param x An object of class \code{"summary.LiLee"}.
#' @inheritParams print.LiLee
#' @keywords internal
#' @export
print.summary.LiLee <- function(x, ...){
  cat('\nFit  :', x$info$name)
  cat('\nModel:', x$info$formula)
  cat('\n\nCoefficients:\n')
  A <- head_tail(x$A, digits = 5, hlength = 6, tlength = 6)
  K <- head_tail(data.frame(. = '|', y = as.integer(x$y), kt = x$K),
                 digits = 5, hlength = 6, tlength = 6)
  print(data.frame(A, K))
  cat('\n')
}


#' Print function
#' @param x An object of class \code{"predict.LiLee"};
#' @inheritParams print_predict_default
#' @keywords internal
#' @export
print.predict.LiLee <- function(x, ...) {
  print_predict_default(x, ...)
  cat('k[t]-ARIMA method:', arima.string1(x$kt.arima, padding = TRUE))
  cat('\n')
}
