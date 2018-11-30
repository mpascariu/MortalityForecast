# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Mon Nov 19 13:43:17 2018
# --------------------------------------------------- #

#' The Lee-Carter Mortality Model
#' 
#' Fit the Lee-Carter mortality model
#' @inheritParams do.MortalityModels
#' @inherit model.Oeppen return
#' @seealso 
#' \code{\link{predict.LeeCarter}}
#' @details \insertNoCite{lee1992}{MortalityForecast}
#' @references \insertAllCited{}
#' @examples 
#' # Data
#' x  <- 0:89
#' y  <- 1985:2014
#' mx <- HMD_male$mx$GBRTENW[paste(x), paste(y)]
#' 
#' # Fit the model
#' M <- model.LeeCarter(data = mx, x = x, y = y)
#' M
#' summary(M)
#' 
#' # Check residuals
#' R <- residuals(M)
#' 
#' plot(R, plotType = "scatter")
#' plot(R, plotType = "colourmap")
#' plot(R, plotType = "signplot")
#' 
#' # Forecast 
#' P <- predict(M, h = 20)
#' P
#' @export
model.LeeCarter <- function(data, 
                            x = NULL, 
                            y = NULL, 
                            verbose = TRUE, 
                            ...){
  
  input <- c(as.list(environment()))
  if (any(data == 0)) {
    stop("The input data contains death rates equal to zero at various ages.")
  }
  x <- x %||% 1:nrow(data)
  y <- y %||% 1:ncol(data)
  
  # Info
  modelLN <- "Lee-Carter Mortality Model"   # long name
  modelSN <- "LC"                           # short name
  modelF  <- "log m[x,t] = a[x] + b[x]k[t]" # formula
  info <- list(name = modelLN, name.short = modelSN, formula = modelF)
  
  # Estimate model parameters: a[x], b[x], k[t]
  ax   <- apply(log(data), 1, mean)
  cmx  <- sweep(log(data), 1, ax, FUN = "-")
  S    <- svd(cmx)
  kt   <- S$v[,1] * sum(S$u[, 1]) * S$d[1]
  bx   <- S$u[,1] / sum(S$u[, 1])
  cf   <- list(ax = as.numeric(ax), bx = as.numeric(bx), kt = as.numeric(kt))
  
  # Variability
  var <- cumsum(S$d^2 / sum(S$d^2))[1]
  
  # Compute fitted values and devinace residuals based on the estimated model
  fv <- sweep(c(bx) %*% t(kt), 1, ax, FUN = "+") # Fitted values
  fv <- exp(fv) # fitted mx
  dimnames(fv) <- list(x, y)
  
  resid <- data - fv # residuals
  
  # Exit
  out <- list(input = input, 
              info = info, 
              call = match.call(), 
              coefficients = cf, 
              fitted.values = fv, 
              observed.values = data,
              residuals = resid, 
              x = x, 
              y = y)
  out <- structure(class = 'LeeCarter', out)
  return(out)
}


#' Forecast age-specific death rates using the Lee-Carter model.
#' 
#' @param object An object of class \code{LeeCarter}.
#' @inheritParams predict.Oeppen
#' @inherit predict.Oeppen return
#' @seealso 
#' \code{\link{model.LeeCarter}}
#' @author Marius D. Pascariu and Marie-Pier Bergeron-Boucher
#' @details \insertNoCite{lee1992}{MortalityForecast}
#' @references \insertAllCited{}
#' @examples # For examples go to ?model.LeeCarter
#' @export
predict.LeeCarter <- function(object, 
                               h, 
                               order = c(0, 1, 0), 
                               include.drift = TRUE,
                               level = c(80, 95),
                               jumpchoice = c("actual", "fit"),
                               method = "ML",
                               verbose = TRUE, ...){

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
  # Note: we have used h + 1 in order to extrapolate 1 more year and use it in 
  # the jump-off adjustment if needed. By rebasing the 1st forecast value to the 
  # last observed values. See the behaviour in get_mx_values(). 
  # The same is done in LL model.
  fkt <- data.frame(tsf$mean, tsf$lower, tsf$upper) # forecast kt
  Cnames <- c('mean', paste0('L', level), paste0('U', level))
  dimnames(fkt) <- list(c(0, fcy), Cnames)
  
  # Get forecast m[x] based on k[t] extrapolation 
  # Here we are also adjusting for the jump-off
  J <- match.arg(jumpchoice)
  m <- get_mx_values(object = object, 
                     kt = fkt, 
                     jumpchoice = J, 
                     y = fcy)
  
  # Exit
  out <- list(call = match.call(), 
              info = object$info,
              kt = fkt, 
              kt.arima = kt.arima, 
              predicted.values = m[[1]], 
              conf.intervals = m[-1], 
              x = object$x, 
              y = fcy)
  out <- structure(class = 'predict.LeeCarter', out)
  return(out)
}


#' Get m[x] values and confidence intervals based on k[t] forecast
#' In this function we compute the m[x] values based on the extrapolation of
#' the k[t] time-series. If necesary an adjustment for the jump-off is 
#' provided.
#' @inheritParams predict.LeeCarter 
#' @inheritParams model.LeeCarter
#' @param kt Predicted k[t] values in the model;
#' @param B.kt Predicted k[t] values of the benchmark model, used in the Li-Lee model only.
#' @keywords internal
get_mx_values <- function(object, jumpchoice, y, kt, B.kt = NULL){
  
  C  <- coef(object)
  OV <- object$observed.values
  N  <- ncol(OV)
  P  <- NULL
    
  for (i in 1:ncol(kt)) {
    
    # This is used only in LiLee model, and it is basically the trend 
    # given by the benchmark population
    if (is.null(B.kt)) { 
      B.cmx <- 0
    } else {
      B.bx <- coef(object$benchmark)$bx
      B.cmx <- c(B.kt[, i]) %*% t(B.bx)
    }
    
    # Compute predicted m[x] values
    p <- c(kt[, i]) %*% t(C$bx) + B.cmx
    p <- sweep(p, 2, C$ax, FUN = "+")
    p <- t(exp(p))
    
    # Here we adjust m[x] for jump-off if needed
    if (jumpchoice == 'actual') {
      J <- as.numeric(OV[, N]/p[, 1]) # jump_off (%)
      p <- sweep(p, 1, J, FUN = "*")
    }
    
    p <- p[, -1]
    dimnames(p) <- list(rownames(OV), y)
    P[[i]] <- p
    remove(p)
  }
  
  names(P) <- colnames(kt)
  return(P)
}


# S3 ----------------------------------------------


#' Residuals of the Lee-Carter Mortality Model
#' @param object An object of class \code{"LeeCarter"}
#' @inheritParams residuals_default
#' @examples # For examples go to ?model.LeeCarter
#' @export
residuals.LeeCarter <- function(object, ...){
  residuals_default(object, ...)
}


#' Print Lee-Carter model
#' @param x An object of class \code{"LeeCarter"}
#' @inheritParams print_default
#' @keywords internal
#' @export
print.LeeCarter <- function(x, ...) {
  print_default(x, ...)
}


#' Summary LeeCarter
#' @param object An object of class \code{"LeeCarter"}.
#' @inheritParams print.LeeCarter
#' @keywords internal
#' @export
summary.LeeCarter <- function(object, ...) {
  axbx <- data.frame(ax = object$coefficients$ax, 
                     bx = object$coefficients$bx,
                     row.names = object$x)
  kt <- data.frame(kt = object$coefficients$kt)
  out = structure(class = 'summary.LeeCarter', 
                  list(A = axbx, K = kt, call = object$call, info = object$info,
                       y = object$y, x_ = object$x))
  return(out)
}


#' Print summary.LeeCarter
#' @param x An object of class \code{"summary.LeeCarter"}.
#' @inheritParams print.LeeCarter
#' @keywords internal
#' @export
print.summary.LeeCarter <- function(x, ...){
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
#' @param x An object of class \code{"predict.LeeCarter"};
#' @inheritParams print_predict_default
#' @keywords internal
#' @export
print.predict.LeeCarter <- function(x, ...) {
  print_predict_default(x, ...)
  cat('k[t]-ARIMA method:', arima.string1(x$kt.arima, padding = TRUE))
  cat('\n')
}
