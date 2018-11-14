
# ------------------------------------------
# Lee-Carter model - 2nd implementation

#' The Lee-Carter Mortality Model
#' 
#' Fit the Lee-Carter mortality model
#' @inheritParams doMortalityModels
#' @seealso 
#' \code{\link{predict.LeeCarter2}}
#' @details \insertNoCite{lee1992}{MortalityForecast}
#' @references \insertAllCited{}
#' @examples 
#' # Data
#' x  <- 0:89
#' y  <- 1985:2014
#' mx <- HMD_male$mx$GBRTENW[paste(x), paste(y)]
#' 
#' # Fit the model
#' M <- fit_LeeCarter2(data = mx, x = x, y = y)
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
fit_LeeCarter2 <- function(data, x = NULL, y = NULL, verbose = TRUE, ...){
  input <- c(as.list(environment()))
  x <- x %||% 1:nrow(data)
  y <- y %||% 1:ncol(data)
  
  # Info
  modelLN <- "Lee-Carter Mortality Model"   # long name
  modelSN <- "LC"                           # short name
  modelF  <- "log m[x,t] = a[x] + b[x]k[t] + e[x,t]" # formula
  info <- list(name = modelLN, name.short = modelSN, formula = modelF)
  
  # Estimate model parameters: a[x], b[x], k[t]
  ax  <- apply(log(data), 1, mean)
  cmx <- sweep(log(data), 1, ax, FUN = "-")
  S   <- svd(cmx)
  kt  <- S$v[,1] * sum(S$u[, 1]) * S$d[1]
  bx  <- S$u[,1] / sum(S$u[, 1])
  cf  <- list(ax = as.numeric(ax), bx = as.numeric(bx), kt = as.numeric(kt))
  
  # Variability
  var <- cumsum(S$d^2 / sum(S$d^2))[1]
  
  # Compute fitted values and devinace residuals based on the estimated model
  fv <- sweep(c(bx) %*% t(kt), 1, ax, FUN = "+") # Fitted values
  fv <- exp(fv) # fitted mx
  dimnames(fv) <- list(x, y)
  
  resid <- data - fv # residuals
  
  # Exit
  out <- list(input = input, info = info, call = match.call(), 
              fitted.values = fv, observed.values = data,
              coefficients = cf, residuals = resid, x = x, y = y)
  out <- structure(class = 'LeeCarter2', out)
  return(out)
}


#' Forecast age-specific death rates using the Lee-Carter model.
#' 
#' @param object An object of class \code{LeeCarter2}.
#' @inheritParams predict.Oeppen
#' @inherit predict.Oeppen return
#' @seealso 
#' \code{\link{fit_LeeCarter2}}
#' @author Marius D. Pascariu and Marie-Pier Bergeron-Boucher
#' @examples 
#' # Data
#' x  <- 0:89
#' y  <- 1985:2014
#' mx <- HMD_male$mx$GBRTENW[paste(x), paste(y)]
#' 
#' M <- fit_LeeCarter2(data = mx, x = x, y = y) # fit
#' P <- predict(M, h = 20)  # forecast
#' P
#' @export
predict.LeeCarter2 <- function(object, 
                               h, 
                               order = c(0, 1, 0), 
                               include.drift = TRUE,
                               level = c(80, 95),
                               jumpchoice = c("actual", "fit"),
                               method = "ML",
                               verbose = TRUE, ...){
  # Data
  mx <- t(object$input$data)
  
  # Timeline
  bop <- max(object$y) + 1
  eop <- bop + h - 1
  fcy <- bop:eop
  
  # Identify the k[t] ARIMA order
  C <- coef(object)
  ts_auto <- auto.arima(C$kt)
  order  <- order %||% arimaorder(ts_auto)
  drift  <- include.drift %||% any(names(coef(ts_auto)) %in% "drift")
  
  # Estimate/fit k[t] time-series model
  tsm <- forecast::Arima(y = C$kt, 
                         order = order, 
                         include.drift = drift, 
                         method = method)
  
  # Forecast k[t] using the time-series model
  tsf <- forecast(tsm, h = h, level = level)  # time series forecast
  fkt <- data.frame(tsf$mean, tsf$lower, tsf$upper) # forecast kt
  Cnames <- c('mean', paste0('L', level), paste0('U', level))
  dimnames(fkt) <- list(fcy, Cnames)
  
  # Get forecast m[x] based on k[t] extrapolation 
  # Here we are also adjusting for the jump-off
  fmx <- get_mx_values(object = object, 
                       kt = fkt, 
                       jumpchoice = match.arg(jumpchoice), 
                       y = fcy)
  pv <- fmx[[1]]
  CI <- fmx[-1]
  
  # Exit
  out <- list(call = match.call(), predicted.values = pv, 
              kt.arima = tsm, kt = fkt, 
              conf.intervals = CI, x = object$x, y = fcy, info = object$info)
  out <- structure(class = 'predict.LeeCarter2', out)
  return(out)
}


#' Get m[x] values based on k[t] forecast
#' In this function we compute the m[x] values based on the extrapolation of
#' the k[t] time-series. If necesary an adjustment for the jump-off is 
#' provided.
#' @inheritParams predict.LeeCarter2 
#' @inheritParams fit_LeeCarter2
#' @param kt Estimated kt vector of parameters in the Lee-Carter model;
#' @keywords internal
get_mx_values <- function(object, kt, jumpchoice, y){
  
  C  <- coef(object)
  OV <- object$observed.values
  FV <- object$fitted.values
  
  if (is.data.frame(kt)) {
    pred <- list()
    for (i in 1:ncol(kt)) {
      pred[[i]] <- get_mx_values(object, kt = kt[, i], jumpchoice, y)
    }
    names(pred) <- colnames(kt)
    return(pred)
    
  } else {
    pv  <- matrix(kt, ncol = 1) %*% C$bx
    pv  <- sweep(pv, 2, C$ax, FUN = "+")
    fmx <- t(exp(pv))
    
    if (jumpchoice == 'actual') {
      N   <- ncol(OV)
      J   <- as.numeric(OV[, N]/FV[, N]) # jump_off (%)
      fmx <- sweep(fmx, 1, J, FUN = "*")
    }
    
    dimnames(fmx) <- list(rownames(OV), y)
    
    return(fmx)
  }
}


# S3 ----------------------------------------------


#' Residuals of the Lee-Carter Mortality Model
#' @param object An object of class \code{"LeeCarter2"}
#' @inheritParams residuals_default
#' @export
residuals.LeeCarter2 <- function(object, ...){
  residuals_default(object, ...)
}


#' Print Lee-Carter model
#' @param x An object of class \code{"LeeCarter2"}
#' @inheritParams print_default
#' @keywords internal
#' @export
print.LeeCarter2 <- function(x, ...) {
  print_default(x, ...)
}


#' Summary LeeCarter2
#' @param object An object of class \code{"LeeCarter2"}.
#' @inheritParams print.LeeCarter2
#' @keywords internal
#' @export
summary.LeeCarter2 <- function(object, ...) {
  axbx <- data.frame(ax = object$coefficients$ax, 
                     bx = object$coefficients$bx,
                     row.names = object$x)
  kt <- data.frame(kt = object$coefficients$kt)
  out = structure(class = 'summary.LeeCarter2', 
                  list(A = axbx, K = kt, call = object$call, info = object$info,
                       y = object$y, x_ = object$x))
  return(out)
}


#' Print summary.LeeCarter2
#' @param x An object of class \code{"summary.LeeCarter2"}.
#' @inheritParams print.LeeCarter2
#' @keywords internal
#' @export
print.summary.LeeCarter2 <- function(x, ...){
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
#' @param x An object of class \code{"predict.LeeCarter2"};
#' @inheritParams print_predict_default
#' @keywords internal
#' @export
print.predict.LeeCarter2 <- function(x, ...) {
  print_predict_default(x, ...)
  cat('k[t]-ARIMA method:', arima.string1(x$kt.arima, padding = TRUE))
  cat('\n')
}
