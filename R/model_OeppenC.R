# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Thu Nov 22 13:31:36 2018
# --------------------------------------------------- #


#' The Coherent Oeppen Mortality Model (Oeppen-C)
#' 
#' @inheritParams model_Oeppen
#' @param B.data A list containing multiple matrices with mortality data 
#' to be used as benchmark. Format: ages \code{x} as row and time \code{y} 
#' as column.
#' @export
model_OeppenC <- function(data, B.data, x = NULL, y = NULL, verbose = TRUE, ...){
  input <- c(as.list(environment()))
  Oeppen.input.check(input)
  x <- x %||% 1:nrow(data)
  y <- y %||% 1:ncol(data)
  
  vsn <- sum(data)/ncol(data) * 1e-10 # very small number
  data[data == 0] <- vsn              # replace zero's with a vsn
  data <- convertFx(x, data, from = "dx", to = "dx", lx0 = 1)
  
  # Info
  modelLN <- "Coherent Oeppen Mortality Model"
  modelSN <- "Oeppen-C"
  modelF  <- "clr d[x,t] = a[x] + b[x]k[t]"
  info <- list(name = modelLN, name.short = modelSN, formula = modelF)
  
  # Fit benchmark model
  B <- model_Oeppen(data = B.data, x = x, y = y, verbose = FALSE)
  B.cdx <- with(coef(B), clrInv(c(kt) %*% t(bx)))
  
  # Estimate model parameters: a[x], b[x], k[t]
  dx  <- t(data)
  ax  <- apply(dx, 2, mean) # general dx pattern
  ax  <- ax/sum(ax)
  A.cdx <- sweep(acomp(dx), 2, ax, FUN = "-") # remove ax
  cdx   <- A.cdx - B.cdx
  ccdx  <- clr(acomp(cdx)) # Centered log ratio transform
  
  S  <- svd(ccdx) # Singular Value Decomposition of a Matrix
  kt <- S$d[1] * S$u[, 1]
  bx <- S$v[,1]
  cf <- list(ax = as.numeric(ax), bx = as.numeric(bx), kt = as.numeric(kt))
  
  # Variability
  var <- cumsum((S$d)^2/sum((S$d)^2))
  
  # Compute fitted values and devinace residuals based on the estimated model
  fv  <- clrInv(c(kt) %*% t(bx)) # Inverse clr
  fv  <- fv + B.cdx 
  fv  <- sweep(fv, 2, ax, FUN = "+")
  fdx <- unclass(t(fv/rowSums(fv)))
  odx <- apply(data, 2, FUN = function(x) x/sum(x)) # observed dx - same scale as fitted dx
  resid <- odx - fdx
  dimnames(fdx) = dimnames(resid) = dimnames(data) <- list(x, y)
  
  # Exit
  out <- list(input = input, info = info, call = match.call(), 
              fitted.values = fdx, observed.values = odx,
              coefficients = cf, residuals = resid, x = x, y = y,
              B.model = B)
  out <- structure(class = 'OeppenC', out)
  return(out)
}



#' #' Forecast the Age at Death Distribution using the Coherent Oeppen model
#' @param object An object of class \code{Oeppen}.
#' @inheritParams predict.Oeppen
#' @inherit predict.Oeppen return
#' @export
predict.OeppenC <- function(object,
                            h,
                            order = NULL,
                            include.drift = NULL,
                            level = c(80, 95),
                            jumpchoice = c("actual", "fit"),
                            method = "ML",
                            verbose = TRUE, ...){
  jumpchoice <- match.arg(jumpchoice)
  
  # Benchmark Oeppen forecast
  B <- object$B.model
  B.pred <- predict(object = B, h, order, include.drift, level, 
                    jumpchoice, method, verbose = FALSE)
  B.cdx <- clrInv(c(B.pred$kt$mean) %*% t(coef(B)$bx))
  
  # Timeline
  bop <- max(object$y) + 1
  eop <- bop + h - 1
  fcy <- bop:eop
  
  # Identify the k[t] ARIMA order
  C <- coef(object)
  A <- find_arima(C$kt)
  
  # forecast kt; ax and bx are time independent.
  kt.arima <- forecast::Arima(y = C$kt, 
                              order = order %||% A$order, 
                              include.drift = include.drift %||% A$drift,
                              method = method)
  
  # Forecast k[t] using the time-series model
  tsf <- forecast(kt.arima, h = h, level = level)  # time series forecast
  fkt <- data.frame(tsf$mean, tsf$lower, tsf$upper) # forecast kt
  Cnames <- c('mean', paste0('L', level), paste0('U', level))
  colnames(fkt) <- Cnames
  
  # Get forecast d[x] based on k[t] extrapolation 
  # Here we are also adjusting for the jump-off
  fdx <- get_dx_values(object = object, 
                       kt = fkt, 
                       y = fcy, 
                       jumpchoice = jumpchoice,
                       adj = B.cdx)
  
  pv <- fdx
  # pv <- fdx[[1]]
  # CI <- fdx[-1]
  # names(CI) <- Cnames[-1]
  
  # Exit
  out <- list(call = match.call(), predicted.values = pv,
              kt.arima = kt.arima, kt = fkt, 
              conf.intervals = NULL, x = object$x, y = fcy, info = object$info,
              benchmark = B.pred)
  out <- structure(class = 'predict.OeppenC', out)
  return(out)
}


# S3 ----------------------------------------------
#' Residuals of the Coherent Oeppen Mortality Model
#' @param object An object of class \code{"OeppenC"}
#' @inheritParams residuals_default
#' @export
residuals.OeppenC <- function(object, ...){
  residuals_default(object, ...)
}


#' Print Coherent Oeppen
#' @param x An object of class \code{"OeppenC"}
#' @inheritParams print_default
#' @keywords internal
#' @export
print.OeppenC <- function(x, ...) {
  print_default(x, ...)
}


#' Summary Coherent Oeppen
#' @param object An object of class \code{"OeppenC"}.
#' @inheritParams print.Oeppen
#' @keywords internal
#' @export
summary.OeppenC <- function(object, ...) {
  axbx <- data.frame(ax = object$coefficients$ax, 
                     bx = object$coefficients$bx,
                     row.names = object$x)
  kt <- data.frame(kt = object$coefficients$kt)
  out = structure(class = 'summary.OeppenC', 
                  list(A = axbx, K = kt, call = object$call, info = object$info,
                       y = object$y, x_ = object$x))
  return(out)
}


#' Print summary.OeppenC
#' @param x An object of class \code{"summary.OeppenC"}.
#' @inheritParams print.OeppenC
#' @keywords internal
#' @export
print.summary.OeppenC <- function(x, ...){
  cat('\nFit  :', x$info$name)
  cat('\nModel:', x$info$formula)
  cat('\n\nCoefficients:\n')
  A <- head_tail(x$A, digits = 5, hlength = 6, tlength = 6)
  K <- head_tail(data.frame(. = '|', y = as.integer(x$y), kt = x$K),
                 digits = 5, hlength = 6, tlength = 6)
  print(data.frame(A, K))
  cat('\n')
}


#' Print predict.OeppenC
#' @param x An object of class \code{"predict.OeppenC"};
#' @inheritParams print.OeppenC
#' @keywords internal
#' @export
print.predict.OeppenC <- function(x, ...) {
  print_predict_default(x, ...)
  cat('k[t]-ARIMA method:', arima.string1(x$kt.arima, padding = TRUE))
  cat('\n')
}
