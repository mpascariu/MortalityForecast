# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Fri Nov 23 16:47:19 2018
# --------------------------------------------------- #


#' The Coherent Oeppen Mortality Model (Oeppen-C)
#' 
#' @inheritParams model_Oeppen
#' @param data.B A data.frame or a matrix containing mortality data for the 
#' benchmark population. Must be the same format as in \code{data}; 
#' @return The output is a list with the components:
#'  \item{input}{List with arguments provided in input. Saved for convenience;}
#'  \item{info}{Short details about the model;}
#'  \item{call}{An unevaluated function call, that is, an unevaluated 
#'  expression which consists of the named function applied to the given 
#'  arguments;}
#'  \item{coefficients}{Estimated coefficients;}
#'  \item{fitted.values}{Fitted values of the estimated model;}
#'  \item{observed.values}{The observed values used in fitting arranged in the 
#'  same format as the fitted.values;}
#'  \item{residuals}{Deviance residuals;} 
#'  \item{x}{Vector of ages used in the fitting;} 
#'  \item{y}{Vector of years used in the fitting;} 
#'  \item{benchmark}{An object of class \code{LeeCarter} containing the fitted
#'  model for the benchmark population.} 
#' @seealso 
#' \code{\link{predict.Oeppen}}
#' \code{\link{plot.Oeppen}}
#' @export
model_OeppenC <- function(data, data.B, x = NULL, y = NULL, verbose = TRUE, ...){
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
  B <- model_Oeppen(data = data.B, x = x, y = y, verbose = FALSE)
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
  fv  <- clrInv(c(kt) %*% t(bx)) + B.cdx # Inverse clr
  fv  <- sweep(fv, 2, ax, FUN = "+")
  fdx <- unclass(t(fv/rowSums(fv)))
  odx <- apply(data, 2, FUN = function(x) x/sum(x)) # observed dx - same scale
  resid <- odx - fdx
  dimnames(fdx) = dimnames(resid) = dimnames(data) <- list(x, y)
  
  # Exit
  out <- list(input = input, 
              info = info, 
              call = match.call(), 
              coefficients = cf, 
              fitted.values = fdx, 
              observed.values = odx,
              residuals = resid, 
              x = x, 
              y = y,
              benchmark = B)
  out <- structure(class = 'OeppenC', out)
  return(out)
}



#' #' Forecast the Age at Death Distribution using the Coherent Oeppen model
#' @param object An object of class \code{Oeppen};
#' @param order.B The ARIMA order for the benchmark population model;
#' @param include.drift.B Logical. Should we include a linear drift 
#' term in the ARIMA of benchmark population model?
#' @inheritParams predict.Oeppen
#' @return The output is a list with the components:
#'  \item{call}{An unevaluated function call, that is, an unevaluated 
#'  expression which consists of the named function applied to the given 
#'  arguments;}
#'  \item{predicted.values}{A list containing the predicted values together
#'  with the associated prediction intervals given by the estimated 
#'  model over the forecast horizon \code{h};}
#'  \item{conf.intervals}{Confidence intervals for the extrapolated \code{kt} 
#'  parameters;}
#'  \item{kt.arima}{An object of class \code{ARIMA} that contains all the
#'  components of the fitted time series model used in \code{kt} prediction;} 
#'  \item{kt}{The extrapolated \code{kt} parameters;}
#'  \item{x}{Vector of ages used in prediction;} 
#'  \item{y}{Vector of years used in prediction;}
#'  \item{info}{Short details about the model;}
#'  \item{benchmark}{An object containing the predicted results for the 
#'  benchmark population.}
#' @export
predict.OeppenC <- function(object,
                            h,
                            order = c(1,0,0),
                            order.B = c(0,1,0),
                            include.drift = FALSE,
                            include.drift.B = TRUE,
                            level = c(80, 95),
                            jumpchoice = c("actual", "fit"),
                            method = "ML",
                            verbose = TRUE, ...){
  
  # Benchmark Oeppen forecast
  B <- object$benchmark
  B.pred <- predict(object = B, h, order.B, include.drift.B, level, 
                    jumpchoice, method, verbose = FALSE)
  
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
  tsf <- forecast(kt.arima, h = h + 1, level = level)  # time series forecast
  fkt <- data.frame(tsf$mean, tsf$lower, tsf$upper) # forecast kt
  Cnames <- c('mean', paste0('L', level), paste0('U', level))
  dimnames(fkt) <- list(c(0, fcy), Cnames)
  
  # Get forecast d[x] based on k[t] extrapolation 
  # Here we are also adjusting for the jump-off
  J <- match.arg(jumpchoice)
  d <- get_dx_values(object = object, 
                     jumpchoice = J,
                     y = fcy, 
                     kt = fkt, 
                     B.kt = B.pred$kt)
  
  # Exit
  out <- list(call = match.call(), 
              info = object$info,
              kt = fkt, 
              kt.arima = kt.arima, 
              predicted.values = d[[1]],
              conf.intervals = d[-1], 
              x = object$x, 
              y = fcy, 
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
