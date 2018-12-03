# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Mon Dec  3 14:49:36 2018
# --------------------------------------------------- #


#' The Coherent Oeppen Mortality Model (Oeppen-C)
#' 
#' @inheritParams model.Oeppen
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
#' @examples 
#' # Example 1 ----------------------
#' # Data
#' x  <- 0:110
#' y  <- 1980:2016
#' dx_M <- HMD_male$dx$USA[paste(x), paste(y)]
#' dx_F <- HMD_female$dx$USA[paste(x), paste(y)]
#' 
#' # Replace zeros 
#' dx_M <- replace.zeros(dx_M)
#' dx_F <- replace.zeros(dx_F)
#' 
#' # Fit model for US males using US females as benchmark
#' M <- model.OeppenC(data = dx_M,
#'                    data.B = dx_F,
#'                    x = x, y = y)
#' M
#' 
#' summary(M)
#' coef(M)
#' coef(M$benchmark)
#' 
#' # Plot observed and fitted values
#' plot(M, plotType = "observed")
#' plot(M, plotType = "fitted")
#' 
#' # Plot residuals
#' R <- residuals(M)
#' plot(R, plotType = "scatter")
#' plot(R, plotType = "colourmap")
#' plot(R, plotType = "signplot")
#' 
#' # Perform forecasts
#' P  <- predict(M, h = 16)
#' P
#' 
#' plot(P, plotType = "mean")
#' plot(P, plotType = "lower")
#' plot(P, plotType = "upper")
#' @seealso 
#' \code{\link{predict.Oeppen}}
#' \code{\link{plot.Oeppen}}
#' @export
model.OeppenC <- function(data, 
                          data.B, 
                          x = NULL, 
                          y = NULL, 
                          verbose = TRUE, 
                          ...){
  
  input <- c(as.list(environment()))
  Oeppen.input.check(input)
  x <- x %||% 1:nrow(data)
  y <- y %||% 1:ncol(data)
  
  data <- convertFx(x, data, from = "dx", to = "dx", lx0 = 1)
  
  # Info
  modelLN <- "Coherent Oeppen Mortality Model"
  modelSN <- "Oeppen-C"
  modelF  <- "clr d[x,t] = {a[x] + b[x]k[t]} + {A[x] + B[x]K[t]}"
  info <- list(name = modelLN, name.short = modelSN, formula = modelF)
  
  # Fit benchmark model
  B <- model.Oeppen(data = data.B, x = x, y = y, verbose = FALSE)
  B.cdx <- with(coef(B), clrInv(c(kt) %*% t(bx)))
  
  # Estimate model parameters: a[x], b[x], k[t]
  dx    <- data %>% t %>% acomp %>% unclass # data close
  ax    <- geometricmeanCol(dx) # geometric mean # general dx pattern
  ax    <- ax/sum(ax)
  A.cdx <- sweep(dx, 2, ax, FUN = "/") # remove ax
  cdx   <- A.cdx - B.cdx
  cdx   <- cdx/rowSums(cdx)
  ccdx  <- clr(cdx) # Centered log ratio transform
  
  S  <- svd(ccdx) # Singular Value Decomposition of a Matrix
  kt <- S$d[1] * S$u[, 1]
  bx <- S$v[,1]
  cf <- list(ax = as.numeric(ax), bx = as.numeric(bx), kt = as.numeric(kt))
  
  # Variability
  var <- cumsum((S$d)^2/sum((S$d)^2))
  
  # Compute fitted values and devinace residuals based on the estimated model
  fv  <- clrInv(c(kt) %*% t(bx)) + B.cdx # Inverse clr
  fv  <- sweep(unclass(fv), 2, ax, FUN = "*")
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



#' Forecast the age-at-death distribution using the Coherent Oeppen model
#' 
#' @param object An object of class \code{Oeppen};
#' @param order.B The ARIMA order for the benchmark population model. This is
#' A specification of the non-seasonal part of the ARIMA model: the three 
#' components (p, d, q) are the AR order, the degree of differencing, and 
#' the MA order. If \code{order.B = NULL}, the ARIMA order will be estimated 
#' automatically using the KPPS algorithm;
#' @param include.drift.B Logical. Should we include a linear drift 
#' term in the ARIMA of benchmark population model?
#' @param order.D The ARIMA order driving the deviation from the benchmark 
#' population model. If \code{order = NULL}, this will be estimated 
#' automatically.
#' @param include.drift.D Logical. Should we include a linear drift 
#' term in the ARIMA driving the deviation from the benchmark?
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
#' @examples # For examples go to ?model.OeppenC
#' @export
predict.OeppenC <- function(object,
                            h,
                            order.B = c(0,1,0),
                            include.drift.B = TRUE,
                            order.D = c(1,0,0),
                            include.drift.D = FALSE,
                            level = c(80, 95),
                            jumpchoice = c("actual", "fit"),
                            method = "ML",
                            verbose = TRUE, 
                            ...){
  
  # Benchmark Oeppen forecast
  B <- object$benchmark
  B.pred <- predict(object = B, 
                    h = h, 
                    order = order.B, 
                    include.drift = include.drift.B, 
                    level = level, 
                    jumpchoice = jumpchoice, 
                    method = method, 
                    verbose = FALSE)
  
  # Timeline
  bop <- max(object$y) + 1
  eop <- bop + h - 1
  fcy <- bop:eop
  
  # Identify the k[t] ARIMA order
  C <- coef(object)
  A <- find_arima(C$kt)
  
  # forecast kt; ax and bx are time independent.
  kt.arima <- forecast::Arima(y = C$kt, 
                              order = order.D %||% A$order, 
                              include.drift = include.drift.D %||% A$drift,
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

#' @rdname residuals.Oeppen
#' @export
residuals.OeppenC <- function(object, ...){
  residuals_default(object, ...)
}


#' @rdname print_default
#' @export
print.OeppenC <- function(x, ...) {
  print_default(x, ...)
}


#' @rdname summary.Oeppen
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


#' @rdname print_default
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


#' @rdname print_default
#' @export
print.predict.OeppenC <- function(x, ...) {
  print_predict_default(x, ...)
  cat('k[t]- Benchmark ARIMA:', 
      arima.string1(x$benchmark$kt.arima, padding = TRUE), '\n')
  cat('k[t]- Deviation ARIMA:', arima.string1(x$kt.arima, padding = TRUE))
  cat('\n')
}




