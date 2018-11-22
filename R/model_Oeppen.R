# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Mon Nov 19 14:16:15 2018
# --------------------------------------------------- #

#' The Oeppen Mortality Model (Oeppen -- Codxa)
#' 
#' Fit the Oeppen model for forecasting the life table 
#' distribution of deaths. This is a Lee-Carter type model adapted to a 
#' compositional-data framework (CoDa). A key difference 
#' between the \insertCite{lee1992;textual}{MortalityForecast}
#' method and the Oeppen model is that the former fits and 
#' forecasts the death rates (mx) while the latter is based on the life table 
#' death distribution (dx). 
#' \insertCite{@See @oeppen2008 and @bergeron2017;textual}{MortalityForecast} 
#' for a detail description and mathematical formulation.
#' 
#' @inheritParams doMortalityModels 
#' @return The output is an object of class \code{Oeppen} with the components:
#'  \item{input}{List with arguments provided in input. Saved for convenience.}
#'  \item{call}{An unevaluated function call, that is, an unevaluated 
#'  expression which consists of the named function applied to the given arguments.}
#'  \item{coefficients}{Estimated coefficients.}
#'  \item{fitted.values}{Fitted values of the estimated CoDa model.}
#'  \item{residuals}{Deviance residuals.} 
#'  \item{x}{Vector of ages used in the fitting.} 
#'  \item{y}{Vector of years used in the fitting.} 
#' @seealso 
#' \code{\link{predict.Oeppen}}
#' \code{\link{plot.Oeppen}}
#' @references \insertAllCited{}
#' @author Marius D. Pascariu, Marie-Pier Bergeron-Boucher and Jim Oeppen 
#' @examples
#' # Data
#' x  <- 0:100
#' y  <- 1980:2014
#' dx <- HMD_male$dx$GBRTENW[paste(x), paste(y)]
#' 
#' # Fit model
#' M <- model_Oeppen(dx, x, y)
#' M
#' R <- residuals(M)
#' 
#' summary(M)
#' coef(M)
#' 
#' plot(M, plotType = "observed")
#' plot(M, plotType = "fitted")
#' 
#' plot(R, plotType = "scatter")
#' plot(R, plotType = "colourmap")
#' plot(R, plotType = "signplot")
#' @export
model_Oeppen <- function(data, x = NULL, y = NULL, verbose = TRUE, ...){
  input <- c(as.list(environment()))
  Oeppen.input.check(input)
  x <- x %||% 1:nrow(data)
  y <- y %||% 1:ncol(data)
  
  vsn <- sum(data)/ncol(data) * 1e-10 # very small number
  data[data == 0] <- vsn              # replace zero's with a vsn
  data <- convertFx(x, data, from = "dx", to = "dx", lx0 = 1)
  
  # Info
  modelLN <- "Compositional-Data Lee-Carter Mortality Model -- Oeppen"
  modelSN <- "Oeppen"
  modelF  <- "clr d[x,t] = a[x] + b[x]k[t]"
  info <- list(name = modelLN, name.short = modelSN, formula = modelF)
  
  # Estimate model parameters: a[x], b[x], k[t]
  dx  <- t(data)
  ax  <- apply(dx, 2, mean) # general dx pattern
  ax  <- ax/sum(ax)
  cdx <- sweep(acomp(dx), 2, ax, FUN = "-") # remove ax
  ccdx <- clr(acomp(cdx)) # Centered log ratio transform
  
  S  <- svd(ccdx) # Singular Value Decomposition of a Matrix
  kt <- S$d[1] * S$u[, 1]
  bx <- S$v[,1]
  cf <- list(ax = as.numeric(ax), bx = as.numeric(bx), kt = as.numeric(kt))
  
  # Variability
  var <- cumsum((S$d)^2/sum((S$d)^2))
  
  # Compute fitted values and devinace residuals based on the estimated model
  fv  <- clrInv(c(kt) %*% t(bx)) # Inverse clr
  fv  <- sweep(fv, 2, ax, FUN = "+")
  fdx <- unclass(t(fv/rowSums(fv)))
  odx <- apply(data, 2, FUN = function(x) x/sum(x)) # observed dx - same scale as fitted dx
  resid <- odx - fdx
  dimnames(fdx) = dimnames(resid) = dimnames(data) <- list(x, y)
  
  # Exit
  out <- list(input = input, info = info, call = match.call(), 
              fitted.values = fdx, observed.values = odx,
              coefficients = cf, residuals = resid, x = x, y = y)
  out <- structure(class = 'Oeppen', out)
  return(out)
}


#' Validate input values
#' @param X A list with input arguments provided in \code{\link{model_Oeppen}} function
#' @keywords internal
Oeppen.input.check <- function(X) {
  # Validate the other arguments
  with(X, {
    if (any(data < 0)) {
      stop("'data' contains negative values. ",
           "The compositions must always be positive or equal to zero.", 
           call. = F)
    }
    if (any(is.na(data))) {
      stop("'data' contains NA values. ",
           "The function does not know how to deal with these yet.", 
           call. = F)
    }
    if (any(is.na(data))) {
      stop("'data' contains NA values", call. = FALSE)
    }
    if (any(is.na(y))) {
      stop("'y' contains NA values", call. = FALSE)
    }
    if (any(is.na(x))) {
      stop("'x' contains NA values", call. = FALSE)
    }
    if ((!is.null(x)) & dim(data)[1] != length(x)) {
      stop("The length of 'x' is not equal to the number or rows in 'data'.", 
           call. = FALSE)
    }
    if ((!is.null(y)) & dim(data)[2] != length(y)) {
      stop("The length of 'y' is not equal to the number or columns in 'data'.", 
           call. = FALSE)
    }
  })
}



#' Forecast the Age at Death Distribution using the Oeppen model.
#' 
#' @param object An object of class \code{Oeppen}.
#' @param order A specification of the non-seasonal part of the ARIMA model: 
#'  the three components (p, d, q) are the AR order, the degree of differencing, 
#'  and the MA order. If \code{order = NULL}, the ARIMA order will be estimated 
#'  automatically using the KPPS algorithm.
#' @param include.drift Logical. Should the ARIMA model include a linear drift term?
#'  If \code{include.drift = NULL}, the model will be estimated automatically.
#' @param method ARIMA fitting method: maximum likelihood or minimize conditional 
#'  sum-of-squares. Options to use: conditional-sum-of-squares (\code{"CSS-ML"}), 
#'  maximum likelihood (\code{"ML"}) and \code{"CSS"}.
#' @param ... Additional arguments to be passed to \code{\link[forecast]{Arima}}
#' @inheritParams doForecasts
#' @return The output is an object of class \code{"predict.Oeppen"} with the components:
#'  \item{call}{An unevaluated function call, that is, an unevaluated 
#'  expression which consists of the named function applied to the given arguments;}
#'  \item{predicted.values}{A list containing the predicted values together
#'  with the associated prediction intervals given by the estimated 
#'  model over the forecast horizon \code{h};}
#'  \item{kt.arima}{An object of class \code{ARIMA} that contains all the
#'  components of the fitted time series model used in \code{kt} prediction;} 
#'  \item{kt}{The extrapolated \code{kt} parameters;}
#'  \item{conf.intervals}{Confidence intervals for the extrapolated \code{kt} 
#'  parameters;}
#'  \item{x}{Vector of ages used in prediction;} 
#'  \item{y}{Vector of years used in prediction;}
#'  \item{info}{Short details about the model.}
#' @author Marius D. Pascariu, Marie-Pier Bergeron-Boucher and Jim Oeppen 
#' @examples 
#' # Example 1 ----------------------
#' x  <- 0:100
#' y  <- 1980:2014
#' dx <- HMD_male$dx$GBRTENW[paste(x), paste(y)]
#' M  <- model_Oeppen(dx, x, y)
#' P  <- predict(M, h = 16)
#' 
#' plot(P, plotType = "mean")
#' plot(P, plotType = "lower")
#' plot(P, plotType = "upper")
#' 
#' #' # Example 2 ----------------------
#' # One can specify manually the ARIMA order, a drift to be included or not 
#' # and the jump choice of the first forecast year.
#' P2 <- predict(M, h = 20, order = c(0,1,0), include.drift = TRUE, jumpchoice = "fit")
#' 
#' \dontrun{
#' # Example 3 ----------------------
#' # Compute life tables using forecast values using the MortalityLaws R package
#' library(MortalityLaws)
#' dx <- P$predicted.values
#' lt <- LifeTable(x = P$x, dx = dx)
#' }
#' @export
predict.Oeppen <- function(object,
                           h, 
                           order = NULL, 
                           include.drift = NULL,
                           level = c(80, 95), 
                           jumpchoice = c("actual", "fit"), 
                           method = "ML", 
                           verbose = TRUE, ...){
  jumpchoice <- match.arg(jumpchoice)
  
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
                       adj = 1)
  
  pv <- fdx[[1]]
  CI <- fdx[-1]
  names(CI) <- Cnames[-1]
  
  # Exit
  out <- list(call = match.call(), predicted.values = pv,
              kt.arima = kt.arima, kt = fkt, 
              conf.intervals = CI, x = object$x, y = fcy, info = object$info)
  out <- structure(class = 'predict.Oeppen', out)
  return(out)
}


#' #' Get d[x] values based on k[t] forecast
#' @inheritParams get_mx_values
#' @keywords internal
get_dx_values <- function(object, kt, jumpchoice, y, adj = 1) {
  
  C  <- coef(object)
  OV <- t(object$observed.values)
  FV <- t(fitted(object))
  
  if (is.data.frame(kt)) {
    pred <- list()
    for (i in 1:ncol(kt)) {
      pred[[i]] <- get_dx_values(object, kt = kt[, i], jumpchoice, y)
    }
    names(pred) <- colnames(kt)
    return(pred)
    
  } else {
    pv  <- clrInv(c(kt) %*% t(C$bx)) + adj
    pdx <- sweep(pv, 2, C$ax, FUN = "+") # predicted dx values
    pdx <- unclass(pdx/rowSums(pdx))
    
    if (jumpchoice == 'actual') {
      N  <- nrow(OV)
      jump_off <- as.numeric(OV[N, ]/FV[N, ])
      pdx <- sweep(pdx, 2, jump_off, FUN = "*")
    }
    
    out <- unclass(t(pdx/rowSums(pdx)))
    dimnames(out) <- list(colnames(OV), y)
    return(out)
  }
} 


# S3 ----------------------------------------------

#' Residuals of the Oeppen Mortality Model
#' @param object An object of class \code{"Oeppen"}
#' @inheritParams residuals_default
#' @export
residuals.Oeppen <- function(object, ...){
  residuals_default(object, ...)
}


#' Print Oeppen
#' @param x An object of class \code{"Oeppen"}
#' @inheritParams print_default
#' @keywords internal
#' @export
print.Oeppen <- function(x, ...) {
  print_default(x, ...)
}


#' Summary Oeppen
#' @param object An object of class \code{"Oeppen"}.
#' @inheritParams print.Oeppen
#' @keywords internal
#' @export
summary.Oeppen <- function(object, ...) {
  axbx <- data.frame(ax = object$coefficients$ax, 
                     bx = object$coefficients$bx,
                     row.names = object$x)
  kt <- data.frame(kt = object$coefficients$kt)
  out = structure(class = 'summary.Oeppen', 
                  list(A = axbx, K = kt, call = object$call, info = object$info,
                       y = object$y, x_ = object$x))
  return(out)
}


#' Print summary.Oeppen
#' @param x An object of class \code{"summary.Oeppen"}.
#' @inheritParams print.Oeppen
#' @keywords internal
#' @export
print.summary.Oeppen <- function(x, ...){
  cat('\nFit  :', x$info$name)
  cat('\nModel:', x$info$formula)
  cat('\n\nCoefficients:\n')
  A <- head_tail(x$A, digits = 5, hlength = 6, tlength = 6)
  K <- head_tail(data.frame(. = '|', y = as.integer(x$y), kt = x$K),
                 digits = 5, hlength = 6, tlength = 6)
  print(data.frame(A, K))
  cat('\n')
}


#' Print predict.Oeppen
#' @param x An object of class \code{"predict.Oeppen"};
#' @inheritParams print.Oeppen
#' @keywords internal
#' @export
print.predict.Oeppen <- function(x, ...) {
  print_predict_default(x, ...)
  cat('k[t]-ARIMA method:', arima.string1(x$kt.arima, padding = TRUE))
  cat('\n')
}


#' ggplot the observed and fitted values of an Oeppen model
#' 
#' @inherit plot.MEM details
#' @inheritParams plot.MEM
#' @examples 
#' # For examples go to ?model_Oeppen
#' @export
plot.Oeppen <- function(x, plotType = c("fitted", "observed"), 
                           ny = 7, level = 80, ...){
  plot.MEM(x, plotType, ny, level, ...)
}


#' ggplot the predicted values of a CoDa-LC mortality model
#' 
#' @param x An object of the class \code{\link{predict.Oeppen}}.
#' @param plotType The type of the plot. The alternatives are 
#' \code{"mean", "lower", "upper"}. Default: \code{"mean"}.
#' @inheritParams plot.predict.MEM
#' @examples 
#' # For examples go to ?predict.Oeppen
#' @export
plot.predict.Oeppen <- function(x, plotType = c("mean", "lower", "upper"), 
                                   ny = 7, level = 80, ...) 
{
  plotType <- match.arg(plotType)
  if (plotType == "mean") {
    mat <- x$predicted.values
    P <- ggplotDistribConvergence(mat, x = x$x, y = x$y, ny, level) + 
      labs(subtitle = "Forecast Values - Best estimate")
    
  } else if (plotType == "lower") {
    mat <- x$conf.intervals[[1]]
    P <- ggplotDistribConvergence(mat, x = x$x, y = x$y, ny, level) + 
      labs(subtitle = "Forecast Values - lower bound")
    
  } else if (plotType == "upper") {
    mat <- x$conf.intervals[[length(x$conf.intervals)/2 + 1]]
    P <- ggplotDistribConvergence(mat, x = x$x, y = x$y, ny, level) + 
      labs(subtitle = "Forecast Values - upper bound")
    
  } 
  suppressMessages(print(P))
}

