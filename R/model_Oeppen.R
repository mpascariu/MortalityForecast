
#' The Compositional-Data Lee-Carter Mortality Model (CoDa-LC)
#' 
#' Fit Compositional Data model for forecasting the life table 
#' distribution of deaths. CoDa is a Lee-Carter type model. A key difference 
#' between the \insertCite{lee1992;textual}{MortalityForecast}
#' method and the Compositional Data model (CoDa-LC) is that the former fits and 
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
#' dx <- MortalityForecast.data$dx[paste(x), paste(y)]
#' 
#' # Fit CoDa-LC model
#' M <- fit_Oeppen(dx, x, y)
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
fit_Oeppen <- function(data, x = NULL, y = NULL, verbose = TRUE, ...){
  input <- c(as.list(environment()))
  Oeppen.input.check(input)
  x <- x %||% 1:nrow(data)
  y <- y %||% 1:ncol(data)
  
  vsn <- sum(data)/ncol(data) * 1e-10 # very small number
  data[data == 0] <- vsn              # replace zero's with a vsn
  data <- convertFx(x, data, from = "dx", to = "dx", lx0 = 1)
  
  close.dx  <- unclass(acomp(t(data)))      # data close
  ax        <- geometricmeanCol(close.dx) # geometric mean
  close.ax  <- ax/sum(ax)
  dxc       <- sweep(close.dx, 2, close.ax, "/") # centering
  close.dxc <- dxc/rowSums(dxc)
  clr_dxc   <- clr(close.dxc) # clr
  
  # SVD: bx and kt
  par <- svd(clr_dxc, nu = 1, nv = 1)
  U   <- par$u
  V   <- par$v
  S   <- diag(par$d)
  bx  <- V[, 1]
  kt  <- S[1, 1] * U[, 1]
  
  var <- cumsum((par$d)^2/sum((par$d)^2)) # variability
  cf  <- list(ax = as.numeric(close.ax), bx = as.numeric(bx), kt = as.numeric(kt))
  clr.proj.fit <- matrix(kt, ncol = 1) %*% bx # projections
  BK.proj.fit  <- unclass(clrInv(clr.proj.fit)) # Inv clr
  proj.fit     <- sweep(BK.proj.fit, 2, close.ax, FUN = "*") # add geometric mean
  fD    <- t(proj.fit/rowSums(proj.fit))
  oD    <- apply(data, 2, FUN = function(x) x/sum(x)) # observed dx - same scale as fitted dx
  resid <- oD - fD
  dimnames(fD) = dimnames(resid) = dimnames(data) <- list(x, y)
  
  modelLN <- "Compositional-Data Lee-Carter Mortality Model"
  modelSN <- "CoDa-LC"
  modelF <- "clr d[x,t] = a[x] + b[x]k[t]"
  info <- list(name = modelLN, name.short = modelSN, formula = modelF)
  
  out <- list(input = input, info = info, call = match.call(), 
              fitted.values = fD, observed.values = oD,
              coefficients = cf, residuals = resid, x = x, y = y)
  out <- structure(class = 'Oeppen', out)
  return(out)
}


#' Validate input values
#' @param X A list with input arguments provided in \code{\link{fit_Oeppen}} function
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
      stop("'data' contains NA values", call. = F)
    }
    if (any(is.na(y))) {
      stop("'y' contains NA values", call. = F)
    }
    if (any(is.na(x))) {
      stop("'x' contains NA values", call. = F)
    }
    if ((!is.null(x)) & dim(data)[1] != length(x)) {
      stop("The length of 'x' is not equal to the number or rows in 'data'.", call. = F)
    }
    if ((!is.null(y)) & dim(data)[2] != length(y)) {
      stop("The length of 'y' is not equal to the number or columns in 'data'.", call. = F)
    }
  })
}


# S3 ----------------------------------------------
#' Residuals of the CoDa-LC Mortality Model
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




#' Forecast the Age at Death Distribution using CoDa-LC model.
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
#' dx <- MortalityForecast.data$dx[paste(x), paste(y)]
#' M  <- fit_Oeppen(dx, x, y)
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
#' dx <- P$predicted.values$mean
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
  
  dx  <- t(object$input$data)
  bop <- max(object$y) + 1
  eop <- bop + h - 1
  fcy <- bop:eop
  jc  <- jumpchoice[1]
  cf  <- coef(object)
  kt  <- cf$kt
  
  # forecast kt; ax and bx are time independent.
  ts_auto = auto.arima(kt)
  
  AO  <- order %||% arimaorder(ts_auto)
  ID  <- include.drift %||% any(names(coef(ts_auto)) %in% "drift")
  tsm <- forecast::Arima(y = kt, order = AO, include.drift = ID, method = method, ...)
  tsf <- forecast(tsm, h = h, level = level)  # time series forecast
  
  fkt <- data.frame(tsf$mean, tsf$lower, tsf$upper) # forecast kt
  Cnames <- c('mean', paste0('L', level), paste0('U', level))
  colnames(fkt) <- Cnames
  
  fdx <- compute_dx(dx = dx, kt = fkt, ax = cf$ax, bx = cf$bx, # forecast dx
                    fit = t(fitted(object)), y = fcy, jumpchoice = jc)
  pv <- fdx[[1]]
  CI <- fdx[-1]
  names(CI) <- Cnames[-1]
  
  out <- list(call = match.call(), predicted.values = pv,
              kt.arima = tsm, kt = fkt, 
              conf.intervals = CI, x = object$x, y = fcy, info = object$info)
  out <- structure(class = 'predict.Oeppen', out)
  return(out)
}


#' Internal function
#' @inheritParams fit_Oeppen
#' @inheritParams predict.Oeppen
#' @param kt Estimated kt vector of parameters
#' @param ax Estimated ax vector of parameters
#' @param bx Estimated bx vector of parameters
#' @param fit Fitted values
#' @keywords internal
compute_dx <- function(dx, kt, ax, bx, fit, y, jumpchoice) {
  
  if (is.data.frame(kt)) {
    pred <- list()
    for (i in 1:ncol(kt)) {
      pred[[i]] <- compute_dx(dx, kt = kt[, i], ax, bx, fit, y, jumpchoice)
      colnames(pred[[i]]) <- y
    }
    return(pred)
    
  } else {
    dx_nrow  <- nrow(dx)
    close.dx <- acomp(dx)
    jump_off <- as.numeric(close.dx[dx_nrow, ]/fit[dx_nrow, ])
    clr_proj <- matrix(kt, ncol = 1) %*% bx
    bk_      <- unclass(clrInv(clr_proj))
    dx_      <- sweep(bk_, 2, ax, FUN = "*")
    if (jumpchoice == 'actual') dx_ <- sweep(dx_, 2, jump_off, FUN = "*")
    out <- t(dx_/rowSums(dx_))
    rownames(out) <- colnames(dx)
    return(out)
  }
} 


# S3 ----------------------------------------------


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


#' Identify ARIMA model - internal function
#' @param object An object generate by Arima function
#' @param padding Logical.
#' @keywords internal
arima.string1 <- function(object, padding = FALSE) {
  order  <- object$arma[c(1, 6, 2, 3, 7, 4, 5)]
  nc     <- names(coef(object))
  result <- paste0("ARIMA(", order[1], ",", order[2], ",", order[3], ")")
  
  if (order[7] > 1 & sum(order[4:6]) > 0) 
    result <- paste0(result, "(", order[4], ",", order[5], 
                     ",", order[6], ")[", order[7], "]")
  if (!is.null(object$xreg)) {
    if (NCOL(object$xreg) == 1 & is.element("drift", nc)) 
      result <- paste(result, "with drift        ")
    else result <- paste("Regression with", result, "errors")
  }
  else {
    if (is.element("constant", nc) | is.element("intercept", nc)) 
      result <- paste(result, "with non-zero mean")
    else if (order[2] == 0 & order[5] == 0) 
      result <- paste(result, "with zero mean    ")
    else result <- paste(result, "                  ")
  }
  if (!padding) 
    result <- gsub("[ ]*$", "", result)
  return(result)
}


#' ggplot the observed and fitted values of a CoDa-LC mortality model
#' 
#' @inherit plot.MEM details
#' @inheritParams plot.MEM
#' @examples 
#' # For examples go to ?fit_Oeppen
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

