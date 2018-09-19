
#' The Maximum-Entropy Mortality Model 
#' 
#' @param n The maximum order of the moments to be used.
#' @inheritParams doMortalityModels
#' @return The output is an object of class \code{MEM} with 
#' the components:
#'  \item{input}{List with arguments provided in input. Saved for convenience.}
#'  \item{call}{An unevaluated function call, that is, an unevaluated 
#'  expression which consists of the named function applied to the given arguments;}
#'  \item{coefficients}{Estimated coefficients;}
#'  \item{fitted.values}{Fitted values of the estimated model;}
#'  \item{observed.values}{Observed values used in fitting;}
#'  \item{fitted.raw.moments}{Fitted raw moments;}
#'  \item{observed.raw.moments}{Observed raw moments of the input data;}
#'  \item{residuals}{Deviance residuals;} 
#'  \item{VAR}{Object containing the components of the fitted time 
#'  series model to the extrapolate moments;}
#'  \item{x}{Vector of ages used in the fitting;} 
#'  \item{y}{Vector of years used in the fitting.} 
#' @seealso 
#' \code{\link{predict.MEM}} 
#' \code{\link{plot.MEM}} 
#' \code{\link{findMoments}}
#' @examples 
#' # Data
#' x  <- 0:110
#' y  <- 1985:2014
#' dx <- MortalityForecast.data$dx[paste(x), paste(y)]
#' 
#' # Fit model the Maximum-Entropy Mortality of order 5,
#' # that is using the first 6 statistical moments (moment 0 to 5).
#' 
#' M <- fit_MEM(data = dx, x = x, y = y, n = 5)
#' M
#' R <- residuals(M)
#' 
#' plot(M, plotType = "observed")
#' plot(M, plotType = "fitted")
#' 
#' plot(R, plotType = "scatter")
#' plot(R, plotType = "colourmap")
#' plot(R, plotType = "signplot")
#' @export
fit_MEM <- function(data, x = NULL, y = NULL, n = 5, verbose = FALSE, ...) {
  input <- c(as.list(environment()))

  AY  <- find_ages_and_years(data, x, y)
  x   <- AY$x
  y   <- AY$y
  M   <- findMoments(data, x, y, n)
  orM <- M$raw.moments
  nM  <- M$normalized.moments
  nMT <- log(abs(nM))
  sg  <- sign(nM)
  
  V   <- fit_MRW(t(nMT), x = NULL, y, include.drift = TRUE)
  fnM <- fnM <- t(exp(fitted(V)) * as.numeric(sg[1,]))       # fitted normalized moments
  frM <- convertMoments(fnM, from = "normalized", to = "raw") # fitted raw moments
  
  fD  <- findDensity(frM[-1, ], x)$density          # fitted dx
  fD  <- cbind(NA, fD) 
  oD  <- apply(data, 2, FUN = function(x) x/sum(x)) # observed dx - same scale as fitted dx
  res <- oD - fD                                    # residuals
  dimnames(oD) = dimnames(fD) = dimnames(res) = list(x = x, y = y)
  
  modelLN <- "Maximum-Entropy Mortality Model "
  modelSN <- "MEM"
  modelF <- "MaxEnt M[n,t] = ..."
  info <- list(name = modelLN, name.short = modelSN, formula = modelF)
  
  out <- list(input = input, info = info, call = match.call(),
              fitted.values = fD, observed.values = oD, coefficients = coef(V),
              residuals = res, fitted.raw.moments = frM, observed.raw.moments = orM, 
              VAR = V, x = x, y = y)
  out <- structure(class = 'MEM', out)
  return(out)
}


#' Find age and year vectors
#' @inheritParams fit_MEM
#' @keywords internal
#' @export
find_ages_and_years <- function(data, x, y) {
  if (!is.null(x) & nrow(data) != length(x)) {
    stop("'data' and 'x' must have the same dimension!\n", 
         "Check: nrow(data) == length(x)", call. = FALSE)
  } 
  if (!is.null(y) & ncol(data) != length(y)) {
    stop("'data' and 'y' must have the same dimension!\n", 
         "Check: ncol(data) == length(y)", call. = FALSE)
  } 
  
  dn  <- dimnames(data)
  dn1 <- suppressWarnings(as.numeric(dn[[1]]))
  dn2 <- suppressWarnings(as.numeric(dn[[2]]))
  if (is.null(x) & !is.null(dn[1]) & !all(is.na(dn1))) {
    x <- dn1 
  } else {
    if (is.null(x) & is.null(dn[1])) {
      x <- 1:nrow(data) 
    } 
  }
  
  if (is.null(y) & !is.null(dn[2]) & !all(is.na(dn2))) {
    y = dn2 } else {
      if (is.null(y) & is.null(dn[2])) {
        y <- 1:ncol(data) 
      } 
    }
  
  out <- list(x = x, y = y)
  return(out)
}

# S3 ----------------------------------------------

#' Residuals of the Maximum-Entropy Mortality Model 
#' 
#' Computed deviance residuals for a fitted Maximum-Entropy Mortality Model.
#' @param object An object of class \code{MEM}.
#' @param ... Further arguments passed to or from other methods.
#' @examples 
#' x  <- 0:110
#' y  <- 1965:2014
#' dx <- MortalityForecast.data$dx[paste(x), paste(y)]
#' M  <- fit_MEM(dx, x, y, n = 5)
#' residuals(M)
#' @export
residuals.MEM <- function(object, ...) {
  structure(class = "residMF", as.matrix(object$residuals))
}


#' Print function for Maximum-Entropy Mortality Model 
#' @param x An object of class \code{MEM}.
#' @inheritParams residuals.MEM
#' @keywords internal
#' @export
print.MEM <- function(x, ...) {
  cat('\nFit  :', x$info$name)
  cat('\nModel:', x$info$formula)
  cat('\nCall : '); print(x$call)
  cat('\nAges in fit   : ', paste0(range(x$x), collapse = ' - '))
  cat('\nYears in fit  : ', paste0(range(x$y), collapse = ' - '))
  cat('\nMoments in fit: ', paste0("M", c(0, x$input$n), collapse = ' - '))
  cat('\n')
}



#' Predict Maximum-Entropy Mortality Model
#' 
#' @param object An object of class \code{\link{fit_MEM}}.
#' @param x.h Numerical vector indicating the ages to be considered in 
#' reconstruction of the density over the forecast horizon. 
#' If \code{NULL}, the number of estimated data points is equal to the number 
#' of fitted data points. This argument can be used for example to estimate a 
#' density between 0 and 130 given the fact that the model was fitted on a 
#' dataset containing values for 0-100 only.
#' @inheritParams doForecasts
#' @seealso \code{\link{fit_MEM}}
#' @examples
#' x  <- 0:110
#' y  <- 1985:2016
#' dx <- MortalityForecast.data$dx[paste(x), paste(y)]
#' M  <- fit_MEM(dx, x, y, n = 5)
#' P  <- predict(M, h = 16, x.h = 0:120)
#' 
#' plot(P, plotType = "mean")
#' plot(P, plotType = "lower")
#' plot(P, plotType = "upper")
#' 
#' plot(P, M, plotType = "raw_moments")
#' plot(P, M, plotType = "normalised_moments")
#' @export
predict.MEM <- function(object, h, x.h = NULL, level = 95,
                        jumpchoice = c("actual", "fit"),
                        verbose = FALSE, ...) {
  jumpchoice <- match.arg(jumpchoice)
  y.h <- max(object$y) + (1:h)
  x.h <- x.h %||% object$x
  
  # Predict ts model
  W <- predict(object$VAR, h, level)
  L <- W$conf.intervals
  L$mean <- W$predicted.values
  frM <- object$fitted.raw.moments
  sg <- sign(convertMoments(tail(frM, 1), from = "raw", to = "normalized"))
  
  fn1 <- function(z) {
    nM <- t(exp(z) * as.numeric(sg))
    rM <- convertMoments(nM, from = "normalized", to = "raw")
    return(rM)
  }
  fn2 <- function(z) {
    px  <- findDensity(z, x = x.h)$density
    out <- correct_jump_off(px, object, jumpchoice, h)
    return(out)
  }
  
  rM <- lapply(L, fn1)
  px <- lapply(rM, fn2)
  
  # Output preparation
  N <- length(px)
  CI  <- list(predicted.values = px[-N], predicted.raw.moments = rM[-N])
  out <- list(call = match.call(),
              predicted.values = px[[N]], 
              predicted.raw.moments = rM[[N]], 
              conf.intervals = CI, VAR = W, x = x.h, y = y.h)
  out <- structure(class = 'predict.MEM', out)
  return(out)
}


#' Jump-off correction
#' @param X matrix of fitted values
#' @inheritParams predict.MEM
#' @keywords internal
#' @export
correct_jump_off <- function(X, object, jumpchoice, h){
  if (jumpchoice == 'actual') {
    n   <- ncol(fitted(object))
    JO  <- with(object, observed.values[, n] / fitted.values[, n]) # jump-off vector
    foo <- function(x) x/sum(x)
    m   <- 40
    L   <- ifelse(h >= m, m, h)
    JM  <- X * 0 + 1  # jump-off matrix
    for (i in 1:length(JO)) JM[i, 1:L] = seq(JO[i], 1, length.out = m)[1:L] 
    
    X.out <- apply(X * JM, 2, FUN = foo)
  } else {
    X.out <- X
  }
  return(X.out)
}  


#' Print function for \code{\link{predict.MEM}}
#' @inherit print.MEM
#' @keywords internal
#' @export
print.predict.MEM <- function(x, ...) {
  cat('\nForecast: Maximum Entropy Mortality Model')
  cat('\nCall    : '); print(x$call)
  cat('\nAges in forecast   : ', paste0(range(x$x), collapse = ' - '))
  cat('\nYears in forecast  : ', paste0(range(x$y), collapse = ' - '))
  cat('\nMoments in forecast: ', paste0("M", c(0, ncol(x$predicted.raw.moments) - 1), 
                                        collapse = ' - '))
  cat('\n')
}


