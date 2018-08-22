
#' Fit Maximum-Entropy Mortality Model 
#' 
#' @param data A data.frame or a matrix containing mortality data (dx) 
#' with ages \code{x} as row and time \code{y} as column.
#' @param x Numerical vector indicating the ages to be considered in the 
#' fitting procedure. If \code{NULL} the number fitted data points is equal 
#' with the observed number of data points.
#' This argument can be used for example to estimate a density between 0 and 130 
#' given the fact that the input dataset contains values for 0-100 only.
#' @param y Numerical vector indicating the years in \code{data}. 
#' Default: \code{NULL}.
#' @param n The maximum order of the moments to be used.
#' @param lag Integer for the lag order in VAR model. Default: 1.
#' @param verbose Show progress bar? Logical. Default: \code{FALSE}.
#' @param ... Further arguments passed to \code{\link[vars]{VAR}} method.
#' @return The output is an object of the class \code{"MEM"} with the components:
#' @return \item{input}{List with arguments provided in input. Saved for convenience.}
#' @return \item{call}{An unevaluated function call, that is, an unevaluated 
#' expression which consists of the named function applied to the given arguments;}
#' @return \item{coefficients}{Estimated coefficients;}
#' @return \item{fitted.values}{Fitted values of the estimated model;}
#' @return \item{observed.values}{Observed values used in fitting;}
#' @return \item{fitted.raw.moments}{Fitted raw moments;}
#' @return \item{observed.raw.moments}{Observed raw moments of the input data;}
#' @return \item{residuals}{Deviance residuals;} 
#' @return \item{VAR}{Object containing the components of the fitted time 
#' series model to extrapolate moments;}
#' @return \item{x}{Vector of ages used in the fitting;} 
#' @return \item{y}{Vector of years used in the fitting.} 
#' @inheritParams vars::VAR
#' @seealso \code{\link{predict.MEM}}, \code{\link{findMoments}}; 
#' \code{\link[vars]{VAR}}
#' @examples 
#' # Data ----------------------------------------
#' x  <- 0:110
#' y  <- 1965:2014
#' dx <- MortalityForecast.data$dx[paste(x), paste(y)]
#' 
#' # Example 1 -----------------------------------
#' # Fit model
#' M1 <- MEM(data = dx, x = x, y = y, n = 5)
#' M1
#' 
#' # Forecast model
#' P1 <- predict(M1, h = 16)
#' P1 
#' @export
MEM <- function(data, x = NULL, y = NULL, n = 5, 
                   exogen = NULL, lag = 1, verbose = FALSE, ...) {
  input <- c(as.list(environment()))
  if (verbose) { # Set progress bar
    pb <- startpb(0, 10)
    on.exit(closepb(pb)); setpb(pb, 1)
  }
  AY  <- find_ages_and_years(data, x, y)
  x   <- AY$x
  y   <- AY$y
  M   <- findMoments(data, x, y, n)
  orM <- M$raw.moments
  nM  <- M$normalized.moments
  nMT <- log(abs(nM))
  sg  <- sign(nM)
  if (verbose) setpb(pb, 3)
  
  V   <- MRW(t(nMT), x = NULL, y, include.drift = TRUE)
  fnM <- fnM <- t(exp(fitted(V)) * as.numeric(sg[1,]))     # fitted normalized moments
  frM <- convertMoments(fnM, In = "normalized", Out = "raw") # fitted raw moments
  if (verbose) setpb(pb, 5)
  
  fD  <- findDensity(frM[-1, ], x)$density                # fitted dx
  fD  <- cbind(NA, fD) #make sure fD has the same dimension as oD data. Not sure it's OK.
  oD  <- apply(data, 2, FUN = function(x) x/sum(x)) # observed dx - same scale as fitted dx
  res <- oD - fD                                    # residuals
  dimnames(oD) = dimnames(fD) = dimnames(res) = list(x = x, y = y)
  
  out <- list(input = input, call = match.call(),
              fitted.values = fD, observed.values = oD, coefficients = coef(V),
              residuals = res, fitted.raw.moments = frM, observed.raw.moments = orM, 
              VAR = V, x = x, y = y)
  out <- structure(class = 'MEM', out)
  if (verbose) setpb(pb, 10)
  return(out)
}


#' Find age and year vectors
#' 
#' @inheritParams MEM
#' @keywords internal
#' @export
#' 
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
#' Computed deviance residuals for a fitted Maximum-Entropy Mortality (MEM) model.
#' @param object An object of class \code{\link{MEM}}.
#' @param ... Further arguments passed to or from other methods.
#' @examples 
#' x  <- 0:110
#' y  <- 1965:2014
#' dx <- MortalityForecast.data$dx[paste(x), paste(y)]
#' M  <- MEM(dx, x, y, n = 5)
#' residuals(M)
#' @export
residuals.MEM <- function(object, ...) {
  X   <- object
  R   <- X$observed.values - X$fitted.values
  out <- structure(class = "residuals.MEM", R)
  return(out)
}

#' Print function for MEM method
#' @param x An object of class \code{\link{MEM}}.
#' @inheritParams residuals.MEM
#' @keywords internal
#' @export
print.MEM <- function(x, ...) {
  cat('\nFit : Lenart-Pascariu Mortality Model')
  cat('\nCall: '); print(x$call)
  cat('\nAges in fit   : ', paste0(range(x$x), collapse = ' - '))
  cat('\nYears in fit  : ', paste0(range(x$y), collapse = ' - '))
  cat('\nMoments in fit: ', paste0("M", c(0, x$input$n), collapse = ' - '))
  cat('\n')
}
