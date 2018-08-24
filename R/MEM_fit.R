
#' The Maximum-Entropy Mortality Model 
#' 
#' @param n The maximum order of the moments to be used.
#' @inheritParams doMortalityModels
#' @return The output is an object of class \code{"fitMaxEntMortality"} with 
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
#' \code{\link{predict.fitMaxEntMortality}} 
#' \code{\link{plot.fitMaxEntMortality}} 
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
#' M <- fitMaxEntMortality(data = dx, x = x, y = y, n = 5)
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
fitMaxEntMortality <- function(data, x = NULL, y = NULL, n = 5, verbose = FALSE, ...) {
  input <- c(as.list(environment()))

  AY  <- find_ages_and_years(data, x, y)
  x   <- AY$x
  y   <- AY$y
  M   <- findMoments(data, x, y, n)
  orM <- M$raw.moments
  nM  <- M$normalized.moments
  nMT <- log(abs(nM))
  sg  <- sign(nM)
  
  V   <- fitMRandomWalk(t(nMT), x = NULL, y, include.drift = TRUE)
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
  out <- structure(class = 'fitMaxEntMortality', out)
  return(out)
}


#' Find age and year vectors
#' @inheritParams fitMaxEntMortality
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
#' @param object An object of class \code{\link{fitMaxEntMortality}}.
#' @param ... Further arguments passed to or from other methods.
#' @examples 
#' x  <- 0:110
#' y  <- 1965:2014
#' dx <- MortalityForecast.data$dx[paste(x), paste(y)]
#' M  <- fitMaxEntMortality(dx, x, y, n = 5)
#' residuals(M)
#' @export
residuals.fitMaxEntMortality <- function(object, ...) {
  structure(class = "residMF", as.matrix(object$residuals))
}


#' Print function for Maximum-Entropy Mortality Model 
#' @param x An object of class \code{\link{fitMaxEntMortality}}.
#' @inheritParams residuals.fitMaxEntMortality
#' @keywords internal
#' @export
print.fitMaxEntMortality <- function(x, ...) {
  cat('\nFit  :', x$info$name)
  cat('\nModel:', x$info$formula)
  cat('\nCall : '); print(x$call)
  cat('\nAges in fit   : ', paste0(range(x$x), collapse = ' - '))
  cat('\nYears in fit  : ', paste0(range(x$y), collapse = ' - '))
  cat('\nMoments in fit: ', paste0("M", c(0, x$input$n), collapse = ' - '))
  cat('\n')
}
