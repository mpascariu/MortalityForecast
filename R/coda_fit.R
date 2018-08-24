
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
#' @return The output is an object of class \code{fitOeppen} with the components:
#'  \item{input}{List with arguments provided in input. Saved for convenience.}
#'  \item{call}{An unevaluated function call, that is, an unevaluated 
#'  expression which consists of the named function applied to the given arguments.}
#'  \item{coefficients}{Estimated coefficients.}
#'  \item{fitted.values}{Fitted values of the estimated CoDa model.}
#'  \item{residuals}{Deviance residuals.} 
#'  \item{x}{Vector of ages used in the fitting.} 
#'  \item{y}{Vector of years used in the fitting.} 
#' @seealso 
#' \code{\link{predict.fitOeppen}}
#' \code{\link{plot.fitOeppen}}
#' @references \insertAllCited{}
#' @examples
#' # Data
#' x  <- 0:100
#' y  <- 1980:2014
#' dx <- MortalityForecast.data$dx[paste(x), paste(y)]
#' 
#' # Fit CoDa-LC model
#' M <- fitOeppen(dx, x, y)
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
fitOeppen <- function(data, x = NULL, y = NULL, verbose = TRUE, ...){
  input <- c(as.list(environment()))
  fitOeppen.input.check(input)
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
  out <- structure(class = 'fitOeppen', out)
  return(out)
}


#' Validate input values
#' @param X A list with input arguments provided in \code{\link{fitOeppen}} function
#' @keywords internal
fitOeppen.input.check <- function(X) {
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
#' @inheritParams summary.fitOeppen
#' @export
residuals.fitOeppen <- function(object, ...){
  structure(class = "residMF", as.matrix(object$residuals))
}


#' Print fitOeppen
#' @param x An object of class \code{"fitOeppen"}
#' @param ... Further arguments passed to or from other methods.
#' @keywords internal
#' @export
print.fitOeppen <- function(x, ...) {
  cat('\nFit  :', x$info$name)
  cat('\nModel:', x$info$formula)
  cat('\nCall : '); print(x$call)
  cat('\nAges  in fit:', paste(range(x$x), collapse = ' - '))
  cat('\nYears in fit:', paste(range(x$y), collapse = ' - '))
  cat('\n')
}


#' Summary fitOeppen
#' @param object An object of class \code{"fitOeppen"}.
#' @inheritParams print.fitOeppen
#' @keywords internal
#' @export
summary.fitOeppen <- function(object, ...) {
  axbx <- data.frame(ax = object$coefficients$ax, 
                     bx = object$coefficients$bx,
                     row.names = object$x)
  kt <- data.frame(kt = object$coefficients$kt)
  out = structure(class = 'summary.fitOeppen', 
                  list(A = axbx, K = kt, call = object$call,
                       y = object$y, x_ = object$x))
  return(out)
}


#' Print summary.fitOeppen
#' @param x An object of class \code{"summary.fitOeppen"}.
#' @inheritParams print.fitOeppen
#' @keywords internal
#' @export
print.summary.fitOeppen <- function(x, ...){
  cat('\nFit  :', x$info$name)
  cat('\nModel:', x$info$formula)
  cat('\nCoefficients:\n')
  A <- head_tail(x$A, digits = 5, hlength = 6, tlength = 6)
  K <- head_tail(data.frame(. = '|', y = as.integer(x$y), kt = x$K),
                 digits = 5, hlength = 6, tlength = 6)
  print(data.frame(A, K))
  cat('\n')
}

