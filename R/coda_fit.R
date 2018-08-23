
#' Fit CoDa Mortality Model
#' 
#' Fit Compositional Data (CoDa) model for forecasting the life table 
#' distribution of deaths. CoDa is a Lee-Carter type model. A key difference 
#' between the \insertCite{lee1992;textual}{MortalityForecast}
#' method and the Compositional Data model (CoDa) is that the former fits and 
#' forecasts the death rates (mx) while the latter is based on the life table 
#' death distribution (dx). 
#' \insertCite{@See @oeppen2008 and @bergeron2017;textual}{MortalityForecast} 
#' for a detail description and mathematical formulation.
#' 
#' @inheritParams doMortalityModels 
#' @return The output is an object of class \code{"coda"} with the components:
#'  \item{input}{List with arguments provided in input. Saved for convenience.}
#'  \item{call}{An unevaluated function call, that is, an unevaluated 
#'  expression which consists of the named function applied to the given arguments.}
#'  \item{coefficients}{Estimated coefficients.}
#'  \item{fitted.values}{Fitted values of the estimated CoDa model.}
#'  \item{residuals}{Deviance residuals.} 
#'  \item{x}{Vector of ages used in the fitting.} 
#'  \item{y}{Vector of years used in the fitting.} 
#' @seealso \code{\link{predict.coda}}
#' @references \insertAllCited{}
#' @examples
#' # Fit CoDa model
#' D <- MortalityForecast.data$dx
#' M <- coda(D)
#' summary(M)
#' 
#' # Forecast life expectancy
#' P <- predict(M, h = 20)
#' @export
coda <- function(data, x = NULL, y = NULL, verbose = TRUE, ...){
  input <- c(as.list(environment()))
  coda.input.check(input)
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
  fit          <- t(proj.fit/rowSums(proj.fit))
  resid        <- data - fit
  dimnames(fit) = dimnames(resid) = dimnames(data) <- list(x, y)
  
  out <- list(input = input, call = match.call(), fitted.values = fit, 
              coefficients = cf, residuals = resid, x = x, y = y)
  out <- structure(class = 'coda', out)
  return(out)
}


#' Validate input values
#' @param X A list with input arguments provided in \code{\link{coda}} function
#' @keywords internal
coda.input.check <- function(X) {
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
#' Residuals of the CoDa Mortality Model
#' @inheritParams summary.coda
#' @export
residuals.coda <- function(object, ...){
  structure(class = "residMF", as.matrix(object$residuals))
}


#' Print coda
#' @param x An object of class \code{"coda"}
#' @param ... Further arguments passed to or from other methods.
#' @keywords internal
#' @export
print.coda <- function(x, ...) {
  cat('\nFit  : Compositional-Data Lee-Carter Mortality Model')
  cat('\nModel: clr d[x,t] = a[x] + b[x]k[t]')
  cat('\nCall : '); print(x$call)
  cat('\nAges  in fit:', paste(range(x$x), collapse = ' - '))
  cat('\nYears in fit:', paste(range(x$y), collapse = ' - '))
  cat('\n')
}


#' Summary coda
#' @param object An object of class \code{"coda"}.
#' @inheritParams print.coda
#' @keywords internal
#' @export
summary.coda <- function(object, ...) {
  axbx <- data.frame(ax = object$coefficients$ax, 
                     bx = object$coefficients$bx,
                     row.names = object$x)
  kt <- data.frame(kt = object$coefficients$kt)
  out = structure(class = 'summary.coda', 
                  list(A = axbx, K = kt, call = object$call,
                       y = object$y, x_ = object$x))
  return(out)
}


#' Print summary.coda
#' @param x An object of class \code{"summary.coda"}.
#' @inheritParams print.coda
#' @keywords internal
#' @export
print.summary.coda <- function(x, ...){
  cat('\nFit  : Compositional-Data Lee-Carter Mortality Model')
  cat('\nModel: clr d[x,t] = a[x] + b[x]k[t]')
  cat('\nCoefficients:\n')
  A <- head_tail(x$A, digits = 5, hlength = 6, tlength = 6)
  K <- head_tail(data.frame(. = '|', y = as.integer(x$y), kt = x$K),
                 digits = 5, hlength = 6, tlength = 6)
  print(data.frame(A, K))
  cat('\n')
}

