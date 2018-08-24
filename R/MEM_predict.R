
#' Predict Maximum-Entropy Mortality Model
#' 
#' @param object An object of class \code{\link{fitMaxEntMortality}}.
#' @param x.h Numerical vector indicating the ages to be considered in 
#' reconstruction of the density over the forecast horizon. 
#' If \code{NULL}, the number of estimated data points is equal to the number 
#' of fitted data points. This argument can be used for example to estimate a 
#' density between 0 and 130 given the fact that the model was fitted on a 
#' dataset containing values for 0-100 only.
#' @inheritParams doForecasts
#' @seealso \code{\link{fitMaxEntMortality}}
#' @examples
#' x  <- 0:110
#' y  <- 1985:2016
#' dx <- MortalityForecast.data$dx[paste(x), paste(y)]
#' M  <- fitMaxEntMortality(dx, x, y, n = 5)
#' P  <- predict(M, h = 16, x.h = 0:120)
#' 
#' plot(P, plotType = "mean")
#' plot(P, plotType = "lower")
#' plot(P, plotType = "upper")
#' 
#' plot(P, M, plotType = "raw_moments")
#' plot(P, M, plotType = "normalised_moments")
#' @export
predict.fitMaxEntMortality <- function(object, h, x.h = NULL, level = 95,
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
  out <- structure(class = 'predict.fitMaxEntMortality', out)
  return(out)
}


#' Jump-off correction
#' @param X matrix of fitted values
#' @inheritParams predict.fitMaxEntMortality
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


#' Print function for \code{\link{predict.fitMaxEntMortality}}
#' @inherit print.fitMaxEntMortality
#' @keywords internal
#' @export
print.predict.fitMaxEntMortality <- function(x, ...) {
  cat('\nForecast: Maximum Entropy Mortality Model')
  cat('\nCall    : '); print(x$call)
  cat('\nAges in forecast   : ', paste0(range(x$x), collapse = ' - '))
  cat('\nYears in forecast  : ', paste0(range(x$y), collapse = ' - '))
  cat('\nMoments in forecast: ', paste0("M", c(0, ncol(x$predicted.raw.moments) - 1), 
                                        collapse = ' - '))
  cat('\n')
}
