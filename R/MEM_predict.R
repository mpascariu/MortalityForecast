
#' Predict Maximum-Entropy Mortality Model
#' 
#' @param object An object of class \code{\link{MEM}}.
#' @param h h Number of units of time (years) to be forecast in the future.
#' @param x.h Numerical vector indicating the ages to be considered in 
#' reconstruction of the density over the forecast horizon. 
#' If \code{NULL}, the number of estimated data points is equal to the number 
#' of fitted data points. This argument can be used for example to estimate a 
#' density between 0 and 130 given the fact that the model was fitted on a 
#' dataset containing values for 0-100 only.
#' @param level Significance level of the confidence interval. Default: 95.
#' @param jumpchoice Method used for computation of jumpchoice. 
#'  Possibilities: \code{"actual"} (use actual rates from final year) 
#'  and \code{"fit"} (use fitted rates).
#' @param verbose Show progress bar? Logical, default \code{TRUE}.
#' @param ... Ignored.
#' @seealso \code{\link{MEM}}
#' @examples
#' y  <- 1965:2014
#' x  <- 0:110
#' dx <- MortalityForecast.data$dx[paste(x), paste(y)]
#' M  <- MEM(dx, x, y, n = 5)
#' 
#' # Forecast the distribution 16 year ahead for ages: 0-110
#' P1 <- predict(object = M, h = 16)
#' 
#' # Forecast for ages: 0-130
#' P2 <- predict(object = M, h = 16, x.h = 0:130)
#' 
#' @export
predict.MEM <- function(object, h, x.h = NULL, level = 95,
                           jumpchoice = c("actual", "fit"),
                           verbose = FALSE, ...) {
  y.h <- max(object$y) + (1:h)
  jumpchoice <- match.arg(jumpchoice)
  if (is.null(x.h)) x.h <- object$x
  
  # Predict ts model
  W <- predict(object$VAR, h, level)
  L <- W$conf.intervals
  L$mean <- W$predicted.values
  frM <- object$fitted.raw.moments
  sg <- sign(convertMoments(tail(frM, 1), In = "raw", Out = "normalized"))
  
  fn1 <- function(z) {
    nM <- t(exp(z) * as.numeric(sg))
    rM <- convertMoments(nM, In = "normalized", Out = "raw")
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
    # lag <- object$input$lag
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
#' 
#' @inherit print.MEM
#' @keywords internal
#' @export
print.predict.MEM <- function(x, ...) {
  cat('\nForecast: Lenart-Pascariu Mortality Model')
  cat('\nCall    : '); print(x$call)
  cat('\nAges in forecast   : ', paste0(range(x$x), collapse = ' - '))
  cat('\nYears in forecast  : ', paste0(range(x$y), collapse = ' - '))
  cat('\nMoments in forecast: ', paste0("M", c(0, ncol(x$predicted.raw.moments) - 1), 
                                        collapse = ' - '))
  cat('\n')
}
