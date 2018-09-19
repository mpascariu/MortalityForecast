
#' Multivariate Random Walk Model
#' 
#' Fit a Multivariate Random Walk Model to the input \code{data}, a
#' multivariate time series.
#' @details 
#' For further information on the Multivariate Random Walk with
#' drift see Appendix B in \insertCite{haberman2011;textual}{MortalityForecast}.
#' @param data Numeric matrix with a multivariate time series.
#' Series are arranged in rows with columns representing time.
#' @param x Numerical vector indicating the age classes in input \code{data}. 
#' Optional. This is used to label the output tables and related plots.
#' Default: \code{NULL}.
#' @param y Numerical vector indicating the years in input \code{data}.
#' Optional. This is used to label the output tables and related plots.
#' Default: \code{NULL}.
#' @param include.drift Should the Random Walk model include a linear drift term? 
#' Default: \code{TRUE}.
#' @inheritParams doMortalityModels
#' @return An object of class \code{MRW} with components:
#'  \item{input}{A list with the input data.}
#'  \item{drift}{A vector with the estimated drift.}
#'  \item{sigma}{A matrix with the estimated variance covariance matrix.}
#'  \item{fitted}{Fitted values.}
#'  \item{residuals}{Residuals from the fitted model. That is 
#'  observed minus fitted values.}
#' @references \insertAllCited{}
#' @source The original implementation of this function was taken from 
#' \code{StMoMo} R package.
#' @author Marius D. Pascariu
#' @examples 
#' # Forecast mortality using a Multivariate Random Walk with Drift model
#' x  <- 0:100
#' y  <- 1980:2000
#' mx <- MortalityForecast.data$mx[paste(x), paste(y)]
#' 
#' M <- fit_MRW(data = log(mx), x = x, y = y, include.drift = TRUE)
#' P <- predict(M, h = 16)
#' R <- residuals(M)
#' 
#' pv <- exp(P$predicted.values)
#' matplot(pv, type = "l", log = "y")
#' matplot(t(pv), type = "l", log = "y")
#' @export
fit_MRW <- function(data, x = NULL, y = NULL, include.drift = TRUE, ...) {
  # Save the input
  input <- as.list(environment())
  call  <- match.call()
  
  # Prepare data
  x <- x %||% 1:nrow(data)
  y <- y %||% 1:ncol(data)
  data  <- as.matrix(data)
  dimnames(data) <- list(x = x, y = y)
  if (ncol(data) == 1L) data <- t(data)
  ny <- ncol(data)
  nx <- nrow(data)
  
  # Drift
  d <- t(colMeans(diff(t(data)), na.rm = TRUE))
  if (!include.drift) d = d * 0
  
  # Fit the model
  fit  <- cbind(array(NA, c(nx, 1)), data[, -ny] + array(d, c(nx, ny - 1)))
  res  <- data - fit
  sigma <- cov(t(res), use = "complete.obs")
  dimnames(fit) <- dimnames(res) <- dimnames(data)
  
  # Exit
  out <- list(input = input, call = call, drift = d, sigma = sigma, 
              fitted = fit, residuals = res, x = x, y = y)
  out <- structure(class = "MRW", out)  
  return(out)
}


# S3 ----------------------------------------------


#' Forecast a Multivariate Random Walk Model
#' 
#' Returns forecasts and other information for a Multivariate 
#' Random Walk Model.
#' 
#' @param object An object of class \code{MRW}.
#' @inheritParams doForecasts
#' @return An object of the class \code{predict.MRW} with components:
#' \item{predicted.values}{data.frame with the central forecast.}
#' \item{conf.intervals}{List with lower and upper limits for prediction intervals.}
#' \item{x}{Numerical vector indicating the age classes in output.}
#' \item{y}{Numerical vector indicating the years in output.}
#' @source The original implementation of this function was taken from 
#' \code{StMoMo} R package.
#' @author Marius D. Pascariu
#' @export
predict.MRW <- function(object, h = 10, level = c(80, 95), ...) {
  
  data <- object$input$data
  x  <- object$x
  y  <- object$y
  nx <- length(x)  
  ny <- length(y)
  nn <- 1:h
  xf <- x
  yf <- max(y) + nn
  
  mean <- data[, ny] + t(array(nn, c(h, nx))) * array(object$drift, c(nx, h))
  dimnames(mean) <- list(x = xf, y = yf)
  
  sigma <- object$sigma
  se    <- sqrt(t(array(nn, c(h, nx))) * array(diag(sigma), c(nx, h)))
  nconf <- length(level)
  z     <- qnorm((1 - level/100)/2)
  lower <- upper <- list()
  
  for (i in 1:nconf) {
    lower[[i]] <- mean - z[i] * se
    upper[[i]] <- mean + z[i] * se
  }
  
  dn <- apply(expand.grid(c("L", "U"), level), 1, paste, collapse = "")
  CI <- c(lower, upper)
  names(CI) <- dn
  out <- list(predicted.values = mean, conf.intervals = CI,
              x = xf, y = yf)
  out <- structure(class = "predict.MRW", out)  
  return(out)
}


#' Residuals of a Multivariate Random-Walk Model 
#' 
#' Computed deviance residuals for a Multivariate Random-Walk model .
#' @param object An object of class \code{MRW}.
#' @param ... Further arguments passed to or from other methods.
#' @seealso \code{\link{fit_MRW}}
#' @export
residuals.MRW <- function(object, ...) {
  structure(class = "residMF", as.matrix(object$residuals))
}


#' Print function for MRW method
#' @param x An object of class \code{MRW}.
#' @inheritParams residuals.MRW
#' @keywords internal
#' @export
print.MRW <- function(x, ...) {
  cat('\nFit  : Multivariate Random-Walk Model')
  cat('\nDrift:', x$input$include.drift)
  cat('\nCall : '); print(x$call)
  cat('\nAges in fit : ', paste0(range(x$x), collapse = ' - '))
  cat('\nYears in fit: ', paste0(range(x$y), collapse = ' - '))
  cat('\n')
}



