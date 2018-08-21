
#' Fit a Multivariate Random Walk Model
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
#' @param ... Ignored.
#' @return An object of class the \code{"MRW"} with components:
#' \item{input}{ A list with the input data.}
#' \item{drift}{ A vector with the estimated drift.}
#' \item{vcov}{ A matrix with the estimated variance covariance matrix.}
#' \item{fitted}{ Fitted values.}
#' \item{residuals}{ Residuals from the fitted model. That is 
#' observed minus fitted values.}
#' @references \insertAllCited{}
#' @examples 
#' # Forecast mortality using a Random Walk with Drift model
#' x <- 0:100
#' y <- 1980:2000
#' D <- MortalityForecast.data$mx[paste(x), paste(y)]
#' logD <- log(D)
#' 
#' M <- MRW(data = log(D), x = x, y = y, include.drift = TRUE)
#' P <- predict(M, h = 16)
#' 
#' pv <- exp(P$predicted.values)
#' matplot(pv, type = "l", log = "y")
#' @export
MRW <- function(data, x = NULL, y = NULL, include.drift = TRUE, ...) {
  # Save the input
  input <- as.list(environment())
  
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
  vcov <- cov(t(res), use = "complete.obs")
  dimnames(fit) <- dimnames(res) <- dimnames(data)
  
  # Exit
  out <- list(input = input, drift = d, vcov = vcov, fitted = fit, 
              residuals = res, x = x, y = y)
  out <- structure(class = "MRW", out)  
  return(out)
}


#' Forecast a Multivariate Random Walk Model
#' 
#' Returns forecasts and other information for a Multivariate 
#' Random Walk Model.
#' 
#' @param object An object of the class \code{"MWR"}.
#' @param h Number of periods for forecasting.
#' @param level Confidence level for prediction intervals.
#' @param ... Ignored.
#' @return An object of the class \code{"predict.MRW"} with components:
#' \item{mean}{ Array with the central forecast.}
#' \item{lower}{ Three dimensional array with lower limits for prediction intervals.}
#' \item{upper}{ Three dimensional array with upper limits for prediction intervals.}
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
  
  vcov  <- object$vcov
  se    <- sqrt(t(array(nn, c(h, nx))) * array(diag(vcov), c(nx, h)))
  nconf <- length(level)
  z     <- qnorm((1 - level/100)/2)
  dn    <- list(rownames(data), yf, paste0(level, "%"))
  lower <- upper <- array(NA, c(nx, h, nconf), dimnames = dn)
  
  for (i in 1:nconf) {
    lower[, , i] <- mean - z[i] * se
    upper[, , i] <- mean + z[i] * se
  }
  out <- list(predicted.values = mean, lower = lower, upper = upper,
              x = xf, y = yf)
  out <- structure(class = "predict.MRW", out)  
  return(out)
}

