
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
#' @inheritParams doMortalityModels
#' @return An object of class \code{"MRW"} with components:
#'  \item{input}{ A list with the input data.}
#'  \item{drift}{ A vector with the estimated drift.}
#'  \item{sigma}{ A matrix with the estimated variance covariance matrix.}
#'  \item{fitted}{ Fitted values.}
#'  \item{residuals}{ Residuals from the fitted model. That is 
#'  observed minus fitted values.}
#' @references \insertAllCited{}
#' @examples 
#' # Forecast mortality using a Multivariate Random Walk with Drift model
#' x  <- 0:100
#' y  <- 1980:2000
#' mx <- MortalityForecast.data$mx[paste(x), paste(y)]
#' 
#' M <- MRW(data = log(mx), x = x, y = y, include.drift = TRUE)
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
  sigma <- cov(t(res), use = "complete.obs")
  dimnames(fit) <- dimnames(res) <- dimnames(data)
  
  # Exit
  out <- list(input = input, drift = d, sigma = sigma, fitted = fit, 
              residuals = res, x = x, y = y)
  out <- structure(class = "MRW", out)  
  return(out)
}

