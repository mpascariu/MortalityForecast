
#' Forecast a Multivariate Random Walk Model
#' 
#' Returns forecasts and other information for a Multivariate 
#' Random Walk Model.
#' 
#' @param object An object of class \code{"fitMRandomWalk"}.
#' @inheritParams doForecasts
#' @return An object of the class \code{"predict.fitMRandomWalk"} with components:
#' \item{predicted.values}{ Data.frame with the central forecast.}
#' \item{conf.intervals}{ List with lower and upper limits for prediction intervals.}
#' \item{x}{ Numerical vector indicating the age classes in output.}
#' \item{y}{ Numerical vector indicating the years in output.}
#' @export
predict.fitMRandomWalk <- function(object, h = 10, level = c(80, 95), ...) {
  
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
  out <- structure(class = "predict.fitMRandomWalk", out)  
  return(out)
}

