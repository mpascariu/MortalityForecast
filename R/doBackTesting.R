

#' Perform in-sample Testing of Mortality Forecasts
#' @param data A data.frame or a matrix containing mortality data 
#' with ages \code{x} as row and time \code{y} as column. This object should
#' contain the data to be used in fitting (traning) and validation process as well.
#' @param y.fit Years to be considered in fitting.
#' @param y.for Years to be forecast.
#' @inheritParams doMortalityModels
#' @inheritParams doForecasts
#' @inheritParams getAccuracy
#' @examples 
#' x = 0:98              # Ages
#' y1 = 1980:2000        # Training period
#' y2 = 2001:2016        # Validation period
#' y  = c(y1, y2)
#' h = max(y2) - max(y1) # Forecasting horizon
#' 
#' D <- MortalityForecast.data$dx[paste(x), paste(y)] # DATA
#' 
#' MM <- c("MRW", "MRWD", "LC", "CoDa", "M6")
#' B.ex <- doBackTesting(data = D, x = x, 
#'                       y.fit = y1, y.for = y2, 
#'                       data.in = "dx", 
#'                       data.out = "ex", 
#'                       models = MM)
#' B.ex$accuracy
#' 
#' plot(B.ex, facet = "x")
#' plot(B.ex, facet = "y")  
#' @export
doBackTesting <- function(data, x, y.fit, y.for, 
                          data.in = c("qx", "mx", "dx", "lx"),
                          data.out = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"), 
                          models = c("MRWD", "LC"), 
                          measures = c("ME", "MAE", "MAPE", "sMAPE", "MRAE", "MASE"), 
                          xa = NULL, ya = NULL,
                          level = 95, 
                          jumpchoice = c("actual", "fit"), 
                          verbose = FALSE,
                          ...) {
  # Prepare data
  data.in    <- match.arg(data.in)
  data.out   <- match.arg(data.out)
  jumpchoice <- match.arg(jumpchoice)
  input      <- as.list(environment())
  h <- max(y.for) - min(y.for) + 1
  training.set   <- data[paste(x), paste(y.fit)]
  validation.set <- data[paste(x), paste(y.for)]
  
  # Fit - Forecast - Check Accuracy
  M <- doMortalityModels(training.set, x, y.fit, data.in, models)
  P <- doForecasts(M, h, level, jumpchoice)
  A <- getAccuracy(P, validation.set, data.in, data.out, measures, xa, ya, ...)
  
  # Datasets used
  vs <- convertFx(x, data = validation.set, from = data.in, to = data.out, lx0 = 1, ...)
  ts <- getObserved(M, data.out)
  fs <- getForecasts(P, data.out)
  D  <- list(training.set = ts, validation.set = vs, forecasts = fs)
  
  # Output
  out <- list(input = input, datasets = D, MortalityModels = M, 
              Forecasts = P, accuracy = A)
  out <- structure(class = "doBackTesting", out)
  return(out)
}

