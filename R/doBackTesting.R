

#' Perform In-Sample Testing of Mortality Forecasts Over One Time Period
#' @param data A data.frame or a matrix containing mortality data 
#' with ages \code{x} as row and time \code{y} as column. This object should
#' contain the data to be used in fitting (traning) and validation process as well.
#' @param y.fit Years to be considered in fitting.
#' @param y.for Years to be forecast.
#' @param xa Ages to be considered in model accuracy evaluation. It can be used 
#' to calculate the measures on a subset of the results. If \code{xa = NULL} 
#' the entire age-range in \code{x} is considered. Default: \code{xa = x[-length(x)]}.
#' In default mode the last row of results (last row in the life table) is 
#' excluded from evaluation. Because of different methods applied to force closing of 
#' the life table (input-specific), the errors can be "false positives" in the last row. 
#' @param ya Years to be considered in accuracy computation. If \code{ya = NULL} 
#' (default) the entire period in \code{y} is considered.
#' @inheritParams doMortalityModels
#' @inheritParams doForecasts
#' @inheritParams getAccuracy
#' @seealso 
#' \code{\link{doBBackTesting}}
#' \code{\link{doMortalityModels}}
#' \code{\link{doForecasts}}
#' \code{\link{getAccuracy}}
#' @examples 
#' x = 0:98              # Ages
#' y1 = 1980:2000        # Training period
#' y2 = 2001:2016        # Validation period
#' y  = c(y1, y2)
#' h = max(y2) - max(y1) # Forecasting horizon
#' 
#' D <- MortalityForecast.data$dx[paste(x), paste(y)] # DATA
#' 
#' MM <- c("MRWD", "LC", "CoDa", "MEM6")
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
                          xa = x[-length(x)], ya = NULL,
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
  M <- doMortalityModels(data = training.set, x = x, y = y.fit, 
                         data.in = data.in, models = models, verbose = FALSE)
  P <- doForecasts(object = M, h = h, level = level, jumpchoice = jumpchoice, 
                   verbose = FALSE)
  A <- getAccuracy(object = P, data = validation.set, 
                   data.in = data.in, data.out = data.out, 
                   measures = measures, xa = xa, ya = ya, ...)
  
  # Datasets used
  vs <- convertFx(x = x, data = validation.set, from = data.in, to = data.out, 
                  lx0 = 1, ...)
  ts <- getObserved(object = M, data.out = data.out)
  fs <- getForecasts(object = P, data.out = data.out)
  D  <- list(training.set = ts, validation.set = vs, forecasts = fs)
  
  # Output
  out <- list(input = input, datasets = D, MortalityModels = M, 
              Forecasts = P, accuracy = A)
  out <- structure(class = "doBackTesting", out)
  return(out)
}
