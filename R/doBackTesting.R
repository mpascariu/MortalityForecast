

#' Perform in-sample testing of mortality forecasts
#' @param data Mortality data for fitting of the models and for validation; 
#' that is the data coresponding to the horizon to be predicted.
#' @param y.fit Years to be considered in fitting.
#' @param y.for Years to be forecast.
#' @inheritParams doMortalityModels
#' @inheritParams doForecasts
#' @inheritParams getAccuracy
#' @examples 
#' x = 0:98
#' y1 = 1980:1999
#' y2 = 2000:2016
#' y  = c(y1, y2)
#' h = max(y2) - max(y1)
#' 
#' D <- dxForecast::dxForecast.data$dx$male[paste(x), paste(y)]
#' ex <- dxForecast::dxForecast.data$ex$male
#' exogen <- ex[paste(y1)]
#' 
#' MM <- c("LC", "CoDa", "M6")
#' B.ex <- doBackTesting(data = D, x = x, 
#'                       y.fit = y1, y.for = y2, 
#'                       data.type = "dx", 
#'                       what = "mx", 
#'                       exogen = exogen,
#'                       models = MM)
#' 
#' plot(B.ex, facet = "x")
#' 
#' plot(B.ex, facet = "y") #+ 
#'   #scale_y_continuous(trans='log10')
#' 
#' @export
doBackTesting <- function(data, x, 
                          y.fit, y.for, data.type, what, 
                          exogen = NULL, ci = 95, jumpchoice = "actual", 
                          xa = NULL, ya = NULL,
                          measures = c("ME", "MAE", "MAPE", "sMAPE", "MRAE", "MASE"), 
                          models = c("LC", "FDM", "PLAT", "CoDa", "M4", "M4X", 
                                     "M5", "M5X", "M6", "M6X"),
                          ...) {
  # Prepare data
  input <- as.list(environment())
  h <- max(y.for) - min(y.for) + 1
  training.set   <- data[paste(x), paste(y.fit)]
  validation.set <- data[paste(x), paste(y.for)]
  if (!is.null(exogen)) exogen <- exogen[paste(y.fit)]
  
  # Fit - Forecast - Check Accuracy
  M <- doMortalityModels(training.set, x, y.fit, data.type, models, exogen = exogen)
  P <- doForecasts(M, h, ci, jumpchoice)
  A <- getAccuracy(P, validation.set, xa, ya, data.type, what, measures, ...)
  
  # Datasets used
  vs <- convertFx(x, data = validation.set, In = data.type, Out = what, lx0 = 1, ...)
  ts <- getObserved(M, what)
  fs <- getForecasts(P, what)
  D  <- list(training.set = ts, validation.set = vs, forecasts = fs)
  
  # Output
  out <- list(input = input, datasets = D, MortalityModels = M, Forecasts = P, 
              accuracy = A)
  out <- structure(class = "doBackTesting", out)
  return(out)
}

