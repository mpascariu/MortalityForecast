

#' Perform in-sample testing of mortality forecasts
#' @param data Mortality data for fitting of the models and for validation; 
#' that is the data coresponding to the horizon to be predicted.
#' @param y.fit Years to be considered in fitting
#' @param y.for Years to be forecast
#' @inheritParams doMortalityModels
#' @inheritParams doForecasts
#' @inheritParams getAccuracy
#' @export
doBackTesting <- function(data, x, y.fit, y.for, data.type, what, exogen, 
                ci = 95, jumpchoice = "actual", 
                measures = c("ME", "MAE", "MAPE", "sMAPE", "MRAE", "MASE"), ...) {
  # Prepare data
  input <- as.list(environment())
  y <- c(y.fit, y.for)
  h <- max(y.for) - min(y.for) + 1
  training.set   <- data[paste(x), paste(y.fit)]
  validation.set <- data[paste(x), paste(y.for)]
  exogen <- exogen[paste(y.fit)]
  
  # Fit - Forecast - Check Accuracy
  M <- doMortalityModels(training.set, x, y.fit, data.type, exogen = exogen)
  P <- doForecasts(M, h, ci, jumpchoice)
  A <- getAccuracy(P, validation.set, x, y, data.type, what, measures, ...)
  
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

