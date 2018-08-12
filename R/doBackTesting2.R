


# y.fit = 1980:1999
# y.for = 2000:2016
# y = c(y.fit, y.for)
# x = 0:100
# D <- dxForecast::dxForecast.data$dx$male[paste(x), paste(y)]
# ex <- dxForecast::dxForecast.data$ex$male
# exogen <- ex[paste(y)]

dbt <- function(data, x, y.fit, y.for, data.type, what, exogen, 
                ci = 95, jumpchoice = "actual", 
                measures = c("ME", "MAE", "MAPE", "sMAPE", "MRAE", "MASE")) {
  input <- as.list(environment())
  y <- c(y.fit, y.for)
  h <- max(y.for) - min(y.for) + 1
  
  
  training.set   <- data[paste(x), paste(y.fit)]
  validation.set <- data[paste(x), paste(y.for)]
  exogen <- exogen[paste(y.fit)]
  
  M <- doMortalityModels(training.set, x, y.fit, data.type, exogen = exogen)
  P <- doForecasts(M, h, ci, jumpchoice)
  A <- getAccuracy(P, validation.set, x, y, data.type, what)
  
  vs <- convertFx(x, validation.set, In = data.type, Out = what, lx0 = 1, ...)
  ts <- getObserved(M, what)
  fs <- getForecasts(P, what)
  out <- list(accuracy = A)
  return(out)
}

# dbt(D, x, y.fit, y.for, data.type = "dx", exogen)
