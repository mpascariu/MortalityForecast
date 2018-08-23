# Thu Aug 23 14:08:14 2018 --------- Marius D. Pascariu ---


#' Perform In-Sample Testing of Mortality Forecasts Over Multiple Time Periods
#' @inheritParams doMortalityModels
#' @inheritParams doBackTesting
#' @inheritParams buildScenarios
#' @seealso 
#' \code{\link{doBackTesting}}
#' \code{\link{doMortalityModels}}
#' \code{\link{doForecasts}}
#' \code{\link{getAccuracy}}
#' @examples 
#' x  <- 0:95
#' y  <- 1970:2016
#' dx <- MortalityForecast.data$dx[paste(x), paste(y)]
#' 
#' BB <- doBBackTesting(data = dx, x, y,
#'                      data.in = "dx", 
#'                      data.out = "ex", 
#'                      models = c("MRWD", "LC"),
#'                      strategy = c(20, 20, 5))
#' @export
doBBackTesting <- function(data, x, y,
                           data.in = c("qx", "mx", "dx", "lx"),
                           data.out = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"),
                           models = c("MRWD", "LC"),
                           measures = c("ME", "MAE", "MAPE", "sMAPE", "MRAE", "MASE"),
                           strategy = c(f = 20, h = 20, s = 2),
                           xa = x[-length(x)], ya = NULL,
                           level = 95,
                           jumpchoice = c("actual", "fit"),
                           verbose = TRUE, ...) {
  
  data.in <- match.arg(data.in)
  data.out <- match.arg(data.out)
  jumpchoice <- match.arg(jumpchoice)
  call <- match.call()
  S <- buildScenarios(y, strategy)
  
  # Do Back-testing
  nc = nrow(S) # no. of cases
  B <- list()
  A <- 0
  for (k in 1:nc) {
    yf <- S[[k, "fit"]]
    yh <- S[[k, "forecast"]]
    y_ <- c(yf, yh)
    if (verbose) cat(paste0("\nTest scenario ", k, "/", nc, ": "))
    Bk <- doBackTesting(data = data[, paste(y_)], 
                        x = x, 
                        y.fit = yf, 
                        y.for = yh, 
                        data.in = "dx", data.out = "mx",
                        models = models, xa = xa, ya = ya, ...)
    Ak <- Bk$accuracy$results
    
    B[[k]] <- Bk
    A <- A + Ak
    remove(k, yf, yh, y_, Bk, Ak)
  }
  
  out <- list(call = call, accuracy = A/nc, scenarios = S, results = B)
  out <- structure(class = "doBackTESTING", out)
  return(out)
}


#' Build Scenarios for doBackTESTING
#' @param strategy Fitting-Forecasting strategy. Format: numerical vector.
#' The strategy \code{c(f, h, s)} consists in \code{f} number of years to use in 
#' fitting, \code{h} number of years to forecast and step \code{s} to change the
#' fitting-forecasting window.
#' @inheritParams doMortalityModels
#' @keywords internal
buildScenarios <- function(y, strategy = c(f = 20, h = 20, s = 2)) {
  # Build scenarios -  method 1
  S <- strategy
  bop_fit = seq(from = max(y) - S[1] - S[2] + 1, to = min(y), by = -1 * S[3])
  eop_fit = bop_fit + S[1] - 1
  bop_fc  = eop_fit + 1
  eop_fc  = bop_fc + S[2] - 1
  
  out <- tibble(
    scenario = paste0("S", 1:length(bop_fit)),
    fit = mapply(":", bop_fit, eop_fit, SIMPLIFY = F),
    forecast = mapply(":", bop_fc, eop_fc, SIMPLIFY = F)
  )
  
  return(out)
}
