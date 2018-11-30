# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Tue Nov 27 15:04:19 2018
# --------------------------------------------------- #

#' Perform In-Sample Testing of Mortality Forecasts Over One Time Period
#' @param data A data.frame or a matrix containing mortality data 
#' with ages \code{x} as row and time \code{y} as column. This object should
#' contain the data to be used in fitting (training) and validation process 
#' as well.
#' @param y.fit Years to be considered in fitting.
#' @param y.for Years to be forecast.
#' @inheritParams do.MortalityModels
#' @inheritParams do.MortalityForecasts
#' @inheritParams evalAccuracy
#' @seealso 
#' \code{\link{do.BBackTesting}}
#' \code{\link{evalAccuracy.BackTesting}}
#' @examples 
#' 
#' x = 0:98              # Ages
#' y1 = 1980:2000        # Training period
#' y2 = 2001:2016        # Validation period
#' y  = c(y1, y2)
#' h = max(y2) - max(y1) # Forecasting horizon
#' 
#' D <- HMD_male$dx$GBRTENW[paste(x), paste(y)] # DATA
#' 
#' # Select various mortality models
#' MM <- c("MRWD", "Oeppen", "MEM6")
#' # Fit & Forecast the models 
#' B <- do.BackTesting(data = D, x = x,
#'                     y.fit = y1, y.for = y2,
#'                     data.in = "dx",
#'                     models = MM)
#' 
#' # Compute accuracy measures
#' # The measures can be computed for different indicators. Even if it is not 
#' # impossible to get a different classification and ranking the 
#' # outcome should be in general the same.
#' evalAccuracy(B, data.out = "mx")
#' evalAccuracy(B, data.out = "qx")
#' evalAccuracy(B, data.out = "dx")
#' 
#' # Rank the model's performance.
#' A <- evalAccuracy(B, data.out = "ex")
#' A
#' R <- do.Ranking(A)
#' R
#' 
#' # Visualize the results
#' plot(B, data.out = "mx", facet = "x")
#' plot(B, data.out = "mx", facet = "y") 
#' @export
do.BackTesting <- function(data, 
                           x, 
                           y.fit, 
                           y.for, 
                           data.in = c("qx", "mx", "dx", "lx"),
                           models, 
                           level = 95, 
                           jumpchoice = c("actual", "fit"), 
                           verbose = FALSE,
                           ...) {
  # Prepare data
  data.in    <- match.arg(data.in)
  jumpchoice <- match.arg(jumpchoice)
  input      <- as.list(environment())
  h          <- max(y.for) - min(y.for) + 1
  
  D <- list(training.set = data[paste(x), paste(y.fit)], 
            validation.set = data[paste(x), paste(y.for)])
  
  # Fit
  M <- do.MortalityModels(data = D[[1]], 
                          x = x, 
                          y = y.fit, 
                          data.in = data.in, 
                          models = models, 
                          verbose = FALSE)
  
  # Forecast
  P <- do.MortalityForecasts(object = M, 
                             h = h, 
                             level = level, 
                             jumpchoice = jumpchoice, 
                             verbose = FALSE)
  
  # Exit
  out <- list(input = input, 
              call = match.call(), 
              Datasets = D, 
              Fitted = M, 
              Forecast = P)
  
  out <- structure(class = "BackTesting", out)
  return(out)
}

