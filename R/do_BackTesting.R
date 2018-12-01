# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Tue Nov 27 15:04:19 2018
# --------------------------------------------------- #

#' Perform Out-of-Sample Testing of Mortality Forecasts Over One Time Period
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
#' library(MortalityForecast)
#' library(magrittr)
#' 
#' x  <- 0:100            # Ages
#' y1 <- 1980:2000        # Training period
#' y2 <- 2001:2016        # Validation period
#' y  <- c(y1, y2)        # entire period
#' K <- "GBRTENW"         # Country to be analysed
#' 
#' # Mortality data for country: K
#' mx.country <- HMD_male$mx[[K]][paste(x), paste(y)] %>% replace.zeros
#' 
#' # Create a mortality data for a benchmark population
#' # to be used in "LiLee" and "OeppenC" models
#' average_mx <- function(w) prod(w, na.rm = TRUE)^(1/(length(w)))
#' 
#' Mx <- HMD_male$mx %>% 
#'   lapply(as.matrix) %>% 
#'   lapply(replace.zeros) %>% 
#'   simplify2array() %>% 
#'   apply(., 1:2, FUN = average_mx) %>% 
#'   as.data.frame()
#' 
#' mx.benchmark <- Mx[paste(x), paste(y)]
#' 
#' # Select various mortality models
#' MM <- c("MRWD", "LeeCarter", "LiLee", "HyndmanUllah", 
#'         "Oeppen", "OeppenC", "MEM5", "MEM6")
#' # Fit & Forecast the models 
#' B <- do.BackTesting(data = mx.country,
#'                     data.B = mx.benchmark,
#'                     x = x,
#'                     y.fit = y1, 
#'                     y.for = y2,
#'                     data.in = "mx",
#'                     models = MM)
#' 
#' # Compute accuracy measures for resulted life exectancies
#' A <- evalAccuracy(B, data.out = "ex")
#' A
#' # Rank the model's performance.
#' R <- do.Ranking(A)
#' R
#' R[, c(1:3, 20)]
#' 
#' # Visualize the results
#' plot(B, data.out = "ex", facet = "x", 
#'      which = c(0, 20, 40, 60, 75, 90))
#' plot(B, data.out = "mx", facet = "y", which = 2016) 
#' @export
do.BackTesting <- function(data, 
                           data.B = NULL,
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
  B <- list(training.set = data.B[paste(x), paste(y.fit)], 
            validation.set = data.B[paste(x), paste(y.for)])
  
  # Fit
  M <- do.MortalityModels(data = D[[1]],
                          data.B = B[[1]],
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

