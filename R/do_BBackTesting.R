# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Tue Nov 27 16:16:04 2018
# --------------------------------------------------- #

#' Perform Out-of-Sample Testing of Mortality Forecasts Over Multiple Time Periods
#' @inheritParams do.MortalityModels
#' @inheritParams do.BackTesting
#' @inheritParams build.scenarios
#' @seealso 
#' \code{\link{do.BackTesting}}
#' \code{\link{evalAccuracy.BBackTesting}}
#' \code{\link{evalRobustness.BBackTesting}}
#' @author Marius D. Pascariu
#' @examples 
#' x  <- 0:95
#' y  <- 1970:2016
#' dx <- HMD_male$dx$GBRTENW[paste(x), paste(y)]
#' 
#' BB <- do.BBackTesting(data = dx, x, y,
#'                       data.in = "dx", 
#'                       models = c("MRWD", "HyndmanUllah"),
#'                       strategy = c(20, 20, 3))
#' 
#' A <- evalAccuracy(BB, data.out = "ex")
#' A
#' 
#' R <- do.Ranking(A)
#' R
#' @export
do.BBackTesting <- function(data, 
                            x, 
                            y,
                            data.in = c("qx", "mx", "dx", "lx"),
                            models,
                            strategy = c(f = 20, h = 20, s = 2),
                            level = 95,
                            jumpchoice = c("actual", "fit"),
                            verbose = TRUE, 
                            ...) {
  
  data.in    <- match.arg(data.in)
  jumpchoice <- match.arg(jumpchoice)
  input <- as.list(environment())
  call  <- match.call()
  S     <- build.scenarios(y, strategy)
  
  # Do Back-testing
  nc <- nrow(S) # no. of cases
  B  <- list()
  for (k in 1:nc) {
    yf <- S[[k, "fit"]]
    yh <- S[[k, "forecast"]]
    y_ <- c(yf, yh)
    if (verbose) cat(paste0("\nTest scenario ", k, "/", nc, ": "))
    Bk <- do.BackTesting(data = data[, paste(y_)], 
                         x = x, 
                         y.fit = yf, 
                         y.for = yh, 
                         data.in = data.in,
                         models = models, 
                         level = level, 
                         jumpchoice = jumpchoice, 
                         verbose = FALSE, 
                         ...)
    B[[k]] <- Bk
    remove(k, yf, yh, y_, Bk)
  }
  names(B) <- paste0("S", 1:nc)
  out <- list(input = input, 
              call = call, 
              scenarios = S, 
              results = B)
  out <- structure(class = "BBackTesting", out)
  return(out)
}


#' Build Scenarios for do.BBackTesting function
#' 
#' This is a utility function used in the \code{\link{do.BBackTesting}} 
#' function to define the rolling evaluation windows.
#' @param strategy Fitting-Forecasting strategy. Format: numerical vector.
#' The strategy \code{c(f, h, s)} consists in \code{f} number of years to use 
#' in fitting, \code{h} number of years to forecast and the step \code{s} to 
#' roll the evaluation window.
#' @inheritParams do.MortalityModels
#' @seealso \code{\link{do.BBackTesting}}
#' @author Marius D. Pascariu
#' @examples 
#' y <- 1900:2016
#' 
#' # Strategy
#' S <- c(f = 20, h = 20, s = 1)
#' # This means that:
#' # - 20 years of data are used to fit the models;
#' # - forecast 20 years and evaluate the results;
#' # - roll the evaluation window 1 year and repeat the process.
#' 
#' # 78 scenarios are created using this years and strategy.
#' W <- build.scenarios(y, S)
#' W
#' 
#' # The first scenario will be:
#' W[[1, "fit"]]
#' W[[1, "forecast"]]
#' 
#' # The last scenario will be:
#' W[[78, "fit"]]
#' W[[78, "forecast"]]
#' @export
build.scenarios <- function(y, strategy = c(f = 20, h = 20, s = 2)) {
  # Build scenarios -  method 1
  S <- strategy
  bop_fit <- rev(seq(from = max(y) - S[1] - S[2] + 1, 
                     to = min(y), 
                     by = -1 * S[3]))
  
  eop_fit <- bop_fit + S[1] - 1
  bop_fc  <- eop_fit + 1
  eop_fc  <- bop_fc + S[2] - 1
  
  out <- tibble(
    scenario = paste0("S", 1:length(bop_fit)),
    fit = mapply(":", bop_fit, eop_fit, SIMPLIFY = FALSE),
    forecast = mapply(":", bop_fc, eop_fc, SIMPLIFY = FALSE)
  )
  
  return(out)
}
