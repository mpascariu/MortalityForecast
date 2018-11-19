# Thu Aug 23 14:08:14 2018 --------- Marius D. Pascariu ---


#' Perform In-Sample Testing of Mortality Forecasts Over Multiple Time Periods
#' @inheritParams doMortalityModels
#' @inheritParams doBackTesting
#' @inheritParams buildScenarios
#' @seealso 
#' \code{\link{doBackTesting}}
#' \code{\link{evalAccuracy.doBBackTesting}}
#' \code{\link{evalRobustness.doBBackTesting}}
#' @author Marius D. Pascariu
#' @examples 
#' x  <- 0:95
#' y  <- 1970:2016
#' dx <- HMD_male$dx$GBRTENW[paste(x), paste(y)]
#' 
#' BB <- doBBackTesting(data = dx, x, y,
#'                      data.in = "dx", 
#'                      models = c("MRWD", "HyndmanUllah"),
#'                      strategy = c(20, 20, 3))
#' 
#' A <- evalAccuracy(BB, data.out = "ex")
#' A
#' 
#' R <- doRanking(A)
#' R
#' @export
doBBackTesting <- function(data, x, y,
                           data.in = c("qx", "mx", "dx", "lx"),
                           models,
                           strategy = c(f = 20, h = 20, s = 2),
                           level = 95,
                           jumpchoice = c("actual", "fit"),
                           verbose = TRUE, ...) {
  
  data.in    <- match.arg(data.in)
  jumpchoice <- match.arg(jumpchoice)
  input <- as.list(environment())
  call  <- match.call()
  S <- buildScenarios(y, strategy)
  
  # Do Back-testing
  nc <- nrow(S) # no. of cases
  B  <- list()
  for (k in 1:nc) {
    yf <- S[[k, "fit"]]
    yh <- S[[k, "forecast"]]
    y_ <- c(yf, yh)
    if (verbose) cat(paste0("\nTest scenario ", k, "/", nc, ": "))
    Bk <- doBackTesting(data = data[, paste(y_)], 
                        x = x, y.fit = yf, y.for = yh, data.in = data.in,
                        models = models, level = level, jumpchoice = jumpchoice, 
                        verbose = FALSE, ...)
    B[[k]] <- Bk
    remove(k, yf, yh, y_, Bk)
  }
  names(B) <- paste0("S", 1:nc)
  out <- list(input = input, call = call, scenarios = S, results = B)
  out <- structure(class = "doBBackTesting", out)
  return(out)
}


#' Build Scenarios for doBBackTesting function
#' 
#' This is a utility function used in the \code{\link{doBBackTesting}} function
#' to define the rolling evaluation windows.
#' @param strategy Fitting-Forecasting strategy. Format: numerical vector.
#' The strategy \code{c(f, h, s)} consists in \code{f} number of years to use in 
#' fitting, \code{h} number of years to forecast and the step \code{s} to roll the
#' evaluation window.
#' @inheritParams doMortalityModels
#' @seealso \code{\link{doBBackTesting}}
#' @author Marius D. Pascariu
#' @examples 
#' y = 1900:2016
#' 
#' # Strategy
#' S = c(f = 20, h = 20, s = 1)
#' # This means that:
#' # - 20 years of data are used to fit the models;
#' # - forecast 20 years and evaluate the results;
#' # - roll the evaluation window 1 year and repeat the process.
#' 
#' # 78 scenarios are created using this years and strategy.
#' W = buildScenarios(y, S)
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
buildScenarios <- function(y, strategy = c(f = 20, h = 20, s = 2)) {
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
