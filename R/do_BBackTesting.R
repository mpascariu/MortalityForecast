# --------------------------------------------------- #
# Author: Marius D. PASCARIU
# License: GNU General Public License v3.0
# Last update: Tue Jan 22 10:03:41 2019 
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
#' M  <- c("MRWD", "LeeCarter", "HyndmanUllah",
#'         "Oeppen", "MEM4", "MEM5")
#' 
#' BB <- do.BBackTesting(data = dx, x = x, y = y,
#'                       data.in = "dx", 
#'                       models = M,
#'                       strategy = c(20, 20, 1))
#' 
#' A <- evalAccuracy(BB, data.out = "ex")
#' A
#' 
#' R <- do.Ranking(A)
#' R
#' @export
do.BBackTesting <- function(data, 
                            data.B = NULL,
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
  
  # Do Back-testing in parallel -----------------------------
  nc <- nrow(S) # no. of cases
  w  <- num_workers(res = nc)
  cl <- makeCluster(w)
  # clusterEvalQ(cl = cl,library(MortalityForecast)) # needed only if the code 
  # is runned outside the package
  registerDoParallel(cl, cores = w)
  
  i <- NULL # hack R CMD Check note
  B <- foreach(i = seq_len(nc)) %dopar% {
    yf <- S[[i, "fit"]]
    yh <- S[[i, "forecast"]]
    Y <- c(yf, yh)
    do.BackTesting(data = data[, paste(Y)],
                   data.B = data.B[, paste(Y)],
                   x = x, 
                   y.fit = yf, 
                   y.for = yh, 
                   data.in = data.in,
                   models = models, 
                   level = level, 
                   jumpchoice = jumpchoice, 
                   verbose = FALSE)
  }
  stopCluster(cl)
  names(B) <- paste0("S", 1:nc)
  
  # Exit -----------------------------------------------------
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
  S <- as.numeric(strategy)
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


#' Determine the number of cores to be used
#' @param n.cores Number of cores to be used in parallel computing;
#' @param res Restrict the number of cores. If there are only 2 operations to 
#' performe even if more cores are available 2 cores are used.
#' @keywords internal
num_workers <- function(n.cores = NULL, res = NULL) {
  if (is.null(n.cores)) {
    
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (nzchar(chk) && chk == "TRUE") {
      # use 2 cores in CRAN/Travis/AppVeyor checks because those are the rules
      n.cores <- 2L
    } else {
      # use all cores in devtools::test()
      n.cores <- detectCores()
    }
  }
  
  if (is.null(res)) res <- 64L
  
  X <- as.integer(min(n.cores, res))
  return(X)
}
