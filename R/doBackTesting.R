

#' Performe back-testing
#' @param Tdata Test data
#' @param by The accuracy measures are available in different output formats.
#' Option 1: \code{by = NULL} the agregated measure over ages and and years will 
#' be return. This is the default option. 
#' Option 2: \code{by = "x"} measure for all ages (agregated over time) are returned.
#' Option 3: \code{by = "y"} measures for all years (aggregated over ages) are returned.
#' @inheritParams doMortalityModels
#' @inheritParams getForecasts
#' @inheritParams computeAccuracy
#' @export
doBackTesting <- function(Tdata, object, 
                          data.type = c("qx", "mx", "dx"),
                          type = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"),
                          by = NULL, na.rm = TRUE, ...) {
  x <- object$x
  O <- convertFx(x, Tdata, In = data.type, Out = type, lx0 = 1)
  H <- getForecasts(object, type)
  B <- H$LC # Benchmark: Lee-Carter for now.
  
  fn <- function(X) core.Accuracy(u = O, u.hat = X, b = B, by, na.rm)
  A  <- lapply(H, fn)
  out <- A
  
  if (is.null(by)) {
    z = NULL
    for (i in 1:length(A)) {
      z <- rbind(z, A[[i]]) 
    }
    rownames(z) <- object$model.names
    
    zz <- z
    zz[, "ME"] <- abs(zz[, "ME"])
    r <- apply(zz, 2, rank)
    s1 <- sort(round(apply(r, 1, mean), 1))
    s2 <- sort(floor(apply(r, 1, median)))
    out <- list(measures = z, rank = r, 
                rankMean = s1, rankMedian = s2)
  }
  return(out)
}


#' Generate different types of output 
#' 
#' @inheritParams doBackTesting
#' @inheritParams computeAccuracy
#' @keywords internal
core.Accuracy <- function(u, u.hat, b, by = NULL, na.rm = TRUE) {
  A <- NULL
  
  if (is.null(by)) {
    A <- computeAccuracy(u, u.hat, b)
  } else {
    
    if (by == "y") {
      for (i in 1:ncol(u.hat)) {
        a <- computeAccuracy(u[, i], u.hat[, i], b[, i])
        A <- rbind(A, a)
      }
      rownames(A) <- colnames(u.hat)
    }
    
    if (by == "x") {
      for (j in 1:nrow(u.hat)) {
        a <- computeAccuracy(u[j, ], u.hat[j, ], b[j, ])
        A <- rbind(A, a)
      }
      rownames(A) <- rownames(u.hat)
    }
  }
  return(A)
}

#' Get Measures of Forecast Accuracy
#' 
#' @param u Observed data
#' @param u.hat In-sample Forecast data
#' @param b Benchmark forecast data. Usualy a naive or rwd forecast.
#' @param what What accurracy measure to compute? If \code{what = NULL} all 
#' 15 measures will be generated.
#' @param na.rm A logical value indicating whether NA values should be stripped 
#' before the computation proceeds.
#' @source Hyndman and Koehler, 2006
#' @keywords internal
#' @export
computeAccuracy <- function(u, u.hat, b, na.rm = TRUE){

  u.hat <- as.matrix(u.hat)
  u <- as.matrix(u)
  b <- as.matrix(b)
  N <- na.rm
  # ---------------------------------------------------------------
  # I. Scale-dependent measures
  E  <- u - u.hat # errors
  AE <- abs(E)    # absolute errors
  
  # 1.Mean Error
  ME  <- mean(E, na.rm = N)
  MdE <- median(E, na.rm = N)
  # 2.Mean Square Error
  MSE  <- mean(E^2, na.rm = N)
  # 3.Root Mean Square Error
  RMSE <- sqrt(MSE)
  # 4.Mean Absolute Error
  MAE  <- mean(AE, na.rm = N)
  # 5.Median Absolute Error
  MdAE <- median(AE, na.rm = N)
  
  # ---------------------------------------------------------------
  # II. Measures based on percentage Error
  PE <- E[u > 0]/u[u > 0] # percentage error
  APE <- abs(PE)
  # 6.Mean Absolute Percentage Error
  MAPE <- mean(APE, na.rm = N)
  # 7.Median Absolute Percentage Error
  MdAPE <- median(APE, na.rm = N)
  # 8.Root Mean Square Percentage Erorr
  RMSPE <- sqrt(mean(PE^2, na.rm = N))
  # 9.Root Median Square Percentage Error
  RMdSPE <- sqrt(median(PE^2, na.rm = N))
  
  # ----------------------------------------------
  # III. Symmetric errors
  # The MAPE and MdAPE have the disadvantage that they
  # put a heavier penalty on positive errors than on negative errors.
  # This observation led to the use of the so called "symmetric"
  # measures (Makridakis, 1993).
  
  # 10.Symmetric Mean Absolute Percentage Error
  sMAPE <- mean(2 * AE/(u + u.hat), na.rm = N)
  # 11.Symmetric Median Absolute Percentage Error
  sMdAPE <- median(2 * AE/(u + u.hat), na.rm = N)
  
  # ---------------------------------------------------------------
  # IV. Measures based on relative errors
  # An alternative way of scaling is to divide each error by
  # the error obtaned using another standard method of forecasting (benchmark method).
  
  bE  <- u - b
  RE  <- E[bE > 0]/bE[bE > 0] # relative errors.
  RAE <- abs(RE)
  # Added a very small number in b.errors to avoid divizion by zero.
  
  # 12.Mean Relative Absolute Error
  MRAE <- mean(RAE, na.rm = N)
  # 13.Median Relative Absolute Error
  MdRAE <- median(RAE, na.rm = N)
  
  # 14.Geometric Mean Relative Absolute Error
  # geometric mean function for positive values
  gm_mean = function(z) exp(mean(log(z[z > 0]), na.rm = N))
  GMRAE <- gm_mean(RAE) 
  
  # ---------------------------------------------------------------
  # V. Scaled measures
  f   <- sum(diff(t(u))) / (nrow(u) * (ncol(u) - 1)) # scaling factor
  SE  <- E / f # scaled errors
  ASE <- abs(SE)
  # 15.Mean Absolute Scaled Error
  MASE <- mean(ASE, na.rm = N)
  # 16.Median Absolute Scaled Error
  MdASE <- median(ASE, na.rm = N)
  
  # ---------------------------------------------------------------
  # Absolute errors and square errors tell the same thing. Because we do not 
  # want to receive the same information multiple times we do not export
  # MSE, RMSE, RMSPE and RMdSPE.
  # Ignore for now GMRAE too.
  out <- data.frame(ME,   MAE,  MAPE,  sMAPE,  MRAE,  MASE, # means
                    MdE, MdAE, MdAPE, sMdAPE, MdRAE, MdASE) # medians
  out <- as.matrix(out)
  return(out)
}



