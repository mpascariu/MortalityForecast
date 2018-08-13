

#' Get Accuracy Measures
#' @param data Validation set of demographic data.
#' @usage getAccuracy(object, data, x, y, data.type, what, measures, na.rm, ...)
#' @inheritParams doMortalityModels
#' @inheritParams getForecasts
#' @inheritParams computeAccuracy
#' @examples 
#' x = 0:110
#' y1 = 1980:1990
#' y2 = 1991:1993
#' h = max(y2) - max(y1)
#' 
#' D1 <- dxForecast::dxForecast.data$dx$male[paste(x), paste(y1)]
#' D2 <- dxForecast::dxForecast.data$dx$male[paste(x), paste(y2)]
#' ex <- dxForecast::dxForecast.data$ex$male
#' exogen <- ex[paste(y1)]
#' 
#' M <- doMortalityModels(data = D1, x, y1, data.type = "dx", exogen = exogen)
#' P <- doForecasts(M, h)
#' A <- getAccuracy(P, D2, x = 0:90, y2, data.type = "dx", what = "ex")
#' A
#' @export
getAccuracy <- function(object, data, x = 0:100, y = NULL,
                        data.type = c("qx", "mx", "dx", "lx"),
                        what = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"),
                        measures = c("ME", "MAE", "MAPE", "sMAPE", "MRAE", "MASE"),
                        na.rm = TRUE, ...) {
  
  O <- convertFx(x = object$x, data, In = data.type, Out = what, lx0 = 1, ...) # observed data
  H <- getForecasts(object, what)                              # forecast data
  B <- H$LC # Benchmark: Lee-Carter for now.
  
  fn <- function(X) computeAccuracy(O, X, B, x, y, measures, na.rm)
  A  <- lapply(H, fn)
  z  <- NULL
  
  for (i in 1:length(A)) {
    z <- rbind(z, A[[i]]) 
  }
  rownames(z) <- object$model.names
  
  zz <- z
  zz[, "ME"] <- abs(zz[, "ME"])
  r  <- apply(zz, 2, rank)
  s1 <- floor(rank(apply(r, 1, mean)))
  s2 <- floor(rank(apply(r, 1, median)))
  out <- list(index = what, results = z, rank = r, rankMean = s1, rankMedian = s2)
  out <- structure(class = "getAccuracy", out)
  return(out)
}


#' Print getAccuracy
#' @param x An object of the class \code{getAccuracy}
#' @inheritParams summary.getResiduals
#' @keywords internal
#' @export
print.getAccuracy <- function(x, digits = max(3L, getOption("digits") - 3L), 
                              ...) {
  cat("\nForecasting Accuracy Measures")
  cat("\nLife Table Index:", x$index, "\n\n")
  print(round(x$results, digits))
  cat("\nRanks - Best performing models in each category:\n")
  print(x$rank)
  cat("\nOverall Classification:\n")
  res <- t(data.frame(MeanRank = x$rankMean, MedianRank = x$rankMedian))
  print(res)
}


#' Get Measures of Forecast Accuracy
#' 
#' @param u Validation dataset.
#' @param u.hat In-sample forecast data.
#' @param b Benchmark forecast data. Usualy a naive or Random-Walk forecast.
#' @param x Ages to be considered in accuracy computation. It can be used to 
#' calculate the measures on a subset of the results. If \code{x = NULL} 
#' (default) the entire age-range in \code{u} is considered.
#' @param y Years to be considered in accuracy computation. Default: \code{NULL}.
#' @param measures What accurracy measure to compute? 
#' @param na.rm A logical value indicating whether NA values should be stripped 
#' before the computation proceeds.
#' @source Hyndman and Koehler, 2006
#' @keywords internal
#' @export
computeAccuracy <- function(u, u.hat, b, x = NULL, y = NULL,
                            measures = c("ME", "MAE", "MAPE", "sMAPE", "MRAE", "MASE"),
                            na.rm = TRUE){
  if (is.null(x)) x <- rownames(u)
  if (is.null(y)) y <- colnames(u)
  L1 <- rownames(u) %in% x
  L2 <- colnames(u) %in% y
    
  u.hat <- as.matrix(u.hat[L1, L2])
  u <- as.matrix(u[L1, L2])
  b <- as.matrix(b[L1, L2])
  N <- na.rm
  # ---------------------------------------------------------------
  # I. Scale-dependent measures
  E  <- u - u.hat # errors
  AE <- abs(E)    # absolute errors
  
  # 1.Mean Error
  ME  <- mean(E, na.rm = N)
  # 2.Median Error
  MdE <- median(E, na.rm = N)
  # 3.Mean Square Error
  MSE  <- mean(E^2, na.rm = N)
  # 4.Root Mean Square Error
  RMSE <- sqrt(MSE)
  # 5.Mean Absolute Error
  MAE  <- mean(AE, na.rm = N)
  # 6.Median Absolute Error
  MdAE <- median(AE, na.rm = N)
  
  # ---------------------------------------------------------------
  # II. Measures based on percentage Error
  PE <- E[u > 0]/u[u > 0] # percentage error
  APE <- abs(PE)
  # 7.Mean Absolute Percentage Error
  MAPE <- mean(APE, na.rm = N)
  # 8.Median Absolute Percentage Error
  MdAPE <- median(APE, na.rm = N)
  # 9.Root Mean Square Percentage Erorr
  RMSPE <- sqrt(mean(PE^2, na.rm = N))
  # 10.Root Median Square Percentage Error
  RMdSPE <- sqrt(median(PE^2, na.rm = N))
  
  # ----------------------------------------------
  # III. Symmetric errors
  # The MAPE and MdAPE have the disadvantage that they
  # put a heavier penalty on positive errors than on negative errors.
  # This observation led to the use of the so called "symmetric"
  # measures (Makridakis, 1993).
  
  # 11.Symmetric Mean Absolute Percentage Error
  sMAPE <- mean(2 * AE/(u + u.hat), na.rm = N)
  # 12.Symmetric Median Absolute Percentage Error
  sMdAPE <- median(2 * AE/(u + u.hat), na.rm = N)
  
  # ---------------------------------------------------------------
  # IV. Measures based on relative errors
  # An alternative way of scaling is to divide each error by
  # the error obtaned using another standard method of forecasting (benchmark method).
  
  bE  <- u - b
  RE  <- E[bE > 0]/bE[bE > 0] # relative errors.
  RAE <- abs(RE)
  # Added a very small number in b.errors to avoid divizion by zero.
  
  # 13.Mean Relative Absolute Error
  MRAE <- mean(RAE, na.rm = N)
  # 14.Median Relative Absolute Error
  MdRAE <- median(RAE, na.rm = N)
  
  # 15.Geometric Mean Relative Absolute Error
  # geometric mean function for positive values
  gm_mean = function(z) exp(mean(log(z[z > 0]), na.rm = N))
  GMRAE <- gm_mean(RAE) 
  
  # ---------------------------------------------------------------
  # V. Scaled measures
  f   <- sum(diff(t(u))) / (nrow(u) * (ncol(u) - 1)) # scaling factor
  SE  <- E / f # scaled errors
  ASE <- abs(SE)
  # 16.Mean Absolute Scaled Error
  MASE <- mean(ASE, na.rm = N)
  # 17.Median Absolute Scaled Error
  MdASE <- median(ASE, na.rm = N)
  
  # ---------------------------------------------------------------
  # Absolute errors and square errors tell the same thing. Because we do not 
  # want to receive the same information multiple times we do not export
  # MSE, RMSE, RMSPE and RMdSPE.
  out <- data.frame(ME,   MAE,  MAPE,  sMAPE,  MRAE,  MASE, # means
                    MdE, MdAE, MdAPE, sMdAPE, MdRAE, MdASE, # medians
                    MSE, RMSE, RMSPE, RMdSPE,               # squared errors
                    GMRAE)                                  # geometric mean
  out <- out[, measures]
  out <- as.matrix(out)
  return(out)
}



