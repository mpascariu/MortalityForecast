

#' Get Accuracy Measures
#' @param object An object of class \code{doForecasts}.
#' @param data Validation set of demographic data.
#' @param data.out Specify the type of data to be returned in output. 
#' Various life table indices are accepted: 
#' \code{"qx", "mx", "dx", "lx", "Lx", "Tx", "ex"}.
#' @inheritParams doMortalityModels
#' @inheritParams computeAccuracy
#' @inherit computeAccuracy details references
#' @examples 
#' x = 0:100             # Ages
#' y1 = 1980:1995        # Training period
#' y2 = 1996:2005        # Validation period
#' h = max(y2) - max(y1) # Forecasting horizon
#' 
#' D1 <- MortalityForecast.data$dx[paste(x), paste(y1)]
#' D2 <- MortalityForecast.data$dx[paste(x), paste(y2)]
#' 
#' MM = c("MRWD", "LC", "CoDa")
#' M <- doMortalityModels(data = D1, x, y1, data.in = "dx", models = MM)
#' P <- doForecasts(M, h)
#' A <- getAccuracy(P, D2, xa = 0:95, data.in = "dx", data.out = "qx")
#' 
#' A
#' plot(A)
#' @export
getAccuracy <- function(object, data, 
                        data.in = c("qx", "mx", "dx", "lx"),
                        data.out = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"),
                        measures = c("ME", "MAE", "MAPE", "sMAPE", "MRAE", "MASE"),
                        xa = NULL, ya = NULL,
                        na.rm = TRUE, ...) {
  
  data.in  <- match.arg(data.in)
  data.out <- match.arg(data.out)
  
  O  <- convertFx(x = object$x, data, from = data.in, to = data.out, lx0 = 1, ...) # observed data
  H  <- getForecasts(object, data.out)                              # forecast data
  B  <- H[[1]] # Benchmark model
  Mn <- object$input$object$input$models # Model names
  
  fn <- function(X) computeAccuracy(O, X, B, xa, ya, measures, na.rm)
  A  <- lapply(H, fn)
  z  <- NULL
  
  for (i in 1:length(A)) {
    z <- rbind(z, A[[i]]) 
  }
  rownames(z) <- Mn
  R  <- doRanking(z)
  
  out <- list(index = data.out, results = z, rank = R$rank, GC = R$GC)
  out <- structure(class = "getAccuracy", out)
  return(out)
}


#' Rank models based on the accuracy results
#' @param z Table containing accuracy measures.
#' @export
doRanking <- function(z) {
  zz <- z
  zz[, "ME"] <- abs(zz[, "ME"])
  r  <- floor(apply(zz, 2, rank))
  s2 <- floor(rank(apply(r, 1, median)))
  out <- list(rank = r, GC = s2)
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
  cat("\nGeneral Classification:\n")
  # res <- t(data.frame(MeanRank = x$rankMean, MedianRank = x$rankMedian))
  print(x$GC)
}


#' Get Measures of Forecast Accuracy
#' 
#' @param u Validation dataset.
#' @param u.hat In-sample forecast data.
#' @param b Benchmark forecast data. Usualy a naive or Random-Walk forecast.
#' @param xa Ages to be considered in model accuracy evaluation. It can be used to 
#' calculate the measures on a subset of the results. If \code{x = NULL} 
#' (default) the entire age-range in \code{u} is considered.
#' @param ya Years to be considered in accuracy computation. Default: \code{NULL}.
#' @param measures What accurracy measure to compute? Various alternatives are 
#' available, \itemize{
#' \item{Mean error measures: } \code{"ME", "MAE", "MAPE", "sMAPE", "MRAE", "MASE"};
#' \item{Median error measures: } \code{"MdE", "MdAE", "MdAPE", "sMdAPE", "MdRAE", "MdASE"};
#' \item{Squared error measures: } \code{"MSE", "RMSE", "RMSPE", "RMdSPE"};
#' \item{Geometric mean measure for positive errors: } \code{"GMRAE"}.}
#' @param na.rm A logical value indicating whether NA values should be stripped 
#' before the computation proceeds.
#' @details See \insertCite{hyndman2006;textual}{MortalityForecast} for a 
#' comprehensive discussion of the accuracy measures.
#' @references \insertAllCited{}
#' @keywords internal
#' @export
computeAccuracy <- function(u, u.hat, b, xa = NULL, ya = NULL,
                            measures, na.rm = TRUE){
  if (is.null(xa)) xa <- rownames(u)
  if (is.null(ya)) ya <- colnames(u)
  L1 <- rownames(u) %in% xa
  L2 <- colnames(u) %in% ya
    
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
  
  bE  <- u - b        # benchmark errors
  RAE <- AE/abs(bE)  # relative absolute errors.
  
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



