

#' Generic for Computing Accuracy Measure from the \code{doBackTesting} and 
#' \code{doBBackTesting} objects
#' @param object An object of the class \code{doBackTesting} or \code{doBBackTesting}.
#' @inheritParams doMortalityModels
#' @keywords internal
#' @author Marius D. Pascariu
#' @examples 
#' # For examples go to ?doBackTesting and ?doBBackTesting
#' @export
evalAccuracy = function(object, ...)
  UseMethod("evalAccuracy")


#' Get Accuracy Measures from \code{doBackTesting}
#' @param object An object of class \code{doBackTesting}.
#' @param data.out Specify the type of data to be returned in output. 
#' Various life table indices are accepted: 
#' \code{"qx", "mx", "dx", "lx", "Lx", "Tx", "ex"}.
#' @inheritParams doMortalityModels
#' @inheritParams computeAccuracy
#' @inherit computeAccuracy details references
#' @seealso \code{\link{doBackTesting}}
#' @author Marius D. Pascariu
#' @examples 
#' # For examples go to ?doBackTesting
#' @export
evalAccuracy.doBackTesting <- function(object, 
                        data.out = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"),
                        measures = NULL,
                        ...) {
  # Data
  x <- object$input$x
  data.in  <- object$input$data.in
  data.out <- match.arg(data.out)
  validation.set <- object$Datasets$validation.set
  
  # Get the observed data in the required format (e.g. e[x], q[x] etc)
  O  <- convertFx(x = x, 
                  data = validation.set, 
                  from = data.in, 
                  to = data.out, 
                  lx0 = 1, ...)
  # Get the forecast results in the same format
  # H is a list of matrices, each matrix coresponding to 1 model
  H  <- getForecasts(object$Forecast, data.out) # forecast data
  B  <- H[[1]]                                  # Benchmark model
  Mn <- object$input$models                     # Model names
  
  # Compute the accuracy measures for all models/forecasts
  fn <- function(X) computeAccuracy(O, X, B, measures)
  A  <- lapply(H, fn)  # this is a list too w accuracy measures
  z  <- NULL
  
  # Create a tibble and exit
  for (i in 1:length(A)) {
    z <- rbind(z, A[[i]]) 
  }
  zt  <- add_column(as.tibble(z), Scenario = "Total", Model = Mn, 
                    LifeTableIndex = data.out, .before = T)
  zt
}


#' Get Accuracy Measures from \code{doBBackTesting}
#' @param object An object of class \code{doBBackTesting}.
#' @inheritParams evalAccuracy.doBackTesting
#' @inherit evalAccuracy.doBackTesting details references
#' @seealso \code{\link{doBBackTesting}}
#' @author Marius D. Pascariu
#' @examples 
#' # For examples go to ?doBBackTesting
#' @export
evalAccuracy.doBBackTesting <- function(object,
                                       data.out = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"),
                                       measures = NULL,
                                       ...) {
  data.out <- match.arg(data.out)
  ns <- nrow(object$scenarios)  # no. of scenarios
  A  <- tibble()
  AA <- 0
  N  <- -c(1:3)
  
  for (s in 1:ns) {
    As <- evalAccuracy(object = object$results[[s]], data.out, measures, ...)
    As$Scenario <- s
    A  <- rbind(A, As)
    AA <- AA + As[, N]/ns  # compute mean values over all scenarios
  }
  
  As[, N] <- AA
  As$Scenario <- "Total"
  out <- rbind(As, A)
  return(out)
}


# #' Print evalAccuracy
# #' @param x An object of the class \code{evalAccuracy}
# #' @inheritParams summary.getResiduals
# #' @keywords internal
# #' @export
# print.evalAccuracy <- function(x, digits = max(3L, getOption("digits") - 3L), 
#                               ...) {
#   cat("\nForecasting Accuracy Measures")
#   cat("\nLife Table Index:", x$index, "\n\n")
#   print(round(x$results, digits))
#   cat("\nRanks - Best performing models in each category:\n")
#   print(x$rank)
#   cat("\nGeneral Classification:\n")
#   # res <- t(data.frame(MeanRank = x$rankMean, MedianRank = x$rankMedian))
#   print(x$GC)
# }


#' Get Measures of Forecast Accuracy
#' 
#' @param u Validation dataset.
#' @param u.hat Out-sample forecast data.
#' @param b Benchmark forecast data. Usualy a naive or Random-Walk w drift forecast.
#' @param measures What accurracy measure to compute? Various alternatives are 
#' available, \itemize{
#'  \item{Mean error measures: } \code{"ME", "MAE", "MAPE", "sMAPE", "sMRAE", "MASE"};
#'  \item{Median error measures: } \code{"MdE", "MdAE", "MdAPE", "sMdAPE", "sMdRAE", "MdASE"};
#'  \item{Squared error measures: } \code{"MSE", "RMSE", "RMSPE", "RMdSPE"};
#'  \item{Geometric mean measure for positive errors: } \code{"GMRAE"}.}
#'  If \code{measures = NULL} all the measures will be computed.
#' @param xa Ages to be considered in model accuracy evaluation. It can be used 
#' to calculate the measures on a subset of the results. If \code{xa = NULL} 
#' (default) the entire age-range in input is considered.
#' @param ya Years to be considered in accuracy computation. Default: \code{ya = NULL}.
#' @param na.rm A logical value indicating whether NA values should be stripped 
#' before the computation proceeds.
#' @details See \insertCite{hyndman2006;textual}{MortalityForecast} for a 
#' comprehensive discussion of the accuracy measures.
#' @references \insertAllCited{}
#' @author Marius D. Pascariu
#' @keywords internal
computeAccuracy <- function(u, u.hat, b, measures = NULL, 
                            xa = NULL, ya = NULL, na.rm = TRUE){
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
  MAPE <- mean(100 * APE, na.rm = N)
  # 8.Median Absolute Percentage Error
  MdAPE <- median(100 * APE, na.rm = N)
  # 9.Root Mean Square Percentage Erorr
  PEs <- (100 * PE)^2
  RMSPE <- sqrt(mean(PEs, na.rm = N))
  # 10.Root Median Square Percentage Error
  RMdSPE <- sqrt(median(PEs, na.rm = N))
  
  # ----------------------------------------------
  # III. Symmetric errors
  # The MAPE and MdAPE have the disadvantage that they
  # put a heavier penalty on positive errors than on negative errors.
  # This observation led to the use of the so called "symmetric"
  # measures (Makridakis, 1993).
  
  # 11.Symmetric Mean Absolute Percentage Error
  sape <- 200 * AE/(u + u.hat)
  sMAPE <- mean(sape, na.rm = N)
  # 12.Symmetric Median Absolute Percentage Error
  sMdAPE <- median(sape, na.rm = N)
  
  # ---------------------------------------------------------------
  # IV. Measures based on relative errors
  # An alternative way of scaling is to divide each error by
  # the error obtaned using another standard method of forecasting (benchmark method).
  
  bE  <- u - b       # benchmark errors
  bAE <- abs(bE)
  # RAE <- AE/bAE      # relative absolute errors.
  sRAE <- 200 * AE/(AE + bAE)  # relative absolute errors.
  
  # 13.Mean Relative Absolute Error
  sMRAE <- mean(sRAE, na.rm = N)
  # 14.Median Relative Absolute Error
  sMdRAE <- median(sRAE, na.rm = N)
  
  # 15.Geometric Mean Relative Absolute Error
  # geometric mean function for positive values
  # gm_mean = function(z) exp(mean(log(z[z > 0]), na.rm = N))
  # GMRAE <- gm_mean(RAE) 
  
  # ---------------------------------------------------------------
  # V. Scaled measures
  meanT <- function(z) mean(z, na.rm = TRUE) # mean() function
  scale <- apply(abs(t(diff(t(u)))), 1, meanT)  # compute a scale factor for each time series
  SE    <- sweep(E, 1, scale, FUN = "/")     # scaled errors
  ASE   <- abs(SE)
  # 16.Mean Absolute Scaled Error
  MASE <- mean(ASE, na.rm = N)
  # 17.Median Absolute Scaled Error
  MdASE <- median(ASE, na.rm = N)
  
  # ---------------------------------------------------------------
  # Absolute errors and square errors tell the same thing. Because we do not 
  # want to receive the same information multiple times we do not export
  # MSE, RMSE, RMSPE and RMdSPE.
  out <- data.frame(ME,   MAE,  MAPE,  sMAPE,  sMRAE,  MASE, # means
                    MdE, MdAE, MdAPE, sMdAPE, sMdRAE, MdASE, # medians
                    MSE, RMSE, RMSPE, RMdSPE)                # squared errors
  
  if(!is.null(measures)) out <- out[, measures]
  out <- as.matrix(out)
  return(out)
}



