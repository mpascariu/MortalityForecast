# --------------------------------------------------- #
# Author: Marius D. PASCARIU
# License: GNU General Public License v3.0
# Last update: Mon Jan 21 17:22:44 2019 
# --------------------------------------------------- #

#' Fit Multiple Stochastic Mortality Models
#' 
#' @param data A data.frame or a matrix containing mortality data 
#' with ages \code{x} as row and time \code{y} as column.
#' @param data.B A data.frame or a matrix containing mortality data for 
#' the benchmark population. This dataset is needed only in the coherent 
#' mortality models (e.g. \code{LiLee, OeppenC}). Must be the same format as 
#' in data;
#' @param x Numerical vector indicating the ages in input \code{data}. 
#' Optional. Default: \code{NULL}.
#' @param y Numerical vector indicating the years in input \code{data}. 
#' Optional. Default: \code{NULL}.
#' @param data.in Specify the type of input \code{data}. Various life table 
#' indices are accepted: \code{"qx", "mx", "dx", "lx"}.
#' @param models One or several mortality models to be estimated. 
#' The following options are available: \itemize{
#'   \item{\code{"MRW"}} -- The Multivariate Random-Walk (w/o Drift);
#'   \item{\code{"MRWD"}} -- The Multivariate Random-Walk with Drift;
#'   \item{\code{"LeeCarter"}} -- The Lee-Carter Mortality Model;
#'   \item{\code{"LiLee"}} -- The Li-Lee Mortality Model;
#'   \item{\code{"HyndmanUllah"}} -- The Hyndman-Ullah Mortality Model;
#'   \item{\code{"RenshawHaberman"}} -- The Renshaw-Haberman Mortality Model;
#'   \item{\code{"Oeppen"}} -- The Oeppen Mortality Model;
#'   \item{\code{"OeppenC"}} -- The Coherent Oeppen Mortality Model;
#'   \item{\code{"MEM2"}} -- The Maximum Entropy Mortality Model of order 2;
#'   \item{\code{"MEM3"}} -- The Maximum Entropy Mortality Model of order 3;
#'   \item{\code{"MEM4"}} -- The Maximum Entropy Mortality Model of order 4;
#'   \item{\code{"MEM5"}} -- The Maximum Entropy Mortality Model of order 5;
#'   \item{\code{"MEM6"}} -- The Maximum Entropy Mortality Model of order 6;
#'   \item{\code{"MEM7"}} -- The Maximum Entropy Mortality Model of order 7;
#' }
#' @param verbose A logical value. Set \code{verbose = FALSE} to silent 
#' the process that take place inside the function and avoid progress messages.
#' @param ... Arguments to be passed to or from other methods.
#' @examples 
#' x  <- 0:100      # Ages
#' y  <- 2005:2016  # Years
#' h  <- 16         # forecasting horizon
#' MM <- c("MRWD", "LeeCarter", "LiLee", "HyndmanUllah", 
#'         "Oeppen", "OeppenC", "MEM5", "RenshawHaberman") # mortality models
#' D  <- HMD_male$dx$GBRTENW[paste(x), paste(y)]  # data
#' B  <- HMD_female$dx$GBRTENW[paste(x), paste(y)]  # benchmark population
#' 
#' # Note: We are fitting various mortality model to E&W males data. Most of them 
#' # are single population model i.e. the estimates are resulted only from data 
#' # specific to that population. However, the coherent models like "OeppenC" or 
#' # "LiLee" require additional data. In these cases, here, female data is used.
#' 
#' 
#' M <- do.MortalityModels(data = D, 
#'                         data.B = B,
#'                         x = x, 
#'                         y = y, 
#'                         data.in = "dx", 
#'                         models = MM)
#' 
#' P <- do.MortalityForecasts(object = M, 
#'                            h = h, 
#'                            level = 95, 
#'                            jumpchoice = "actual")
#' 
#' oex <- get.Observed(M, data.out = "ex")
#' fex <- get.Fitted(M, data.out = "ex")
#' rex <- get.Residuals(M, data.out = "ex")
#' pex <- get.Forecasts(P, data.out = "ex")
#' @export
do.MortalityModels <- function(data,
                               data.B = NULL,
                               x = NULL, 
                               y = NULL, 
                               data.in = c("qx", "mx", "dx", "lx"),
                               models = c("MRWD"),
                               verbose = TRUE, 
                               ...) {
  
  data.in <- match.arg(data.in)
  input <- as.list(environment())
  call  <- match.call()
  M     <- which_model_is_in(models)
  x <- x %||% 1:nrow(data)
  y <- y %||% 1:ncol(data)
  
  
  ## SINGLE POPULATION MODELS
  # ------------------------------------------
  mx.data <- convertFx(x, data, from = data.in, to = "mx", lx0 = 1, ...)
  dx.data <- convertFx(x, data, from = data.in, to = "dx", lx0 = 1, ...)

  # The Naive model - Multivariate Random-Walk
  if (M["MRW"]) {
    MRW <- model.MRW(data = log(mx.data), x = x, y = y, include.drift = FALSE)
  }
  # Random Walk with drift
  if (M["MRWD"]) {
    MRWD <- model.MRW(data = log(mx.data), x = x, y = y, include.drift = TRUE)
  }
  # Lee-Carter (1992)
  if (M["LeeCarter"]) {
    LeeCarter <- model.LeeCarter(data = mx.data, x = x, y = y, verbose = FALSE)
  }
  # Hyndman-Ullah (1992)
  if (M["HyndmanUllah"]) {
    HyndmanUllah <- model.HyndmanUllah(data = mx.data, x = x, y = y, 
                                       verbose = FALSE)
  }
  # Renshaw-Haberman (2006)
  if (M["RenshawHaberman"]) {
    RenshawHaberman <- model.RenshawHaberman(data = mx.data, x = x, y = y, 
                                             verbose = FALSE)
  }
  # Oeppen (2008)
  if (M["Oeppen"]) Oeppen <- model.Oeppen(data = dx.data, x = x, y = y, 
                                                 verbose = FALSE)
  # Maximum Entropy Mortality Models - PLC (2019)
  if (M["MEM2"]) MEM2 <- model.MEM(data = dx.data, x = x, y = y, 
                                          n = 2, verbose = FALSE)
  if (M["MEM3"]) MEM3 <- model.MEM(data = dx.data, x = x, y = y, 
                                          n = 3, verbose = FALSE)
  if (M["MEM4"]) MEM4 <- model.MEM(data = dx.data, x = x, y = y, 
                                          n = 4, verbose = FALSE)
  if (M["MEM5"]) MEM5 <- model.MEM(data = dx.data, x = x, y = y, 
                                          n = 5, verbose = FALSE)
  if (M["MEM6"]) MEM6 <- model.MEM(data = dx.data, x = x, y = y, 
                                          n = 6, verbose = FALSE)
  if (M["MEM7"]) MEM7 <- model.MEM(data = dx.data, x = x, y = y, 
                                          n = 7, verbose = FALSE)
  
  ## COHERENT MODELS
  # ------------------------------------------
  if (!is.null(data.B)) {
    mx.data.B <- convertFx(x, data.B, from = data.in, to = "mx", lx0 = 1, ...)
    dx.data.B <- convertFx(x, data.B, from = data.in, to = "dx", lx0 = 1, ...)
  
    if (M["LiLee"]) {
      LiLee <- model.LiLee(data = mx.data, data.B = mx.data.B, 
                           x = x, y = y, verbose = FALSE)
    }
    if (M["OeppenC"]) {
      OeppenC <- model.OeppenC(data = dx.data, data.B = dx.data.B, 
                               x = x, y = y, verbose = FALSE)
    }
    remove(data.B, dx.data.B, mx.data.B)
  }
  
  remove(data, data.in, models, dx.data, mx.data, M)
  out <- as.list(environment())
  out <- structure(class = "MortalityModels", out)
  return(out)
}


#' Determine which model is included in input
#' @param X Models in input
#' @param M All models
#' @keywords internal
which_model_is_in <- function(X,
                              M = c("MRW", 
                                    "MRWD", 
                                    "LeeCarter", 
                                    "HyndmanUllah", 
                                    "RenshawHaberman", 
                                    "Oeppen",
                                    "MEM2", 
                                    "MEM3", 
                                    "MEM4", 
                                    "MEM5", 
                                    "MEM6", 
                                    "MEM7",
                                    "LiLee", 
                                    "OeppenC")) {
   out <- M %in% X
   names(out) <- M
   return(out)
}


#' Print function for do.MortalityModels
#' @param x An object of class \code{"MortalityModels"}
#' @inheritParams do.MortalityModels
#' @keywords internal
#' @export
print.MortalityModels <- function(x, ...) {
  cat("Stochastic Mortality Models - FIT")
  cat("\nCall : "); print(x$call)
  cat("\nModels in fit:", paste(x$input$models, collapse = ", "))
  cat("\nAges  in fit :", paste(range(x$x), collapse = " - "))
  cat("\nYears in fit :", paste(range(x$y), collapse = " - "))
  cat("\n")
}

