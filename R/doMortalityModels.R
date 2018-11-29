# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Tue Nov 27 15:50:43 2018
# --------------------------------------------------- #


#' Fit Multiple Stochastic Mortality Models
#' 
#' @param data A data.frame or a matrix containing mortality data 
#' with ages \code{x} as row and time \code{y} as column.
#' @param x Numerical vector indicating the ages in input \code{data}. Optional.
#' Default: \code{NULL}.
#' @param y Numerical vector indicating the years in input \code{data}. Optional.
#' Default: \code{NULL}.
#' @param data.in Specify the type of input \code{data}. Various life table 
#' indices are accepted: \code{"qx", "mx", "dx", "lx"}.
#' @param models Mortality models to be evaluated.
#' @param verbose A logical value. Set \code{verbose = FALSE} to silent 
#' the process that take place inside the function and avoid progress messages.
#' @param ... Arguments to be passed to or from other methods.
#' @examples 
#' x  <- 0:100
#' y  <- 2005:2016
#' D  <- HMD_male$dx$GBRTENW[paste(x), paste(y)]
#' MM <- c("MRWD", "HyndmanUllah", "Oeppen")
#' 
#' M <- doMortalityModels(data = D, x, y, data.in = "dx", models = MM)
#' 
#' oex <- get.Observed(M, data.out = "ex")
#' fex <- get.Fitted(M, data.out = "ex")
#' rex <- get.Residuals(M, data.out = "ex")
#' @export
doMortalityModels <- function(data, 
                              x = NULL, 
                              y = NULL, 
                              data.in = c("qx", "mx", "dx", "lx"),
                              models = c("MRWD"),
                              verbose = TRUE, 
                              ...) {
  
  data.in <- match.arg(data.in)
  input   <- as.list(environment())
  call    <- match.call()
  x <- x %||% 1:nrow(data)
  y <- y %||% 1:ncol(data)
  
  mx.data <- convertFx(x, data, from = data.in, to = "mx", lx0 = 1, ...)
  dx.data <- convertFx(x, data, from = data.in, to = "dx", lx0 = 1, ...)
  # qx.data <- convertFx(x, data, from = data.in, to = "qx", lx0 = 1) # not needed now
  # lx.data <- convertFx(x, data, from = data.in, to = "lx", lx0 = 1)

  # The Naive model - Multivariate Random-Walk
  if ("MRW" %in% models) {
    MRW <- model.MRW(data = log(mx.data), x = x, y = y, include.drift = FALSE)
  }
  # Random Walk with drift
  if ("MRWD" %in% models) {
    MRWD <- model.MRW(data = log(mx.data), x = x, y = y, include.drift = TRUE)
  }
  # Lee-Carter (1992)
  if ("LC" %in% models) {
    LC <- LC(data = mx.data, x = x, y = y, link = "log")
  }
  if ("LeeCarter" %in% models) {
    LeeCarter <- model.LeeCarter(data = mx.data, x = x, y = y)
  }
  # Hyndman-Ullah (1992)
  if ("HyndmanUllah" %in% models) {
    HyndmanUllah <- model.HyndmanUllah(data = mx.data, x = x, y = y)
  }
  # Plat (2009)
  if ("PLAT" %in% models) PLAT <- PLAT(data = mx.data, x = x, y = y)
  # Oeppen (2008)
  if ("Oeppen" %in% models) Oeppen <- model.Oeppen(data = dx.data, x = x, y = y)
  # Maximum Entropy Mortality Models - PLC (2018)
  if ("MEM2" %in% models)  MEM2 <- model.MEM(data = dx.data, x = x, y = y, n = 2)
  if ("MEM3" %in% models)  MEM3 <- model.MEM(data = dx.data, x = x, y = y, n = 3)
  if ("MEM4" %in% models)  MEM4 <- model.MEM(data = dx.data, x = x, y = y, n = 4)
  if ("MEM5" %in% models)  MEM5 <- model.MEM(data = dx.data, x = x, y = y, n = 5)
  if ("MEM6" %in% models)  MEM6 <- model.MEM(data = dx.data, x = x, y = y, n = 6)
  if ("MEM7" %in% models)  MEM7 <- model.MEM(data = dx.data, x = x, y = y, n = 7)

  remove(data, data.in, models, dx.data, mx.data)
  out <- as.list(environment())
  out <- structure(class = "MortalityModels", out)
  return(out)
}



#' Print function for doMortalityModels
#' @param x An object of class \code{"MortalityModels"}
#' @inheritParams doMortalityModels
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

