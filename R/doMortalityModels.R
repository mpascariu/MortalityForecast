

#' Wrapper for Function Fitting Stochastic Mortality Models
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
#' D  <- MortalityForecast.data$dx[paste(x), paste(y)]
#' MM <- c("LC", "FDM", "CoDa", "M6")
#' 
#' M <- doMortalityModels(data = D, x, y, data.in = "dx", models = MM)
#' 
#' oex <- getObserved(M, data.out = "ex")
#' fex <- getFitted(M, data.out = "ex")
#' rex <- getResiduals(M, data.out = "ex")
#' @export
doMortalityModels <- function(data, x = NULL, y = NULL, 
                              data.in = c("qx", "mx", "dx", "lx"),
                              models = c("MRWD","LC"),
                              verbose = TRUE, ...) {
  
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
  if ("MRW" %in% models) MRW <- MRW(data = log(mx.data), x, y, include.drift = FALSE)
  # Random Walk with drift
  if ("MRWD" %in% models) MRWD <- MRW(data = log(mx.data), x, y, include.drift = TRUE)
  # LC (1992)
  if ("LC" %in% models) LC <- LC(data = mx.data, x, y)
  # FDM (1992)
  if ("FDM" %in% models) FDM <- FDM(data = mx.data, x, y)
  # Plat Model (2009)
  if ("PLAT" %in% models) PLAT <- PLAT(data = mx.data, x, y)
  # CoDa-LC (2008)
  if ("CoDa" %in% models) CoDa <- coda(data = dx.data, x, y)
  # Maximum Entropy Mortality Models - PLC (2018)
  if ("M4" %in% models)  M4 <- MEM(data = dx.data, x, y, n = 4)
  if ("M5" %in% models)  M5 <- MEM(data = dx.data, x, y, n = 5)
  if ("M6" %in% models)  M6 <- MEM(data = dx.data, x, y, n = 6)
  if ("M7" %in% models)  M7 <- MEM(data = dx.data, x, y, n = 7)

  remove(data, data.in, models, dx.data, mx.data)
  out <- as.list(environment())
  out <- structure(class = "MortalityModels", out)
  return(out)
}


#' Get Fitted Values
#' @param object An object of class \code{MortalityModels}.
#' @inheritParams getAccuracy
#' @export
getFitted <- function(object, 
                      data.out = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"),
                      ...) {
  data.out <- match.arg(data.out)
  Mn   <- object$input$models # Model names
  x    <- object$x
  
  DX <- list()
  for (i in 1:length(Mn)) {
    M <- with(object, get(Mn[i]))
    
    if (Mn[i] %in% c("MRW", "MRWD", "LC", "PLAT")) {
      mx <- exp(fitted(M))
      dx <- convertFx(x, mx, from = "mx", to = "dx", lx0 = 1, ...)
      
    } else if (Mn[i] %in% c("FDM")) {
      mx <- exp(fitted(M)$y)
      dx <- convertFx(x, mx, from = "mx", to = "dx", lx0 = 1, ...)
      
    } else {
      dx <- fitted(M)
    }
    DX[[i]] <- dx
  }
  
  fn  <- function(Z) convertFx(x, Z, from = "dx", to = data.out, lx0 = 1, ...)
  out <- lapply(DX, fn)
  names(out) <- Mn
  out <- structure(class = "getFitted", out)
  return(out)
}


#' Get Observed Values
#' @inheritParams getFitted
#' @export
getObserved <- function(object, 
                        data.out = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"),
                        ...) {
  x        <- object$x
  data     <- object$input$data
  data.in  <- object$input$data.in 
  data.out <- match.arg(data.out)
  
  if (data.in == data.out) {
    out <- data
  } else {
    out <- convertFx(x, data, from = data.in, to = data.out, lx0 = 1, ...)
  }
  return(out)
}






