

#' Wrapper for function dealing with mortality models
#' 
#' @param data Mortality data
#' @param x Vector of ages.
#' @param y Vector of years.
#' @param data.type Type of data in 'data' argument. Options: 
#' \code{"qx", "mx", "dx", "lx"}.
#' @param models Mortality models to be evaluated.
#' @param ... Arguments to be passed to or from other methods.
#' @inheritParams dxForecast::lenart
#' @examples 
#' x = 0:100
#' y = 2005:2016
#' D  <- MortalityForecast.data$dx[paste(x), paste(y)]
#' ex <- MortalityForecast.data$ex
#' exogen <- ex[paste(y)]
#' 
#' MM = c("LC", "FDM", "CoDa", "M6", "M6X")
#' M <- doMortalityModels(data = D, x, y, data.type = "dx",
#'                        models = MM, exogen = exogen)
#' 
#' oex <- getObserved(M, what = "ex")
#' fex <- getFitted(M, what = "ex")
#' rex <- getResiduals(M, what = "ex")
#' @export
doMortalityModels <- function(data, x, y, 
                              data.type = c("qx", "mx", "dx", "lx"),
                              models = c("MRW", "MRWD","LC", "FDM", "PLAT", "CoDa", 
                                         "M4", "M4X", "M5", "M5X", "M6", "M6X"),
                              exogen = NULL, ...) {
  input <- as.list(environment())
  data.type <- match.arg(data.type)
  call <-  match.call()
  
  mx.data <- convertFx(x, data, In = data.type, Out = "mx", lx0 = 1, ...)
  dx.data <- convertFx(x, data, In = data.type, Out = "dx", lx0 = 1, ...)
  # qx.data <- convertFx(x, data, In = data.type, Out = "qx", lx0 = 1) # not needed now
  # lx.data <- convertFx(x, data, In = data.type, Out = "lx", lx0 = 1)

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
  # Mortality Moments Model - PLC (2018)
  if ("M4" %in% models)  M4 <- lenart(data = dx.data, x, y, n = 4)
  if ("M4X" %in% models) M4X <- lenart(data = dx.data, x, y, n = 4, exogen = exogen)
  if ("M5" %in% models)  M5 <- lenart(data = dx.data, x, y, n = 5)
  if ("M5X" %in% models) M5X <- lenart(data = dx.data, x, y, n = 5, exogen = exogen)
  if ("M6" %in% models)  M6 <- lenart(data = dx.data, x, y, n = 6)
  if ("M6X" %in% models) M6X <- lenart(data = dx.data, x, y, n = 6, exogen = exogen)
  if ("M7" %in% models)  M7 <- lenart(data = dx.data, x, y, n = 7)
  if ("M7X" %in% models) M7X <- lenart(data = dx.data, x, y, n = 7, exogen = exogen)

  remove(data, exogen, data.type, models, dx.data, mx.data)
  out <- as.list(environment())
  out <- structure(class = "MortalityModels", out)
  return(out)
}


#' Get Fitted Values
#' 
#' @param object An object of the class \code{MortalityModels}.
#' @param what Specify the type of values to be extracted. 
#' @inheritParams doMortalityModels
#' @export
getFitted <- function(object, 
                      what = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"),
                      ...) {
  what <- match.arg(what)
  Mn   <- object$input$models # Model names
  x    <- object$x
  
  DX <- list()
  for (i in 1:length(Mn)) {
    M <- with(object, get(Mn[i]))
    
    if (Mn[i] %in% c("MRW", "MRWD", "LC", "PLAT")) {
      mx <- exp(fitted(M))
      dx <- convertFx(x, mx, In = "mx", Out = "dx", lx0 = 1, ...)
      
    } else if (Mn[i] %in% c("FDM")) {
      mx <- exp(fitted(M)$y)
      dx <- convertFx(x, mx, In = "mx", Out = "dx", lx0 = 1, ...)
      
    } else {
      dx <- fitted(M)
    }
    DX[[i]] <- dx
  }
  
  fn  <- function(Z) convertFx(x, Z, In = "dx", Out = what, lx0 = 1, ...)
  out <- lapply(DX, fn)
  names(out) <- Mn
  out <- structure(class = "getFitted", out)
  return(out)
}


#' Get Observed Values
#' 
#' @inheritParams getFitted
#' @export
getObserved <- function(object, 
                        what = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"),
                        ...) {
  what <- match.arg(what)
  x    <- object$x
  data <- object$input$data
  In   <- object$input$data.type 
  if (In == what) {
    out <- data
  } else {
    out <- convertFx(x, data, In, Out = what, lx0 = 1, ...)
  }
  return(out)
}






