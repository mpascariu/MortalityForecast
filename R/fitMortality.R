

#' Wrapper for function dealing with mortality models
#' 
#' @param data Mortality data
#' @param x Vector of ages.
#' @param y Vector of years.
#' @param data.type Type of data in 'data' argument. Options: 
#' \code{"qx", "mx", "dx"}.
#' @param ... Arguments to be passed to or from other methods.
#' @examples 
#' x = 0:110
#' y = 1960:1980
#' dxm <- dxForecast::dxForecast.data$dx$male[paste(x), paste(y)]
#' dx  <- apply(dxm, 2, function(x) x/sum(x))
#' 
#' M <- doMortalityModels(data = dx, x, y, data.type = "dx")
#' P <- doForecasts(M, h = 16, ci = 95, jumpchoice = "actual")
#' 
#' oex <- getObserved(M, type = "ex")
#' fex <- getFitted(M, type = "ex")
#' rex <- getResiduals(M, type = "ex")
#' pex <- getForecasts(P, type = "ex")
#' @export
#' 
doMortalityModels <- function(data, x, y, 
                               data.type = c("qx", "mx", "dx"), ...) {
  input <- as.list(environment())
  data.type <- match.arg(data.type)
  
  mx.data <- fx2gx(x, data, In = data.type, Out = "mx", lx0 = 1)
  dx.data <- fx2gx(x, data, In = data.type, Out = "dx", lx0 = 1)
  qx.data <- fx2gx(x, data, In = data.type, Out = "qx", lx0 = 1)

  # Mortality Moments Model - PLC (2018)
  M1 <- dxForecast::lenart(data = dx.data, x, y, n = 5)
  
  # CoDa-LC (2008)
  M2 <- CoDa::coda(data = dx.data, x, y)
  
  # LC (1992)
  lx0 <- 1e5
  Dx  <- mx.data * lx0
  Ex  <- Dx * 0 + lx0
  wxt <- genWeightMat(ages = x, years = y, clip = 3) # weighting matrix
  M3  <- StMoMo::fit(lc(), Dxt = Dx, Ext = Ex, ages = x, 
             years = y, ages.fit = x, wxt = wxt, verbose = FALSE)
  
  out <- list(call = match.call(), input = input,
              PLC = M1, codaLC = M2, LC = M3, x = x, y = y)
  out <- structure(class = "MortalityModels", out)
  return(out)
}


#' Get Fitted Values
#' 
#' @param object An object of the class \code{MortalityModels}.
#' @param type Specify the type of values to be extracted. 
#' @inheritParams doMortalityModels
#' @export
getFitted <- function(object, 
                      type = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"),
                      ...) {
  type <- match.arg(type)
  x    <- object$x
  
  dx1 <- fitted(object$PLC)
  dx2 <- fitted(object$codaLC)
  mx3 <- exp(fitted(object$LC))
  dx3 <- fx2gx(x, mx3, In = "mx", Out = "dx", lx0 = 1)
  
  dx  <- list(dx1, dx2, dx3)
  fn  <- function(Z) fx2gx(x, Z, In = "dx", Out = type, lx0 = 1)
  out <- lapply(dx, fn)
  names(out) <- c("PLC", "codaLC", "LC")
  return(out)
}

#' Get Observed Values
#' 
#' @inheritParams getFitted
#' @export
getObserved <- function(object, 
                        type = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"),
                        ...) {
  type <- match.arg(type)
  x    <- object$x
  data <- object$input$data
  In   <- object$input$data.type 
  if (In == type) {
    out <- data
  } else {
    out <- fx2gx(x, data, In, Out = type, lx0 = 1)
  }
  return(out)
}


#' Get Deviance Residuals
#' 
#' @inheritParams getFitted
#' @export
getResiduals <- function(object,
                         type = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"),
                         ...) {
  ov  <- getObserved(object, type, ...)
  fv  <- getFitted(object, type, ...)
  fn  <- function(X) ov - X
  out <- lapply(fv, fn)
  return(out)
}





