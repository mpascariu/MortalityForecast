

#' Wrapper for function dealing with mortality models
#' 
#' @param data Mortality data
#' @param x Vector of ages.
#' @param y Vector of years.
#' @param data.type Type of data in 'data' argument. Options: 
#' \code{"qx", "mx", "dx"}.
#' @param ... Arguments to be passed to or from other methods.
#' @inheritParams dxForecast::lenart
#' @examples 
#' x = 0:110
#' y = 1980:1999
#' h = 17
#' dxm <- dxForecast::dxForecast.data$dx$male[paste(x), paste(y)]
#' ex <- dxForecast::dxForecast.data$ex$male
#' exogen <- ex[paste(y)]
#' 
#' M <- doMortalityModels(data = dxm, x, y, data.type = "dx", exogen = exogen)
#' P <- doForecasts(M, h, ci = 95, jumpchoice = "actual")
#' 
#' 
#' oex <- getObserved(M, type = "ex")
#' fex <- getFitted(M, type = "ex")
#' rex <- getResiduals(M, type = "ex")
#' pex <- getForecasts(P, type = "ex")
#' 
#' 
#' y2 <- max(y) + 1:h
#' Tdata <- dxForecast::dxForecast.data$dx$male[paste(x), paste(y2)]
#' doBackTesting(Tdata, P, data.type = "dx", type = "ex")
#' @export
#' 
doMortalityModels <- function(data, x, y, 
                              data.type = c("qx", "mx", "dx"), exogen = NULL, ...) {
  input <- as.list(environment())
  data.type <- match.arg(data.type)
  
  mx.data <- fx2gx(x, data, In = data.type, Out = "mx", lx0 = 1)
  dx.data <- fx2gx(x, data, In = data.type, Out = "dx", lx0 = 1)
  qx.data <- fx2gx(x, data, In = data.type, Out = "qx", lx0 = 1)

  # # The Naive model
  # M0 <- list(fitted.values = dx.data * 0 + dx.data[, length(y)])

  # Random Walk with drift
  # M1 <- StMoMo::mrwd(mx.data)
  
  # LC (1992)
  lx0 <- 1e5
  Dx  <- mx.data * lx0
  Ex  <- Dx * 0 + lx0
  wxt <- genWeightMat(ages = x, years = y, clip = 3) # weighting matrix
  M2  <- StMoMo::fit(lc(), Dxt = Dx, Ext = Ex, ages = x, 
             years = y, ages.fit = x, wxt = wxt, verbose = FALSE)
  
  # CoDa-LC (2008)
  M3 <- CoDa::coda(data = dx.data, x, y)
  
  # Mortality Moments Model - PLC (2018)
  M4 <- dxForecast::lenart(data = dx.data, x, y, n = 5)
  
  # Mortality Moments Model - PLC (2018)
  M5 <- dxForecast::lenart(data = dx.data, x, y, n = 5, exogen = exogen)
  

  # Mn <- c("RWD", "Lee-Carter", "CoDa-LC", "PLC")  
  # Mn <- c("Lee-Carter", "PLC", "CoDa-LC")  
  Mn <- c("M2", "M3", "M4", "M5")  
  out <- list(call = match.call(), input = input, x = x, y = y,
              # M1 = M1, M2 = M2, M3 = M3, M4 = M4, model.names = Mn)
              M2 = M2, M3 = M3, M4 = M4, M5 = M5, model.names = Mn)
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
  
  # mx1 <- fitted(object$M1)
  # dx1 <- fx2gx(x, mx1, In = "mx", Out = "dx", lx0 = 1)
  mx2 <- exp(fitted(object$M2))
  dx2 <- fx2gx(x, mx2, In = "mx", Out = "dx", lx0 = 1)
  dx3 <- fitted(object$M3)
  dx4 <- fitted(object$M4)
  dx5 <- fitted(object$M5)
  
  # dx  <- list(dx1, dx2, dx3, dx4)
  dx  <- list(dx2, dx3, dx4, dx5)
  fn  <- function(Z) fx2gx(x, Z, In = "dx", Out = type, lx0 = 1)
  out <- lapply(dx, fn)
  names(out) <- object$model.names
  out <- structure(class = "getFitted", out)
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






