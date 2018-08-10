

#' Wrapper for function dealing with mortality models
#' 
#' @param data Mortality data
#' @param x Vector of ages.
#' @param y Vector of years.
#' @param data.type Type of data in 'data' argument. Options: 
#' \code{"qx", "mx", "dx", "lx"}.
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
                              data.type = c("qx", "mx", "dx", "lx"), exogen = NULL, ...) {
  input <- as.list(environment())
  data.type <- match.arg(data.type)
  call <-  match.call()
  
  mx.data <- convertFx(x, data, In = data.type, Out = "mx", lx0 = 1, ...)
  dx.data <- convertFx(x, data, In = data.type, Out = "dx", lx0 = 1, ...)
  # qx.data <- convertFx(x, data, In = data.type, Out = "qx", lx0 = 1) # not needed now
  # lx.data <- convertFx(x, data, In = data.type, Out = "lx", lx0 = 1)

  # # The Naive model
  # M0 <- list(fitted.values = dx.data * 0 + dx.data[, length(y)])

  # Random Walk with drift
  # M1 <- StMoMo::mrwd(mx.data)
  
  # LC (1992)
  lcm <- function() {
    lx0 <- 1e5
    Dx  <- mx.data * lx0
    Ex  <- Dx * 0 + lx0
    wxt <- StMoMo::genWeightMat(ages = x, years = y, clip = 3) # weighting matrix
    
    StMoMo::fit(StMoMo::lc(), Dxt = Dx, Ext = Ex, ages = x, years = y, 
                ages.fit = x, wxt = wxt, verbose = FALSE)
  }
  M2 <- lcm()
  # CoDa-LC (2008)
  M3 <- CoDa::coda(data = dx.data, x, y)
  # Mortality Moments Model - PLC (2018)
  M4 <- dxForecast::lenart(data = dx.data, x, y, n = 4)
  # Mortality Moments Model - PLC (2018)
  M5 <- dxForecast::lenart(data = dx.data, x, y, n = 4, exogen = exogen)
  # Mortality Moments Model - PLC (2018)
  M6 <- dxForecast::lenart(data = dx.data, x, y, n = 5)
  # Mortality Moments Model - PLC (2018)
  M7 <- dxForecast::lenart(data = dx.data, x, y, n = 5, exogen = exogen)
  # Mortality Moments Model - PLC (2018)
  M8 <- dxForecast::lenart(data = dx.data, x, y, n = 6)
  # Mortality Moments Model - PLC (2018)
  M9 <- dxForecast::lenart(data = dx.data, x, y, n = 6, exogen = exogen)
  
  model.names <- c("LC", "CoDa", "M4", "M4X", "M5", "M5X", "M6", "M6X")
  remove(data, exogen, data.type, dx.data, mx.data, lcm)
  out <- as.list(environment())
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
  Mn   <- object$model.names
  x    <- object$x
  
  # mx1 <- fitted(object$M1)
  # dx1 <- convertFx(x, mx1, In = "mx", Out = "dx", lx0 = 1)
  mx2 <- exp(fitted(object$M2))
  dx2 <- convertFx(x, mx2, In = "mx", Out = "dx", lx0 = 1)
  
  dx3 = dx4 = dx5 = dx6 = dx7 = dx8 = dx9 <- NULL
  for (i in 3:9) {
    Mi <- with(object, get(paste0("M", i)))
    assign(paste0("dx", i), fitted(Mi))
  }
  
  dx  <- list(dx2, dx3, dx4, dx5, dx6, dx7, dx8, dx9)
  fn  <- function(Z) convertFx(x, Z, In = "dx", Out = type, lx0 = 1)
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
    out <- convertFx(x, data, In, Out = type, lx0 = 1)
  }
  return(out)
}






