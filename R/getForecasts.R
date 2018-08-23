# Thu Aug 23 17:46:12 2018 --------- Marius D. Pascariu ---


#' Get Predicted Values
#' @param object An object of class \code{doForecasts}.
#' @inheritParams getAccuracy
#' @export
getForecasts <- function(object, 
                         data.out = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"),
                         ...) {
  data.out <- match.arg(data.out)
  Mn <- object$input$object$input$models # Model names
  x  <- object$x
  y  <- object$y  
  
  DX <- list()
  for (i in 1:length(Mn)) {
    M <- with(object, get(Mn[i]))
    
    if (Mn[i] %in% c("LC", "PLAT")) {
      mx <- M$rates
      dx <- convertFx(x, mx, from = "mx", to = "dx", lx0 = 1, ...)
      
    } else if (Mn[i] %in% c("HyndmanUllah", "LeeCarter")) {
      mx <- M$rate$mean
      dimnames(mx) <- list(x, y)
      dx <- convertFx(x, mx, from = "mx", to = "dx", lx0 = 1, ...)
      
    } else if (Mn[i] %in% c("MRW", "MRWD")) {
      mx <- exp(M$predicted.values)
      dx <- convertFx(x, mx, from = "mx", to = "dx", lx0 = 1, ...)
      
    } else {
      dx <- M$predicted.values
    }
    DX[[i]] <- dx
  }
  
  fn  <- function(Z) convertFx(x, Z, from = "dx", to = data.out, lx0 = 1)
  out <- lapply(DX, fn)
  names(out) <- Mn
  out <- structure(class = "getForecasts", out)
  return(out)
}
