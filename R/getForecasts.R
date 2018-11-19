# Thu Aug 23 17:46:12 2018 --------- Marius D. Pascariu ---


#' Get Predicted Values
#' @param object An object of class \code{doForecasts}.
#' @inheritParams evalAccuracy.doBackTesting
#' @seealso \code{\link{doForecasts}}
#' @author Marius D. Pascariu
#' @export
getForecasts <- function(object, 
                         data.out = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"),
                         ...) {
  data.out <- match.arg(data.out)
  Mn <- object$input$object$input$models # Model names
  x  <- object$x
  y  <- object$y  
  
  MX <- list()
  for (i in 1:length(Mn)) {
    M <- with(object, get(Mn[i]))
    
    if (Mn[i] %in% c("LC", "PLAT")) {
      mx <- M$rates
      
    } else if (Mn[i] %in% c("MRW", "MRWD")) {
      mx <- exp(M$predicted.values)
      
    } else if (Mn[i] %in% c("LeeCarter2", "HyndmanUllah")) {
      mx <- M$predicted.values
      
    } else {
      dx <- M$predicted.values
      mx <- convertFx(x, dx, from = "dx", to = "mx", lx0 = 1)
    }
    MX[[i]] <- mx
  }
  
  fn <- function(mx, x_max = 110) {
    xx <- x
    Z <- convertFx(x = xx, data = mx, from = "mx", to = data.out, lx0 = 1)
    Z[paste(object$x), ]
  }
  
  out <- lapply(MX, fn)
  names(out) <- Mn
  out <- structure(class = "getForecasts", out)
  return(out)
}
