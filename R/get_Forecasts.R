# --------------------------------------------------- #
# Author: Marius D. PASCARIU
# License: GNU General Public License v3.0
# Last update: Mon Jan 21 17:22:23 2019 
# --------------------------------------------------- #
# 
#' Get Predicted Values
#' 
#' @inherit get.Fitted description
#' @param object An object of class \code{MortalityForecasts}.
#' @inheritParams evalAccuracy.BackTesting
#' @seealso \code{\link{do.MortalityForecasts}}
#' @author Marius D. Pascariu
#' @export
get.Forecasts <- function(object, 
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
      
    } else if (Mn[i] %in% c("LeeCarter", "HyndmanUllah", "LiLee", 
                            "RenshawHaberman")) {
      mx <- M$predicted.values
      
    } else {
      if ("x.h" %in% names(list(...))) x <- list(...)$x.h # this is for MEM models only
      
      dx <- M$predicted.values
      mx <- convertFx(x = x, data = dx, from = "dx", to = "mx", lx0 = 1)
    }
    MX[[i]] <- mx
  }
  
  fn <- function(mx) {
    xx <- x
    Z <- convertFx(x = xx, data = mx, from = "mx", to = data.out, lx0 = 1)
    Z[paste(object$x), ]
  }
  
  out <- lapply(MX, fn)
  names(out) <- Mn
  out <- structure(class = "get.Forecasts", out)
  return(out)
}
