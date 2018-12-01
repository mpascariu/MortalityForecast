# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Thu Nov 29 14:25:29 2018
# --------------------------------------------------- #


#' Get Fitted Values
#' 
#' Extract the values from an object containing an estimated mortality model
#' (e.g. \code{LeeCarter, Oeppen} etc.). This function allows the user to read 
#' the values already converted to a life table index (e.g. life expectancy) 
#' by specifying the requiered format in the argument \code{data.out}.
#' @param object An object of class \code{MortalityModels}.
#' @inheritParams evalAccuracy.BackTesting
#' @seealso \code{\link{do.MortalityModels}}
#' @author Marius D. Pascariu
#' @examples # For examples go to ?do.MortalityModels
#' @export
get.Fitted <- function(object, 
                       data.out = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"),
                       ...) {
  
  data.out <- match.arg(data.out)
  Mn   <- object$input$models # Model names
  x    <- object$x
  
  MX <- list()
  for (i in seq_along(Mn)) {
    M <- with(object, get(Mn[i]))
    
    if (Mn[i] %in% c("MRW", "MRWD", "LC", "PLAT")) {
      mx <- exp(fitted(M))
      
    } else if (Mn[i] %in% c("LeeCarter", "LiLee", "HyndmanUllah")) {
      mx <- fitted(M)
      
    } else {
      dx <- fitted(M)
      mx <- convertFx(x = x, data = dx, from = "dx", to = "mx", lx0 = 1)
    }
    MX[[i]] <- mx
  }
  
  fn <- function(mx, x_max = 110) {
    xx <- x
    Z <- convertFx(x = xx, data = mx, from = "mx", to = data.out, lx0 = 1, ...)
    Z[paste(object$x), ]
  }
  
  out <- lapply(MX, fn)
  names(out) <- Mn
  out <- structure(class = "get.Fitted", out)
  return(out)
}


#' @rdname summary.get.Residuals
summary.get.Fitted <- function(object, ..., digits = NULL) {
  digits <- digits %||% max(4L, getOption("digits") - 2L)
  summary.get.Residuals(object, ..., digits)
}
