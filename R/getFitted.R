

#' Get Fitted Values
#' @param object An object of class \code{MortalityModels}.
#' @inheritParams evalAccuracy.doBackTesting
#' @seealso \code{\link{doMortalityModels}}
#' @author Marius D. Pascariu
#' @examples 
#' # For examples go to ?doMortalityModels
#' @export
getFitted <- function(object, 
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
      
    } else if (Mn[i] %in% c("HyndmanUllah")) {
      mx <- exp(fitted(M)$y)
      
    } else if (Mn[i] %in% c("LeeCarter")) {
      mx <- fitted(M)$y
      
    } else if (Mn[i] %in% c("LeeCarter2")) {
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
  out <- structure(class = "getFitted", out)
  return(out)
}


#' Summary for getFitted
#' @param object An object of class \code{\link{getFitted}}.
#' @inheritParams summary.getResiduals
#' @keywords internal
#' @export
summary.getFitted <- function(object, ..., digits = NULL) {
  digits <- digits %||% max(4L, getOption("digits") - 2L)
  summary.getResiduals(object, ..., digits)
}
