

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
  
  DX <- list()
  for (i in 1:length(Mn)) {
    M <- with(object, get(Mn[i]))
    
    if (Mn[i] %in% c("MRW", "MRWD", "LC", "PLAT")) {
      mx <- exp(fitted(M))
      dx <- convertFx(x, mx, from = "mx", to = "dx", lx0 = 1, ...)
      
    } else if (Mn[i] %in% c("HyndmanUllah", "LeeCarter")) {
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


#' Summary for getFitted
#' @param object An object of class \code{\link{getFitted}}.
#' @inheritParams summary.getResiduals
#' @keywords internal
#' @export
summary.getFitted <- function(object, ..., digits = max(4L, getOption("digits") - 2L)) {
  summary.getResiduals(object, ..., digits)
}
