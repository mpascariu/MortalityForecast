

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
      
    } else if (Mn[i] %in% c("HyndmanUllah", "LeeCarter")) {
      mx <- exp(fitted(M)$y)
      
    } else {
      dx <- fitted(M)
      mx <- convertFx(x, dx, from = "dx", to = "mx", lx0 = 1)
    }
    MX[[i]] <- mx
  }
  
  fn <- function(mx, x_max = 110) {
    xx <- x
    data.in <- "mx"
    # if (max(x) < x_max) {
    #   x_fit <- 80:max(x - 2)
    #   x_extr<- max(x - 2):x_max
    #   mx    <- extra_mortality(mx, x, x_fit, x_extr, law = "kannisto")$values
    #   xx    <- min(x):x_max
    # }
    Z <- convertFx(xx, mx, from = data.in, to = data.out, lx0 = 1, ...)
    # is_zero <- apply(Z, 2, function(x) all(x == 0))
    # Z[, is_zero] <- NA
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
summary.getFitted <- function(object, ..., digits = max(4L, getOption("digits") - 2L)) {
  summary.getResiduals(object, ..., digits)
}
