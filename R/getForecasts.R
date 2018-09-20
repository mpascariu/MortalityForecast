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
      
    } else if (Mn[i] %in% c("HyndmanUllah", "LeeCarter")) {
      mx <- M$rate$mean
      dimnames(mx) <- list(x, y)
      
    } else if (Mn[i] %in% c("MRW", "MRWD")) {
      mx <- exp(M$predicted.values)
      
    } else {
      dx <- M$predicted.values
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
    Z <- convertFx(xx, mx, from = data.in, to = data.out, lx0 = 1)
    # is_zero <- apply(Z, 2, function(x) all(x == 0))
    # Z[, is_zero] <- NA
    Z[paste(object$x), ]
  }
  
  out <- lapply(MX, fn)
  names(out) <- Mn
  out <- structure(class = "getForecasts", out)
  return(out)
}
