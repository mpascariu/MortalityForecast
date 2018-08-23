

#' Get Observed Values
#' @inheritParams getFitted
#' @export
getObserved <- function(object, 
                        data.out = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"),
                        ...) {
  x        <- object$x
  data     <- object$input$data
  data.in  <- object$input$data.in 
  data.out <- match.arg(data.out)
  
  if (data.in == data.out) {
    out <- data
  } else {
    out <- convertFx(x, data, from = data.in, to = data.out, lx0 = 1, ...)
  }
  return(out)
}


