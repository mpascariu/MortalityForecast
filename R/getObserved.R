

#' Get Observed Values
#' @inheritParams getFitted
#' @seealso \code{\link{doMortalityModels}}
#' @author Marius D. Pascariu
#' @examples 
#' # For examples go to ?doMortalityModels
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


