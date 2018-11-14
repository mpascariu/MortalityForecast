

#' Get Observed Values
#' @inheritParams getFitted
#' @inheritParams doMortalityModels
#' @seealso \code{\link{doMortalityModels}}
#' @author Marius D. Pascariu
#' @examples 
#' # For examples go to ?doMortalityModels
#' @export
getObserved <- function(object = NULL, 
                        data.out = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"),
                        data.in = NULL,
                        data = NULL,
                        x = NULL,
                        ...) {
  
  x <- x %||% object$x
  data <- data %||% object$input$data
  data.in  <- data.in %||% object$input$data.in
  data.out <- match.arg(data.out)
  
  if (data.in == data.out) {
    out <- data
    
  } else {
    xx <- x
    out <- convertFx(x = xx, data = data, 
                     from = data.in, to = data.out, 
                     lx0 = 1, ...)
    out <- out[paste(x), ]
  }
  
  return(out)
}

