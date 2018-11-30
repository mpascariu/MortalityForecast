# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Tue Nov 27 16:15:10 2018
# --------------------------------------------------- #

#' Get Observed Values
#' @inherit get.Fitted description
#' @inheritParams get.Fitted
#' @inheritParams do.MortalityModels
#' @seealso \code{\link{do.MortalityModels}}
#' @author Marius D. Pascariu
#' @examples 
#' # For examples go to ?do.MortalityModels
#' @export
get.Observed <- function(object = NULL, 
                         data.out = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"),
                         data.in = NULL,
                         data = NULL,
                         x = NULL,
                         ...) {
  
  data     <- data %||% object$input$data
  data.in  <- data.in %||% object$input$data.in
  data.out <- match.arg(data.out)
  
  x <- x %||% object$x
  
  if (data.in == data.out) {
    out <- data
    
  } else {
    xx <- x
    out <- convertFx(x = xx, 
                     data = data, 
                     from = data.in, 
                     to = data.out, 
                     lx0 = 1, 
                     ...)
    
    out <- out[paste(x), ]
  }
  
  return(out)
}
