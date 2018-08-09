

#' Performe back-testing
#' @param Tdata Test data
#' @inheritParams getAccuracy
#' @inheritParams doMortalityModels
#' @export
doBackTesting <- function(Tdata, object, 
                          data.type = c("qx", "mx", "dx"),
                          type = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"),
                          ...) {
  x   <- object$x
  O   <- fx2gx(x, Tdata, In = data.type, Out = type, lx0 = 1)
  out <- getAccuracy(data = O, object, type)
  return(out)
}




