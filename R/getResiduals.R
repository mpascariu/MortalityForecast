

#' Get Deviance Residuals
#' 
#' @inheritParams getFitted
#' @export
getResiduals <- function(object,
                         what = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"),
                         ...) {
  ov  <- getObserved(object, what, ...)
  fv  <- getFitted(object, what, ...)
  fn  <- function(X) ov - X
  out <- lapply(fv, fn)
  out <- structure(class = "getResiduals", out)
  return(out)
}


#' Summary for getResiduals
#' 
#' @param object An object of the class \code{\link{getResiduals}}.
#' @param digits Number of digits to display.
#' @inheritParams doMortalityModels
#' @export
summary.getResiduals <- function(object, ..., 
                                 digits = max(3L, getOption("digits") - 3L)) {
  a <- lapply(object, FUN = function(k) summary(unlist(k)))
  A <- matrix(unlist(a), nrow = 4, byrow = T)
  # models <- c("RWD", "Lee-Carter", "CoDa-LC", "PLC")
  models <- c("Lee-Carter", "CoDa-LC", "PLC")
  st <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
  dimnames(A) <- list(Models = models, st)
  out <- list(resid = A, digits = digits)
  class(out) <- "summary.getResiduals"
  return(out)
}


#' Print summary.getResiduals
#' @param x An object of the class \code{summary.getResiduals}
#' @inheritParams doMortalityModels
#' @keywords internal
#' @export
print.summary.getResiduals <- function(x, ...) {
  cat('\nDeviance Residuals:\n')
  print(round(x$resid, x$digits))
}
