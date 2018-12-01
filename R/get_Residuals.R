# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Thu Nov 29 14:13:30 2018
# --------------------------------------------------- #


#' Get Deviance Residuals
#' @inherit get.Fitted description
#' @inheritParams get.Fitted
#' @seealso \code{\link{do.MortalityModels}}
#' @author Marius D. Pascariu
#' @examples 
#' # For examples go to ?do.MortalityModels
#' @export
get.Residuals <- function(object,
                          data.out = c("qx", "mx", "dx", "lx", 
                                       "Lx", "Tx", "ex"),
                          ...) {
  
  data.out <- match.arg(data.out)
  
  ov  <- get.Observed(object, data.out, ...)
  fv  <- get.Fitted(object, data.out, ...)
  fn  <- function(X) ov - X
  out <- lapply(fv, fn)
  out <- structure(class = "get.Residuals", out)
  
  return(out)
}


#' Generic Summary
#' @param object An object of class \code{\link{get.Residuals}}.
#' @param digits Number of digits to display.
#' @inheritParams do.MortalityModels
#' @keywords internal
#' @export
summary.get.Residuals <- function(object, ..., 
                                  digits = max(4L, getOption("digits") - 2L)) {
  rn <- names(object)
  cn <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
  N <- length(rn)
  
  fn <- function(Z) {
    z <- unlist(Z)
    z <- z[!is.na(z)]
    summary(z)
  } 
  
  O <- lapply(object, fn)
  A <- matrix(unlist(O), nrow = N, byrow = T)
  colnames(A) <- cn
  out <- data.frame(model = rn, round(A, digits))
  
  out <- as.tibble(out)
  return(out)
}

