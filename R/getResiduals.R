

#' Get Deviance Residuals
#' @inheritParams getFitted
#' @seealso \code{\link{doMortalityModels}}
#' @author Marius D. Pascariu
#' @examples 
#' # For examples go to ?doMortalityModels
#' @export
getResiduals <- function(object,
                         data.out = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"),
                         ...) {
  
  data.out <- match.arg(data.out)
  
  ov  <- getObserved(object, data.out, ...)
  fv  <- getFitted(object, data.out, ...)
  fn  <- function(X) ov - X
  out <- lapply(fv, fn)
  out <- structure(class = "getResiduals", out)
  return(out)
}


#' Summary for getResiduals
#' @param object An object of class \code{\link{getResiduals}}.
#' @param digits Number of digits to display.
#' @inheritParams doMortalityModels
#' @keywords internal
#' @export
summary.getResiduals <- function(object, ..., 
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

