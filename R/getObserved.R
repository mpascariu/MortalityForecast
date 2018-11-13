

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



# Fri Aug 31 17:42:02 2018 --------- Marius D. Pascariu ---
#' Extrapolate old-age human mortality curve using mortality laws
#' 
#' @inheritParams MortalityLaws::MortalityLaw
#' @param mx Vector or matrix of age specific death-rates.
#' @param x_fit Ages to be considered in estimating the mortality model parameters.
#' \code{x_fit} can be a subset of \code{x}. However, after the model is identifies
#' fitted values and residuals are computed for all ages in \code{x}.
#' @param x_extr Ages for which to extrapolate the death-rates. 
#' @param ... Other arguments to be passed on to the 
#' \code{\link[MortalityLaws]{MortalityLaw}} function.
#' @seealso 
#' \code{\link[MortalityLaws]{MortalityLaw}}
#' \code{\link[MortalityLaws]{predict.MortalityLaw}}
#' @return An object of class \code{extra_mortality} with the following components:
#'  \item{input}{List with arguments provided in input. Saved for convenience.}
#'  \item{call}{An unevaluated function call, that is, an unevaluated expression 
#'  which consists of the named function applied to the given arguments.}
#'  \item{fitted.model}{An object of class \code{\link[MortalityLaws]{MortalityLaw}}.
#'  Here one can find fitted values, residuals, goodness of fit measures etc.}
#'  \item{values}{A vector or matrix containing the complete mortality data, that is
#'  the modified input data following the extrapolation procedure.}
#' 
#' @examples 
#' # Example 1 - abridged data
#' 
#' # Age-specific death rates
#' mx <- c(.0859, .0034, .0009, .0007, .0016, .0029, .0036, .0054, 
#'         .0053, .0146, .0127, .0269, .0170, .0433, .0371, .0784, 
#'         .0930, .1399, .1875, .2250, .2500, .3000)
#' # Vector of ages
#' x <- c(0, 1, seq(5, 100, by = 5))
#' names(mx) <- x
#' 
#' # Fit the models / Extrapolate the mortality curve
#' x_fit  = c(80, 85, 90, 95, 100)
#' x_extr = 90:110
#' f1 <- extra_mortality(mx, x, x_fit, x_extr, law = "kannisto")
#' 
#' # ----------------------------------------------
#' # Example 2 - 1-year age data
#' 
#' # Age-specific death rates
#' mx1 <- c(.0070, .0082, .0091, .0096, .0108, .0122, .0141, .0150, .0165, .0186, .0205, 
#'          .0229, .0259, .0294, .0334, .0379, .0426, .0482, .0550, .0628, .0716, .0806, 
#'          .0897, .1003, .1149, .1264, .1558, .1563, .1812, .2084, .2298, .2536, .2813, 
#'          .3143, .3352, .3651, .4128)
#' # Vector of ages
#' x1 <- 65:101
#' names(mx1) <- x1
#' 
#' # Fit the models / Extrapolate the mortality curve
#' x_fit = 80:95
#' x_extr = 80:125
#' g1 <- extra_mortality(mx1, x1, x_fit, x_extr, law = "kannisto")
#' 
#' # ----------------------------------------------
#' # Example 3 - Extrapolate mortality for multiple years at once
#' 
#' # Create some data
#' mx_matrix <- matrix(rep(mx1, 3), ncol = 3) %*% diag(c(1, 1.05, 1.1))
#' dimnames(mx_matrix) <- list(age = x1, year = c("year1", "year2", "year3"))
#' 
#' F1 <- extra_mortality(mx_matrix, x = x1, x_fit, x_extr, law = "kannisto")
#' @author Marius D. Pascariu <rpascariu@@outlook.com>
#' @keywords internal
#' @export
extra_mortality <- function(mx, x, x_fit = x, x_extr,
                            law = c("kannisto", "gompertz", "ggompertz", "makeham", 
                                    "beard", "beard_makeham", "quadratic"), 
                            opt.method = c("LF2", "LF1", "LF3", "LF4", "LF5", "LF6",
                                           "poissonL", "binomialL"), ...) {
  # Save the input
  input <- as.list(environment())
  
  # Fit the mortality model
  M <- MortalityLaw(x = x, mx = mx, 
                    fit.this.x = x_fit, 
                    law = match.arg(law), 
                    opt.method = match.arg(opt.method), ...)
  pv <- predict(M, x = x_extr)
  
  # which ages are not to be replaced with fitted values?
  L  <- !(x %in% x_extr)  
  
  # Create the output object
  if (is.vector(mx)) {
    names(mx) <- x
    values <- c(mx[L], pv)
    
  } else {
    rownames(mx) <- x
    values <- rbind(mx[L,], pv)
  }
  
  # Exit
  out <- list(input = input, call = match.call(), fitted.model = M, 
              values = values)
  out <- structure(class = "extra_mortality", out)
  return(out)
}
