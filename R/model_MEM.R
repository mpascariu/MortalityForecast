# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Sat Dec  1 15:04:54 2018
# --------------------------------------------------- #


#' The Maximum-Entropy Mortality Model 
#' 
#' @param n The maximum order of the moments to be used.
#' @inheritParams do.MortalityModels
#' @return The output is an object of class \code{MEM} with 
#' the components:
#'  \item{input}{List with arguments provided in input. Saved for convenience;}
#'  \item{info}{Short details about the model;}
#'  \item{call}{An unevaluated function call, that is, an unevaluated 
#'  expression which consists of the named function applied to the given 
#'  arguments;}
#'  \item{coefficients}{Estimated coefficients;}
#'  \item{fitted.values}{Fitted values of the estimated model;}
#'  \item{observed.values}{Observed values used in fitting;}
#'  \item{fitted.raw.moments}{Fitted raw moments;}
#'  \item{observed.raw.moments}{Observed raw moments of the input data;}
#'  \item{residuals}{Deviance residuals;} 
#'  \item{random.walk.model}{Object containing the components of the fitted time 
#'  series model to the extrapolate moments;}
#'  \item{x}{Vector of ages used in the fitting;} 
#'  \item{y}{Vector of years used in the fitting.} 
#' @seealso 
#' \code{\link{predict.MEM}} 
#' \code{\link{plot.MEM}} 
#' \code{\link{find.moments}}
#' @details \insertNoCite{pascariu_phd2018}{MortalityForecast}
#' @references \insertAllCited{}
#' @examples 
#' # Data
#' x  <- 0:110
#' y  <- 1985:2014
#' dx <- HMD_male$dx$GBRTENW[paste(x), paste(y)]
#' 
#' # Fit model the Maximum-Entropy Mortality of order 5,
#' # that is using the first 6 statistical moments (moment 0 to 5).
#' 
#' M <- model.MEM(data = dx, x = x, y = y, n = 5)
#' M
#' R <- residuals(M)
#' 
#' plot(M, plotType = "observed")
#' plot(M, plotType = "fitted")
#' 
#' plot(R, plotType = "scatter")
#' plot(R, plotType = "colourmap")
#' plot(R, plotType = "signplot")
#' 
#' # Perform forecasts
#' P <- predict(M, h = 16, x.h = 0:120)
#' 
#' plot(P, plotType = "mean")
#' plot(P, plotType = "lower")
#' plot(P, plotType = "upper")
#' 
#' plot(P, M, plotType = "raw_moments")
#' plot(P, M, plotType = "normalised_moments")

#' @export
model.MEM <- function(data, 
                      x = NULL, 
                      y = NULL, 
                      n = 5, 
                      verbose = FALSE, 
                      ...) {
  
  input <- c(as.list(environment()))
  
  # Info
  modelLN <- "Maximum-Entropy Mortality Model "
  modelSN <- "MEM"
  modelF  <- "d[x,t] = MaxEnt a[n] + M[n, t]"
  info <- list(name = modelLN, name.short = modelSN, formula = modelF)
  
  AY  <- find_ages_and_years(data, x, y)
  x   <- AY$x
  y   <- AY$y
  M   <- find.moments(data, x, y, n)
  orM <- M$raw.moments
  nM  <- M$normalized.moments
  nMT <- log(abs(nM))
  sg  <- sign(nM)
  
  V   <- model.MRW(t(nMT), x = NULL, y, include.drift = TRUE)
  fnM <- fnM <- t(exp(fitted(V)) * as.numeric(sg[1,]))       # fitted normalized moments
  frM <- convert.moments(fnM, from = "normalized", to = "raw") # fitted raw moments
  
  fD  <- find.density(frM[-1, ], x)$density          # fitted dx
  fD  <- cbind(NA, fD) 
  oD  <- apply(data, 2, FUN = function(x) x/sum(x)) # observed dx - same scale as fitted dx
  res <- oD - fD                                    # residuals
  dimnames(oD) = dimnames(fD) = dimnames(res) = list(x = x, y = y)
  
  # Exit
  out <- list(input = input, 
              info = info, 
              call = match.call(),
              coefficients = coef(V),
              fitted.values = fD, 
              observed.values = oD, 
              residuals = res, 
              fitted.raw.moments = frM, 
              observed.raw.moments = orM, 
              random.walk.model = V, 
              x = x, 
              y = y)
  out <- structure(class = 'MEM', out)
  return(out)
}


#' Find age and year vectors
#' @inheritParams model.MEM
#' @keywords internal
#' @export
find_ages_and_years <- function(data, x, y) {
  if (!is.null(x) & nrow(data) != length(x)) {
    stop("'data' and 'x' must have the same dimension!\n", 
         "Check: nrow(data) == length(x)", call. = FALSE)
  } 
  if (!is.null(y) & ncol(data) != length(y)) {
    stop("'data' and 'y' must have the same dimension!\n", 
         "Check: ncol(data) == length(y)", call. = FALSE)
  } 
  
  dn  <- dimnames(data)
  dn1 <- suppressWarnings(as.numeric(dn[[1]]))
  dn2 <- suppressWarnings(as.numeric(dn[[2]]))
  if (is.null(x) & !is.null(dn[1]) & !all(is.na(dn1))) {
    x <- dn1 
  } else {
    if (is.null(x) & is.null(dn[1])) {
      x <- 1:nrow(data) 
    } 
  }
  
  if (is.null(y) & !is.null(dn[2]) & !all(is.na(dn2))) {
    y = dn2 } else {
      if (is.null(y) & is.null(dn[2])) {
        y <- 1:ncol(data) 
      } 
    }
  
  out <- list(x = x, y = y)
  return(out)
}


#' Forecast the age-at-death distribution using the Maximum-Entropy Mortality Model
#' 
#' @param object An object of class \code{MEM};
#' @param x.h Numerical vector indicating the ages to be considered in 
#' reconstruction of the density over the forecast horizon.
#' If \code{NULL}, the number of estimated data points is equal to the number 
#' of fitted data points. This argument can be used for example to estimate a 
#' density between 0 and 130 given the fact that the model was fitted on a 
#' dataset containing values for 0-100 only.
#' @inheritParams do.MortalityForecasts
#' @return The output is a list with the components:
#'  \item{call}{An unevaluated function call, that is, an unevaluated 
#'  expression which consists of the named function applied to the given 
#'  arguments;}
#'  \item{info}{Short details about the model;}
#'  \item{predicted.values}{A list containing the predicted d[x] values given 
#'  by the estimated model over the forecast horizon \code{h};}
#'  \item{conf.intervals}{Confidence intervals for the predicted values;}
#'  \item{predicted.raw.moments}{Predicted raw moments over the forecast horizon 
#'  \code{h};}
#'  \item{random.walk.model}{The estimated multivariate Random Walk model used 
#'  in the extrapolation of the moments;} 
#'  \item{x}{Vector of ages used in prediction;} 
#'  \item{y}{Vector of years used in prediction.}
#' @seealso \code{\link{model.MEM}}
#' @details \insertNoCite{pascariu_phd2018}{MortalityForecast}
#' @references \insertAllCited{}
#' @examples # For examples go to ?model.MEM
#' @export
predict.MEM <- function(object, 
                        h, 
                        x.h = NULL, 
                        level = 95,
                        jumpchoice = c("actual", "fit"),
                        verbose = FALSE, 
                        ...) {
  
  jumpchoice <- match.arg(jumpchoice)
  y.h <- max(object$y) + (1:h)
  x.h <- x.h %||% object$x
  
  # Predict ts model
  W <- predict(object$random.walk.model, h + 1, level)
  L <- W$conf.intervals
  L$mean <- W$predicted.values
  sg  <- object$fitted.raw.moments %>% tail(1) %>% 
    convert.moments(from = "raw", to = "normalized") %>% sign %>% as.numeric
  
  compute_rM <- function(X) {  # Compute raw moments
    t(exp(X) * sg) %>% convert.moments(from = "normalized", to = "raw")
  }
  
  compute_pdf <- function(rM) {
    px <- find.density(rM, x = x.h)$density    # estimate distribution
    
    if (jumpchoice == "actual") {
      ov <- object$observed.values %>% t %>% tail(1) %>% t %>% c
      n <- seq_along(ov)
      J <- rep(1, length(x.h))
      J[n] <- ov / px[n, 1]                          # jumpoff adjustment
      px <- sweep(px, 1, J, FUN = "*")               # adjust
      px <- apply(px, 2, FUN = function(z) z/sum(z)) # close composition
    }
    
    px <- px[, -1]
    dimnames(px) <- list(x.h, y.h)
    return(px)
  }
  
  rM <- lapply(L, compute_rM)
  px <- lapply(rM, compute_pdf)
  
  # Output preparation
  N  <- length(px)
  CI <- list(predicted.values = px[-N], 
             predicted.raw.moments = rM[-N])
  
  out <- list(call = match.call(),
              info = object$info,
              predicted.values = px[[N]], 
              conf.intervals = CI, 
              predicted.raw.moments = rM[[N]], 
              random.walk.model = W, 
              x = x.h, 
              y = y.h)
  out <- structure(class = 'predict.MEM', out)
  return(out)
}


# S3 ----------------------------------------------
#' @rdname residuals.Oeppen
#' @export
residuals.MEM <- function(object, ...) {
  residuals_default(object, ...)
}


#' @rdname print_default
#' @export
print.MEM <- function(x, ...) {
  print_default(x, ...)
  cat("Moments in fit: ", paste0("M", c(0, x$input$n), collapse = " - "))
  cat("\n")
}


#' @rdname print_default
#' @export
print.predict.MEM <- function(x, ...) {
  print_predict_default(x, ...)
  cat("Moments in forecast: ", paste0("M", c(0, ncol(x$predicted.raw.moments) - 1), 
                                      collapse = " - "))
  cat("\n")
}


