# --------------------------------------------------- #
# Author: Marius D. PASCARIU
# License: GNU General Public License v3.0
# Last update: Mon Jan 21 11:16:48 2019 
# --------------------------------------------------- #


#' The Renshaw-Haberman Mortality Model imported from StMoMo package
#' 
#' The Renshaw-Haberman mortality model is a Lee-Carter model with cohort effects.
#' @inheritParams do.MortalityModels
#' @inheritParams StMoMo::rh
#' @inheritParams StMoMo::fit.rh
#' @param radix Radix.
#' @details \insertNoCite{renshaw2006}{MortalityForecast}
#' @seealso 
#' \code{\link{predict.RenshawHaberman}}
#' @references \insertAllCited{}
#' @examples 
#' # Data
#' x  <- 0:89
#' y  <- 2008:2014
#' mx <- HMD_male$mx$GBRTENW[paste(x), paste(y)]
#' 
#' # Fit the model
#' M <- model.RenshawHaberman(data = mx, x = x, y = y)
#' M
#' summary(M)
#' 
#' # Check residuals
#' R <- residuals(M)
#' 
#' plot(R, plotType = "scatter")
#' plot(R, plotType = "colourmap")
#' plot(R, plotType = "signplot")
#' 
#' # Forecast
#' P <- predict(M, h = 5)
#' @export
model.RenshawHaberman <- function(data, 
                                  x, 
                                  y, 
                                  link = c("log", "logit"), 
                                  cohortAgeFun = c("1", "NP"),
                                  approxConst = FALSE,
                                  radix = 1e5, 
                                  verbose = FALSE) {
  link <- match.arg(link)
  cohortAgeFun <- match.arg(cohortAgeFun)
  
  input <- c(as.list(environment()))
  # Info
  modelLN <- "Renshaw-Haberman Mortality Model"                # long name
  modelSN <- "RH"                                              # short name
  modelF  <- paste(link, "m[x,t] = a[x] + b[x] k[t] + g[t-x]") # formula
  info <- list(name = modelLN, name.short = modelSN, formula = modelF)
  
  M <- StMoMo::fit(object = rh(link = link, 
                               cohortAgeFun = cohortAgeFun, 
                               approxConst = approxConst), 
                   Dxt = data * radix, 
                   Ext = data * 0 + radix, 
                   ages = x, 
                   years = y, 
                   verbose = verbose)
  
  cf <- lapply(coef(M), as.numeric)
  fv <- fitted(M)  # fitted values
  fv <- switch(link,
               log = exp(fv),
               logit = 1 / (1 + exp(-fv))
               )
  dimnames(fv) <- list(x, y)
  
  resid <- data - fv # residuals
  
  # Exit
  out <- list(input = input, 
              info = info, 
              call = match.call(), 
              coefficients = cf, 
              fitted.values = fv, 
              observed.values = data,
              residuals = resid, 
              x = x, 
              y = y,
              StMoMo.object = M)
  out <- structure(class = 'RenshawHaberman', out)
  return(out)
}


#' Forecast age-specific death rates using the Renshaw-Haberman model
#' 
#' @param object An object of class \code{RenshawHaberman}.
#' @param order A specification of the ARIMA model for the cohort effect: 
#' the three components (p, d, q) are the AR order, the degree of differencing, 
#' and the MA. Default: \code{c(1, 1, 0)}.
#' @param include.drift a logical value indicating if the ARIMA model should 
#' include a constant value. Default: \code{TRUE}.
#' @param ... Other arguments to be passed on to 
#' \code{\link[StMoMo]{forecast.fitStMoMo}}
#' @inheritParams predict.Oeppen
#' @details \insertNoCite{renshaw2006}{MortalityForecast}
#' @seealso 
#' \code{\link{model.RenshawHaberman}}
#' @references \insertAllCited{}
#' @examples # For examples go to ?model.RenshawHaberman
#' @export
predict.RenshawHaberman <- function(object, 
                                    h, 
                                    order = c(2, 1, 0), 
                                    include.drift = TRUE,
                                    level = c(80, 95),
                                    jumpchoice = c("actual", "fit"),
                                    verbose = TRUE, 
                                    ...) {
  
  # Timeline
  bop <- max(object$y) + 1
  eop <- bop + h - 1
  fcy <- bop:eop
  
  jumpchoice <- match.arg(jumpchoice)
  
  P <- forecast(object = object$StMoMo.object, 
                h = h, 
                gc.order = order,
                gc.include.constant = include.drift,
                jumpchoice = jumpchoice, 
                level = level, 
                ...)
  
  lw <- paste0('L', level)
  up <- paste0('U', level)
  Cnames <- c('mean', lw, up)
  G <- K <- matrix(NA, nrow = h, ncol = 1 + 2 * length(level), 
                   dimnames = list(fcy, Cnames))
  
  K[, "mean"] <- P$kt.f$mean
  K[, lw] <- P$kt.f$lower[1,,]
  K[, up] <- P$kt.f$upper[1,,]
  
  G[, "mean"] <- P$gc.f$mean
  G[, lw] <- P$gc.f$lower
  G[, up] <- P$gc.f$upper
  
  # Exit
  out <- list(call = match.call(), 
              info = object$info,
              kt = K,
              kt.arima = NULL,
              gc = G,
              gc.arima = NULL,
              predicted.values = P$rates, 
              conf.intervals = NULL, # not implemented yet.
              x = object$x, 
              y = fcy,
              StMoMo.object = P)
  out <- structure(class = 'predict.RenshawHaberman', out)
  return(out)
}


# S3 ----------------------------------------------
#' @rdname residuals.Oeppen
#' @export
residuals.RenshawHaberman <- function(object, ...){
  residuals_default(object, ...)
}


#' @rdname print_default
#' @export
print.RenshawHaberman <- function(x, ...) {
  print_default(x, ...)
}


#' @rdname summary.Oeppen
#' @export
summary.RenshawHaberman <- function(object, ...) {
  axbx <- data.frame(x = object$x,
                     ax = object$coefficients$ax, 
                     bx = object$coefficients$bx,
                     b0x = object$coefficients$b0x)
  
  out <- structure(class = 'summary.RenshawHaberman', 
                   list(A = axbx, 
                        kt = object$coefficients$kt, 
                        gc = object$coefficients$gc,
                        call = object$call, 
                        info = object$info,
                        y = object$y, 
                        x = object$x,
                        cohorts = object$StMoMo.object$cohorts))
  return(out)
}


#' @rdname print_default
#' @export
print.summary.RenshawHaberman <- function(x, ...){
  cat('\nFit  :', x$info$name)
  cat('\nModel:', x$info$formula)
  cat('\n\nCoefficients:\n')
  A <- head_tail(x$A, digits = 4, hlength = 4, tlength = 4)
  K <- head_tail(data.frame('|', as.integer(x$y), x$kt),
                 digits = 4, hlength = 4, tlength = 4)
  C <- head_tail(data.frame('|', as.integer(x$cohorts), x$gc),
                 digits = 4, hlength = 4, tlength = 4)
  D <- data.frame(A, K, C)
  colnames(D) <- c("x (age)", "a[x]", "b[x]", "b0[x]", 
                   "|", "t (year)", "k[t]", 
                   ".|", "c (cohort)", "g[c]")
  print(D, row.names = FALSE)
  cat('\n')
}


#' @rdname print_default
#' @export
print.predict.RenshawHaberman <- function(x, ...) {
  print_predict_default(x, ...)
  cat('k[t] method: Random-Walk with drift')
  cat('\ng[c]-ARIMA method:', 'ARIMA(0, 1, 1) w drift')
  cat('\n')
}


# ---------------------------------------------------------
#' Hack gnm::Mult weird dependency
#' The way this function is called in StMoMo functions makes it necesary 
#' to have the gnm dependecy. A workaround is redefine it in the package and 
#' export it from here, and import the gnm package. 
#' @inheritParams gnm::Mult
#' @keywords internal
#' @source \code{\link[gnm]{Mult}} or from here: 
#' \url{https://github.com/hturner/gnm/blob/master/R/Mult.R}
#' @keywords internal
#' @export
Mult <- function(..., inst = NULL){
  if ("multiplicity" %in% names(match.call()[-1]))
    stop("multiplicity argument of Mult has been replaced by",
         "\"inst\" argument.")
  dots <- match.call(expand.dots = FALSE)[["..."]]
  list(predictors = dots,
       term = function(predLabels, ...) {
         paste("(", paste(predLabels, collapse = ")*("), ")", sep = "")
       },
       call = as.expression(match.call()),
       match = seq(dots))
}
class(Mult) <- "nonlin"
