
#' Predict distribution of deaths using CoDa model.
#' 
#' @param object coda object
#' @param h Number of years to be forecast in the future
#' @param order A specification of the non-seasonal part of the ARIMA model: 
#'  the three components (p, d, q) are the AR order, the degree of differencing, 
#'  and the MA order. If \code{order = NULL}, the ARIMA order will be estimated 
#'  automatically using the KPPS algorithm.
#' @param include.drift Logical. Should the ARIMA model include a linear drift term?
#'  If \code{include.drift = NULL}, the model will be estimated automatically.
#' @param method Fitting method: maximum likelihood or minimize conditional 
#'  sum-of-squares. Options to use:
#'  conditional-sum-of-squares (\code{"CSS-ML"}), maximum likelihood (\code{"ML"}) 
#'  and \code{"CSS"}.
#' @param level Confidence level for prediction intervals.
#' @param ... Additional arguments to be passed to \code{\link[forecast]{Arima}}
#' @param jumpchoice Method used for computation of jumpchoice. 
#'  Possibilities: \code{"actual"} (use actual rates from final year) 
#'  and \code{"fit"} (use fitted rates).
#' @return The output is an object of class \code{"predict.coda"} with the components:
#' @return \item{call}{An unevaluated function call, that is, an unevaluated 
#' expression which consists of the named function applied to the given arguments.}
#' @return \item{predicted.values}{A list containing the predicted values together
#' with the associated prediction intervals given by the estimated \code{link{coda}} 
#' model over the forecast horizon \code{h}.}
#' @return \item{kt}{The extrapolated kt parameters.}
#' @return \item{conf.intervals}{The extrapolated kt parameters.}
#' @return \item{deep}{An object of class \code{ARIMA} that contains all the
#' components of the fitted time series model used in \code{kt} prediction.} 
#' @return \item{x}{Vector of ages used in prediction.} 
#' @return \item{y}{Vector of years used in prediction.} 
#' @examples 
#' # Example 1 ----------------------
#' # Fit CoDa Mortality Model
#' D <- MortalityForecast.data$dx
#' M <- coda(D)
#' 
#' # Predict life expectancy 20 years in the future using CoDa model
#' P <- predict(M, h = 20)
#' 
#' # Example 2 ----------------------
#' # One can specify manually the ARIMA order, a drift to be included or not 
#' # and the jump choice of the first forecast year.
#' P2 <- predict(M, h = 20, order = c(0,1,0), include.drift = TRUE, jumpchoice = "fit")
#' 
#' \dontrun{
#' # Example 3 ----------------------
#' # Compute life tables using forecast values using the MortalityLaws R package
#' library(MortalityLaws)
#' dx <- P$predicted.values$mean
#' lt <- LifeTable(x = P$x, dx = dx)
#' }
#' @export
#' 
predict.coda <- function(object, h, order = NULL, include.drift = NULL,
                         method = "ML", level = c(80, 95), 
                         jumpchoice = c("actual", "fit"), ...){
  dx  <- t(object$input$data)
  bop <- max(object$y) + 1
  eop <- bop + h - 1
  fcy <- bop:eop
  jc  <- jumpchoice[1]
  cf  <- coef(object)
  kt  <- cf$kt
  
  # forecast kt; ax and bx are time independent.
  ts_auto = auto.arima(kt)
  
  AO  <- order %||% arimaorder(ts_auto)
  ID  <- include.drift %||% any(names(coef(ts_auto)) %in% "drift")
  tsm <- forecast::Arima(y = kt, order = AO, include.drift = ID, method = method, ...)
  tsf <- forecast(tsm, h = h, level = level)  # time series forecast
  
  fkt <- data.frame(tsf$mean, tsf$lower, tsf$upper) # forecast kt
  Cnames <- c('mean', paste0('L', level), paste0('U', level))
  colnames(fkt) <- Cnames
  
  fdx <- compute_dx(dx = dx, kt = fkt, ax = cf$ax, bx = cf$bx, # forecast dx
                    fit = t(fitted(object)), y = fcy, jumpchoice = jc)
  pv <- fdx[[1]]
  CI <- fdx[-1]
  names(CI) <- Cnames[-1]
  
  out <- list(call = match.call(), predicted.values = pv, 
              kt = fkt, conf.intervals = CI,
              deep = tsm, x = object$x, y = fcy)
  out <- structure(class = 'predict.coda', out)
  return(out)
}


#' Internal function
#' 
#' @inheritParams coda
#' @inheritParams predict.coda
#' @param kt Estimated kt vector of parameters
#' @param ax Estimated ax vector of parameters
#' @param bx Estimated bx vector of parameters
#' @param fit Fitted values
#' @keywords internal
#' 
compute_dx <- function(dx, kt, ax, bx, fit, y, jumpchoice) {
  
  if (is.data.frame(kt)) {
    pred <- list()
    for (i in 1:ncol(kt)) {
      pred[[i]] <- compute_dx(dx, kt = kt[, i], ax, bx, fit, y, jumpchoice)
      colnames(pred[[i]]) <- y
    }
    return(pred)
    
  } else {
    dx_nrow  <- nrow(dx)
    close.dx <- acomp(dx)
    jump_off <- as.numeric(close.dx[dx_nrow, ]/fit[dx_nrow, ])
    clr_proj <- matrix(kt, ncol = 1) %*% bx
    bk_      <- unclass(clrInv(clr_proj))
    dx_      <- sweep(bk_, 2, ax, FUN = "*")
    if (jumpchoice == 'actual') dx_ <- sweep(dx_, 2, jump_off, FUN = "*")
    out <- t(dx_/rowSums(dx_))
    rownames(out) <- colnames(dx)
    return(out)
  }
} 


# S3 ----------------------------------------------


#' Print predict.coda
#' @param x An object of class \code{"predict.coda"}
#' @inheritParams print.coda
#' @keywords internal
#' @export
#' 
print.predict.coda <- function(x, ...) {
  cat('\nForecast: Compositional-Data Lee-Carter Mortality Model')
  cat('\nModel   : clr d[x,t] = a[x] + b[x]k[t]')
  cat('\nCall    : '); print(x$call)
  cat('\nkt TS method     :', arima.string1(x$deep, padding = TRUE))
  cat('\nAges  in forecast:', paste(range(x$x), collapse = ' - '))
  cat('\nYears in forecast:', paste(range(x$y), collapse = ' - '))
  cat('\n')
}


#' Identify ARIMA model - internal function
#' @param object An object generate by Arima function
#' @param padding Logical.
#' @keywords internal
#' 
arima.string1 <- function(object, padding = FALSE) {
  order  <- object$arma[c(1, 6, 2, 3, 7, 4, 5)]
  nc     <- names(coef(object))
  result <- paste0("ARIMA(", order[1], ",", order[2], ",", order[3], ")")
  
  if (order[7] > 1 & sum(order[4:6]) > 0) 
    result <- paste0(result, "(", order[4], ",", order[5], 
                     ",", order[6], ")[", order[7], "]")
  if (!is.null(object$xreg)) {
    if (NCOL(object$xreg) == 1 & is.element("drift", nc)) 
      result <- paste(result, "with drift        ")
    else result <- paste("Regression with", result, "errors")
  }
  else {
    if (is.element("constant", nc) | is.element("intercept", nc)) 
      result <- paste(result, "with non-zero mean")
    else if (order[2] == 0 & order[5] == 0) 
      result <- paste(result, "with zero mean    ")
    else result <- paste(result, "                  ")
  }
  if (!padding) 
    result <- gsub("[ ]*$", "", result)
  return(result)
}

