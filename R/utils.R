# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Thu Nov 29 13:11:11 2018
# --------------------------------------------------- #


#' IF `a` is not NULL choose `a` else take `b`
#' @name AorB
#' @keywords internal
"%||%" <- function(A, B) {
  if (!is.null(A)) A else B  
}


#' Summary function - display head and tail in a single data.frame
#' The original code for this function was first written for 'psych' R package
#' here we have modified it a bit
#' @param x A matrix or data frame or free text
#' @param hlength The number of lines at the beginning to show
#' @param tlength The number of lines at the end to show
#' @param digits Round off the data to digits
#' @param ellipsis separate the head and tail with dots
#' @keywords internal
#' @author William Revelle (\email{revelle@@northwestern.edu})
#' @keywords internal
#' 
head_tail <- function(x, 
                      hlength = 4, 
                      tlength = 4, 
                      digits = 2, 
                      ellipsis = TRUE) {
  
  if (is.data.frame(x) | is.matrix(x)) {
    if (is.matrix(x)) x = data.frame(unclass(x))
    nvar <- dim(x)[2]
    dots <- rep("...", nvar)
    h    <- data.frame(head(x, hlength))
    t    <- data.frame(tail(x, tlength))
    
    for (i in 1:nvar) {
      if (is.numeric(h[1, i])) {
        h[i] <- round(h[i], digits)
        t[i] <- round(t[i], digits)
      } else {
        dots[i] <- NA
      }
    }
    out <- if (ellipsis) rbind(h, ... = dots, t) else rbind(h, t)
    
  } else {
    h   <- head(x, hlength)
    t   <- tail(x, tlength)
    out <- if (ellipsis) rbind(h, "...       ...", t) else as.matrix(rbind(h, t))
  }
  return(out)
}


#' Wide table to long table
#' @param data Wide table
#' @param x vector of ages
#' @param which.x subset x vector
#' @param which.y subset y vector
#' @inheritParams do.MortalityModels
#' @keywords internal
#' @export
wide2long <- function(data, 
                      x = NULL, 
                      y = NULL, 
                      which.x = NULL, 
                      which.y = NULL, 
                      ...) {
  
  x <- x %||% 1:nrow(data)
  y <- y %||% 1:ncol(data)
  which.x <- which.x %||% x
  which.y <- which.y %||% y
  
  D <- data %>% as.data.frame %>% gather(key = "y")
  D$y <- as.numeric(D$y)
  D$x <- x
  
  out <- D[D$x %in% which.x & D$y %in% which.y, ]
  return(out)
}


#' List containing wide tables to 1 long table
#' @inheritParams wide2long
#' @keywords internal
#' @export
wide.list.2.long.df <- function(data, 
                                x, 
                                y, 
                                which.x = NULL, 
                                which.y = NULL, 
                                ...) {
  
  fn  <- function(Z) wide2long(data = Z, x, y, which.x, which.y)
  B   <- lapply(data, fn)
  out <- NULL
  
  for (i in 1:length(B)) {
    Bi <- B[[i]]
    Bi$Name <- names(B)[i]
    out <- rbind(out, Bi)
  }
  return(out)
}


#' Return the name of the object as.character
#' @param x Any object.
#' @keywords internal
#' @export
whatsYourName <- function(x) {
  deparse(substitute(x))
}


#' Find ARIMA order of a time series
#' @param a A numerical vector.
#' @keywords internal
find_arima <- function(a) {
  W <- auto.arima(a)
  O <- arimaorder(W)
  D <- any(names(coef(W)) %in% "drift")
  return(list(order = O, drift = D))
}


#' Replace zero's in input matrix
#' 
#' If a matrix or data.frame contains death-rates equal to zero,
#' they will be replaced with positive values.
#' @param mx Matrix of death rates
#' @param method Method of replacement: \itemize{
#'  \item{"min"} -- Replace with minimum observed value in the input dataset;
#'  \item{"mult"} -- Replace using the multiplicative replacement strategy.
#' }
#' @param radix radix. Default: 1e5.
#' @examples 
#' x  <- 0:25
#' y  <- 2005:2016
#' mx <- HMD_male$mx$DNK[paste(x), paste(y)]
#' mx == 0
#' 
#' new.mx <- replace.zeros(mx)
#' new.mx == 0
#' 
#' sum(mx) == sum(new.mx)
#' (mx - new.mx)/mx * 100
#' @export
replace.zeros <- function(mx, method = c("min", "mult"), radix = 1e5) {
  method <- match.arg(method)
  
  if (method == "min") {
    mx[mx == 0] <- min(mx[mx > 0], na.rm = TRUE)/2
    out <- mx
  }
  
  if (method == "mult") {
    Dx <- mx * radix
    p <- sweep(Dx, 2, colSums(Dx), FUN = "/")
  
    for (i in 1:nrow(p)) {
      sdx <- 0.5 / sum(Dx[i, ], na.rm = TRUE)
      L <- as.numeric(p[i, ] == 0) * sdx
      p[i, ] <- p[i, ] + L
    }
  
    z <- sweep(p, 2, colSums(p), FUN = "/")
    out <- sweep(z, 2, colSums(Dx), FUN = "*") / radix
  }

  return(out)
}


#' Identify ARIMA model - internal function
#' @param object An object generate by Arima function
#' @param padding Logical.
#' @keywords internal
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
    
  } else {
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


# S3 - default ---------------------------------

#' Print function for mortality objects
#' @param x A mortality object;
#' @param ... Further arguments passed to or from other methods.
#' @keywords internal
print_default <- function(x, ...) {
  cat("\nFit  :", x$info$name)
  cat("\nModel:", x$info$formula)
  cat("\nCall : "); print(x$call)
  cat("Ages  in fit:", paste(range(x$x), collapse = " - "))
  cat("\nYears in fit:", paste(range(x$y), collapse = " - "))
  cat("\n")
}


#' @rdname print_default
print_predict_default <- function(x, ...) {
  cat("\nForecast:", x$info$name)
  cat("\nModel   :", x$info$formula)
  cat("\nCall    : "); print(x$call)
  cat("Ages  in forecast:", paste(range(x$x), collapse = " - "))
  cat("\nYears in forecast:", paste(range(x$y), collapse = " - "))
  cat("\n")
}


#' Extract Model Residuals
#' @param object A fitted mortality model
#' @inheritParams print_default
#' @keywords internal
#' @export
residuals_default <- function(object, ...){
  structure(class = "residMF", as.matrix(object$residuals))
}


