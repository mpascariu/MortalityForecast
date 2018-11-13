# Mon Aug 20 16:02:25 2018 --------- Marius D. Pascariu ---


#' IF `a` is not NULL choose `a` else take `b`
#' @name AorB
#' @keywords internal
"%||%" <- function(A, B) {
  if (!is.null(A)) A else B  
}


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
head_tail <- function(x, hlength = 4, tlength = 4, 
                      digits = 2, ellipsis = TRUE) {
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
  }
  else {
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
#' @inheritParams doMortalityModels
#' @keywords internal
#' @export
wide2long <- function(data, x = NULL, y = NULL, 
                      which.x = NULL, which.y = NULL, ...) {
  if (is.null(x)) x = 1:nrow(data)
  if (is.null(y)) y = 1:ncol(data)
  if (is.null(which.x)) which.x = x
  if (is.null(which.y)) which.y = y
  
  D <- gather(data, key = "y")
  D$x <- x
  D$y <- as.numeric(D$y)
  out <- D[D$x %in% which.x & D$y %in% which.y, ]
  return(out)
}


#' List containing wide tables to 1 long table
#' @inheritParams wide2long
#' @keywords internal
#' @export
wide.list.2.long.df <- function(data, x, y, which.x = NULL, which.y = NULL, ...) {
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


#' Return the name of the object as.charater
#' @param x Any object.
#' @keywords internal
#' @export
whatsYourName <- function(x) {
  deparse(substitute(x))
}


# S3 - default ---------------------------------

#' Print function for mortality models
#' @param x A fitted mortality model
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

#' Default print function for predict methods
#' @param x An predict object;
#' @param ... Further arguments passed to or from other methods.
#' @keywords internal
print_predict_default <- function(x, ...) {
  cat("\nForecast:", x$info$name)
  cat("\nModel   :", x$info$formula)
  cat("\nCall    : "); print(x$call)
  cat("Ages  in forecast:", paste(range(x$x), collapse = " - "))
  cat("\nYears in forecast:", paste(range(x$y), collapse = " - "))
  cat("\n")
}

#' Residuals function for mortality models
#' @param object A fitted mortality model
#' @inheritParams print_default
#' @export
residuals_default <- function(object, ...){
  structure(class = "residMF", as.matrix(object$residuals))
}

