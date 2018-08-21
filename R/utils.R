


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







