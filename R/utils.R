

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







