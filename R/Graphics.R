

#' @keywords internal
#' @export
#' 
wide2long <- function(data, x = NULL, filter.x = NULL) {
  if (is.null(x)) x = 1:nrow(data)
  if (is.null(filter.x)) filter.x = x
  
  D <- gather(data, key = "y")
  D$x <- x
  D$y <- as.numeric(D$y)
  out <- D[D$x %in% filter.x, ]
  return(out)
}


#' @keywords internal
#' @export
#' 
wide.list.2.long.df <- function(data, x, filter.x = NULL, ...) {
  fn  <- function(Z) wide2long(data = Z, x, filter.x)
  B   <- lapply(data, fn)
  out <- NULL
  
  for (i in 1:length(B)) {
    Bi <- B[[i]]
    Bi$Name <- names(B)[i]
    out <- rbind(out, Bi)
  }
  return(out)
}


#' @keywords internal
#' @export
#'
whatsYourName <- function(x) {
  deparse(substitute(x))
}







