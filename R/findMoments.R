#' Compute Statistical Moments
#' 
#' @description Compute raw and central statistical moments of a distribution of deaths.
#' @details 
#' If \eqn{f(x)} is a probability density function, then \eqn{\mu(n)} is called 
#' the n-th moment of the probability distribution, where: \cr
#' \deqn{\mu_{n} = \int_{-\infty}^{\infty} (x-c)^n f(x)dx.}{\mu(n) = \int (x-c)^n f(x)dx.}
#' @param data A data.frame or matrix containing 3 columns: Year, Age, Dx.
#' @param x Vector of ages.
#' @param n The maximum order of the moments to be computed. 
#' The order should be at least 2.  
#' @param na.rm Logical value. If \code{TRUE}, remove \code{NA} values. 
#' Otherwise, keep \code{NA} values.
#' @inheritParams fitMaxEntMortality
#' @return An object containing:
#' \itemize{
#' \item \code{central.moments} --- Moments about the mean (c = mean). The zeroth moment is 
#' the total probability (i.e. one), the first moment is zero, 
#' the second central moment is the variance, the third central moment 
#' is the skewness (with normalization), and the fourth central moment 
#' is the kurtosis (with normalization).
#' \item \code{raw.moments} --- Moments about zero (c = 0). The first moment is the mean.
#' \item \code{normalized.moments} --- Normalized moments. 
#' } 
#' @seealso 
#' \code{\link{convertMoments}}
#' \code{\link[moments]{all.moments}} 
#' \code{\link[moments]{raw2central}}
#' @examples 
#' x  <- 0:110
#' y  <- 1960:2016
#' dx <- MortalityForecast.data$dx
#' findMoments(data = dx, x = x, y = y, n = 4)
#' @export
findMoments <- function(data, x, y = NULL, n, na.rm = TRUE) {
  N  <- ncol(data)
  rM <- matrix(NA, nrow = N, ncol = n + 1)
  if (is.null(y)) y <- 1:N
  dimnames(rM) <- list(y, paste0("M", 0:n))
  
  w   <- mean(diff(x)) # width of the age interval
  ax1 <- if (min(x) == 0) w/5 else w/2  
  a   <- as.numeric(x) + c(ax1, rep(w/2, length(x) - 1))
  if (x[1] > 0) {
    a <- a - min(a)
  }
  
  data <- as.matrix(data)
  f    <- 1e5 / colSums(data)
  D    <- round(sweep(data, 2, f, FUN = "*"))
  for (i in 1:N) {
    vect <- rep(a, D[, i])
    rM[i, ] <- moments::all.moments(x = vect, order.max = n, na.rm = na.rm)
  }
  
  cM <- convertMoments(rM, from = "raw", to = "central")
  nM <- convertMoments(rM, from = "raw", to = "normalized")
  dimnames(cM) <- dimnames(rM)
  
  out <- list(raw.moments = rM, central.moments = cM, normalized.moments = nM)
  return(out)
}


#' Convert Statistical Moments
#' 
#' Transform the raw, central or normalized statistical moments between them.
#' @param data data.frame containing statistical moments.
#' @param from What type of statistical moments do we have in input \code{data}?
#' Three types of moments are accepted: \code{"raw", "central", "normalized"}.
#' @param to What type of statistical moments do you want to obtain? 
#' Three types of moments can be obtained: \code{"raw", "central", "normalized"}.
#' @param eta A numeric vector of the expected values. This is required ONLY 
#' is we convert central moments into raw or normalized. Default: \code{NULL}.
#' 
#' @details Wikipedia: In probability theory and statistics, 
#' the standardized moment of a probability distribution is a moment 
#' (normally a higher degree central moment) that is normalized. 
#' The normalization is typically a division by an expression of the 
#' standard deviation which renders the moment scale invariant. 
#' This has the advantage that such normalized moments differ only in other 
#' properties than variability, facilitating e.g. comparison of shape of 
#' different probability distributions.
#' @examples 
#' # raw moments
#' RM <- c(1, 68.75099, 4991.724, 371531.9, 28199680,  
#'         2176435499, 170477697491)
#' 
#' 
#' CM1 <- convertMoments(RM, from = "raw", to = "central")    # raw to central
#' NM1 <- convertMoments(RM, from = "raw", to = "normalized") # raw to normalized
#' 
#' CM2 <- convertMoments(NM1, from = "normalized", to = "central")
#' RM2 <- convertMoments(NM1, from = "normalized", to = "raw")
#' 
#' RM3 <- convertMoments(CM2, from = "central", to = "raw", eta = 68.75099)
#' NM3 <- convertMoments(CM2, from = "central", to = "normalized", eta = 68.75099)
#' 
#' # The resulted error following multiple conversions is negligible
#' sum(RM - RM3)
#' @export
convertMoments <- function(data, 
                           from = c("raw", "central", "normalized"), 
                           to = c("raw", "central", "normalized"), 
                           eta = NULL) {
  if (is.vector(data)) data = matrix(data, nrow = 1)
  from <- match.arg(from)
  to <- match.arg(to)
  out <- data
  
  if (from == "raw") out <- convertRM(data, to)
  if (from == "central") out <- convertCM(data, to, eta)
  if (from == "normalized") out <- convertNM(data, to)
  
  return(out)
}


#' Convert Central Moments to Raw or Normalized moments
#' @param data Data.frame with central moments.
#' @inheritParams convertMoments
#' @keywords internal
convertCM <- function(data, to = c("raw", "normalized"), eta) {
  cM  <- data 
  n   <- ncol(cM)
  Mnames <- paste0("M", 1:n - 1)
  rM <- t(moments::central2raw(t(cM), eta))
  colnames(rM) = colnames(cM) <- Mnames
  out <- rM
  
  if (to == "normalized") {
    nM <- convertRM(rM, to)
    out <- nM
  }
  rownames(out) <- rownames(cM)
  return(out)
}


#' Convert Normalized Moments to Central or Raw moments
#' @param data Data.frame with normalized moments.
#' @inheritParams convertMoments
#' @keywords internal
convertNM <- function(data, to = c("central", "raw")) {
  nM <- data                                 # Normalized moments
  n  <- ncol(nM)
  if (n < 3L) {
    stop("'data' must have at least 3 columns. That is M0, M1 and M2.", call. = F)
  }
  Mnames <- paste0("M", 1:n - 1)
  colnames(nM) <- Mnames
  
  # Compute Mixed Moments
  mM <- nM
  if (n > 3) {
    p <- seq(1.5, 0.5 * (n - 1), by = 0.5)
    
    for (i in 4:n) {
      mM[, i] <- nM[, i] * nM[, 'M2']^p[i - 3]
    }
  }
  colnames(mM) <- Mnames[1:n]
  
  cM <- mM
  cM[, "M1"] <- 0
  out <- cM
  
  if (to == "raw") {
    rM <- t(moments::central2raw(t(cM), eta = mM[, "M1"]))
    dimnames(rM) <- dimnames(mM)
    out <- rM
  }
  # Exit
  return(out)
} 


#' Convert Raw Moments to Central or Normalized moments
#' @param data Data.frame with raw moments.
#' @inheritParams convertMoments
#' @keywords internal
convertRM <- function(data, to = c("central", "normalized")) {
  rM <- data                                 # Raw moments
  n  <- ncol(rM)
  if (n < 3L) {
    stop("'data' must have at least 3 columns. That is M0, M1 and M2.", call. = F)
  }
  
  mM <- cM <- t(moments::raw2central(t(rM)))  # Central moments
  mM[, 2]  <- rM[, 2]                         # Mixed moments
  colnames(rM) = colnames(cM) = colnames(mM) <- paste0("M", 1:n - 1)
  out <- cM
  
  if (to == "normalized") {
    # Compute Normalised Moments
    nM <- mM
    if (n > 3) {
      p <- seq(1.5, 0.5 * (n - 1), by = 0.5)
      
      for (i in 4:n) {
        nM[, i] <- mM[, i] / mM[, 'M2']^p[i - 3] 
      }
    }
    
    # Assign names
    Mnames <- c('M0', 'Mean (M1)', 'Variance (M2)', 'Skewness (M3)',
                'Kurtosis (M4)', 'Hyperskewness (M5)', 'Hyperflatness (M6)')
    if (n > 7) {
      for (k in 8:n) Mnames[k] = paste0('M', k - 1) 
    }
    colnames(nM) <- Mnames[1:n]
    out <- nM
  }
  
  return(out)
} 


