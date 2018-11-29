
#' Estimate density function based on Maximum Entropy method (MaxEnt)
#' 
#' Function for estimating the probability density function starting from 
#' a limited number of raw statistical moments. The method is based on 
#' Maximum Entropy method \insertCite{mead1984}{MortalityForecast}.
#' @param data A numeric vector or a data.frame/matrix with raw moments.
#' @param x A numeric vector with ages. The density frequencies are 
#' returned only for the specified ages. 
#' @param omega The maximum age to be used in MaxEnt estimation. The 
#' algorithm is reliable and robust when the underlying distribution is closing 
#' on the right-hand side. In the case of age-at-death distribution, \code{omega} 
#' should be an old age; otherwise the MaxEnt might fail to converge. 
#' Numeric scalar.
#' @inheritParams model_MEM
#' @return An object with the estimated densities.
#' @author Adam Lenart and Marius D. Pascariu.
#' @references \insertAllCited{}
#' @examples 
#' # Example 1 -- simple case ---------------
#' 
#' x <- 0:120
#' # Raw moments M0 - M7
#' raw_moments <- c(1, 68.75099, 4991.724, 371531.9, 28199680, 
#'                  2176435499, 170477697491, 1.353288e+13)
#' ed <- findDensity(raw_moments, x)
#' 
#' # Example 2 -- estimat for entire table --
#' 
#' x   <- 0:110
#' y   <- 1965:2016
#' dx  <- HMD_male$dx$GBRTENW[paste(x), paste(y)]
#' mom <- findMoments(dx, x, y, n = 7)$raw.moments
#' 
#' dx.hat <- findDensity(mom, x)
#' 
#' # Visual check --
#' yr = "1990"
#' observed.dx <- dx[, yr]/sum(dx[, yr])
#' estimated.dx <- dx.hat$density[, yr]
#' plot(observed.dx, pch = 16, main = yr)
#' lines(estimated.dx, col = 2)
#' legend("topleft", legend = c("observed", "estimated"), col = 1:2, 
#'        pch = c(16, NA), lty = c(NA, 1), lwd = 2, bty = "n")
#' @export
findDensity <- function(data, 
                        x, 
                        omega = 110, 
                        verbose = FALSE) {
  
  if (is.matrix(data) || is.data.frame(data)) {
    # If matrix or data.frame do the for-loop below
    
    if (verbose) { # Set progress bar
      pb <- startpb(0, ncol(data))
      on.exit(closepb(pb)); setpb(pb, 1)
    }
    
    N <- nrow(data)
    P <- matrix(NA, nrow = length(x), ncol = N)
    M <- L <- I <- NULL
    
    for (k in 1:N) { # THIS LOOP
      maxEnt <- findDensity(data[k, ], x, omega)
      P[, k] <- maxEnt$density
      M      <- rbind(M, maxEnt$fitted.raw.moments)
      L      <- rbind(L, maxEnt$lambda)
      I      <- c(I, maxEnt$iterations)
      if (verbose) setpb(pb, k + 1)
    }
    names(I) = rownames(L) = rownames(M) <- rownames(data)
    dimnames(P) <- list(x = x, y = rownames(data))
    out <- list(x = x, 
                density = P, 
                input.raw.moments = data,
                fitted.raw.moments = M, 
                lambda = L, 
                iterations = I)
    
  } else {

    # THIS IS THE MaxEnt ALGORITHM
    In     <- data[-1]
    N      <- length(In)
    L      <- max(omega, max(x)) - min(x)                     
    mom2   <- In/L^(1:N)               # scale moments
    lambda <- numeric(N)
    xx     <- seq(0, 1, length.out = (L + 1))
    x0     <- NULL
    for (i in 1:N) {
      x0 <- as.matrix(cbind(x0, xx^i)) 
    }
    k      <- 0
    a      <- 1
    
    while (all(round(a, 7) == 0) != TRUE ) {
      fx         <- exp(-apply(x0, 1, FUN = function(z) sum(z * lambda)) ) # step 1
      px         <- fx / sum(fx)                                           # step 2-3
      mom.est    <- sapply(X = 1:(2*N), FUN = function(n) sum(px*(xx^n)) ) # step 4
      mom.star   <- mom.est[1:N]                                           # step 5
      mom.star.t <- mom.star * L^(1:N)
      
      W <- matrix(NA, ncol = N, nrow = N)                         # step 6
      for (i in 1:N) {
        for (j in 1:N) {
          W[i, j] <- mom.est[i + j] - mom.est[i]*mom.est[j]
        }
      }
      
      invW   <- solve(W)         # step 7
      nu     <- mom2 - mom.star  # step 8
      a      <- invW %*% nu      # step 9
      lambda <- lambda - a       # step 10
      k      <- k + 1
    }
    
    lambda <- c(0, as.numeric(t(lambda)))
    mom <- c(1, mom.star.t)
    names(data) = names(mom) = names(lambda) <- paste0("M", 0:length(In))
    
    out <- list(x = x, 
                density = px[1:length(x)], 
                input.raw.moments = data, 
                fitted.raw.moments = mom, 
                lambda = lambda, 
                iterations = as.numeric(k - 1))
  }
  out$call <- match.call()
  out <- structure(class = "findDensity", out)
  return(out)
}

