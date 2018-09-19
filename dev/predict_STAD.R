set.seed(2018) 


## -- FORECASTING -------------
## uncertainty interval & seed for reproducibility

predict.fit_BaselliniCamarda <- function(object, h, level = 80, nsim = 200, ...) {
  
  x <- object$x
  y <- object$y
  n <- length(y)
  y_for <- max(y) + 1:h
  
  xx <- min(x):120 ## expanded ages
  m <- length(xx)
  ss <- object$coefficients$s
  bsLO <- object$coefficients$bl
  bsUP <- object$coefficients$bu
  deg <- 3
  xA  <- object$xA
  PLO <- object$PLO
  PUP <- object$PUP
  XLO <- object$XLO
  XUP <- object$XUP
  mA  <- length(xA)
  Mstand <- object$Mstand
  coeff_Stand <- object$coeff_Stand
  
  level.p <- level/100
  ## create ts dataframe for s, bL and bU
  s.ts <- ts(ss, start = y[1])
  df.par <- data.frame(bsLO = bsLO, bsUP = bsUP)
  df.par.ts <- ts(df.par, start = y[1])
  
  ## Forecasting s: univariate best ARIMA ------------
  ## Testing for unit root in original series with augmented DF test
  df.test.s <- ur.df(s.ts, lags = 2, 'trend')
  ds <- if (df.test.s@teststat[1] > df.test.s@cval[1,][1]) 1 else 0
  
  ## ts model for s
  s.mod1 <- auto.arima(s.ts, d = ds, max.p = 3, max.q = 3, trace = FALSE)
  # pred.s <- forecast(s.mod1, h = h, level = level)
  
  ## Forecasting bL and bU: VAR model 
  var <- VAR(df.par.ts, p = 1, type = "both")
  # pred.var <- forecast(var, h = h, level = level)
  
  ## BOOTSTRAP FOR CONSTRUCTING PI  --------------------
  ## define three bootstapping matrices
  bootstrap.s <- bootstrap.bsLO <- bootstrap.bsUP <- 
    matrix(NA, nrow = h, ncol = nsim)
  
  ## bootsrap S
  for (i in 1:nsim) {
    ## generate simulation with bootsrapping
    s.sim <- simulate(s.mod1, nsim = h, future = TRUE, bootstrap = TRUE)
    ## derive the bootsrap values
    bootstrap.s[,i] <- s.sim
  }
  
  ## bootsrap BL and BU
  ## estimated residuals
  res.bL <- resid(var)[, 1]
  res.bU <- resid(var)[, 2]
  for (i in 1:nsim) {
    ## var dimension
    n.var <- n - 1
    ## sampling residuals
    res.bL.S <- sample(res.bL, n.var, replace = TRUE)
    res.bU.S <- sample(res.bU, n.var, replace = TRUE)
    ## fictitious responses
    bL.S <- c(fitted(var)[, 1] + res.bL.S)
    bU.S <- c(fitted(var)[, 2] + res.bU.S)
    ## re-fit the model
    df.par.S <- data.frame(bsLO = bL.S, bsUP = bU.S)
    df.par.ts.S <- ts(df.par.S, start = y[1 + n.var])
    var.S <- VAR(df.par.ts.S, p = 1, type = "both")
    ## forecast
    pred.var.S <- forecast(var.S, h = h, level = level)
    ## derive the bootsrap values
    bootstrap.bsLO[,i] <- pred.var.S$forecast$bsLO$mean
    bootstrap.bsUP[,i] <- pred.var.S$forecast$bsUP$mean
  }
  
  ## matrix and arrays to store bootsrap results
  dwb.boot <- array(NA, c(m, nsim, h))
  
  ## for each simulation and year, compute age-at-death distribution
  for (j in 1:h) {
    for (i in 1:nsim) {
      ## get values of s, bsLO, bsUP for each forecasted year
      s.boot <- bootstrap.s[j,i]
      bsLO.boot <- bootstrap.bsLO[j,i]
      bsUP.boot <- bootstrap.bsUP[j,i]
      ## compute new mode
      M.boot <- s.boot + Mstand
      ## divide forecast distribution in two pieces
      PLO[[n + j]] <- which(xA <= floor(M.boot))
      PUP[[n + j]] <- which(xA > floor(M.boot))
      XLO[[n + j]] <- xA[PLO[[n + j]]]
      XUP[[n + j]] <- xA[PUP[[n + j]]]
      
      ## segment a linear transformation function
      ## below the mode of the forecast year
      aL <- Mstand - bsLO.boot * (s.boot + Mstand)
      wL <- aL + bsLO.boot * XLO[[n + j]]
      ## above the mode of the forecast year
      aU <- Mstand - bsUP.boot * (s.boot + Mstand)
      wU <- aU + bsUP.boot * XUP[[n + j]]
      ## unique transformation function
      wbhat <- c(wL, wU)
      wbhat <- sort(wbhat)
      ## evaluating B on shifted x
      xlA <- min(xA)
      xrA <- max(xA)
      xminA <- round(xlA - 0.01 * (xrA - xlA), 3)
      xmaxA <- round(xrA + 0.01 * (xrA - xlA), 3)
      ndxA  <- floor((xA[mA] - xA[1])/3)
      Bs.fore <- MortSmooth_bbase(wbhat, xl = xminA, xr = xmaxA, ndx = ndxA, deg)
      ## shifted density
      fwb.boot <- exp(Bs.fore %*% coeff_Stand)
      fwb.boot <- fwb.boot[xA %in% xx]
      ## check if density is greater than one
      if (sum(fwb.boot) > 1) fwb.boot <- fwb.boot / sum(fwb.boot)
      ## save density,
      dwb.boot[, i, j] <- fwb.boot
    }
  }
  
  ## derive median and PI of distribution
  dwb.fore.MEAN <- dwb.fore.UPS <- dwb.fore.LOS <- matrix(NA, nrow = m, ncol = h)
  for (j in 1:h) {
    dwb.fore.MEAN[,j] <- apply(dwb.boot[,,j], 1, median)
    dwb.fore.UPS[,j] <- apply(dwb.boot[,,j], 1, quantile, 
                              prob = 1 - (1 - level.p)/2, na.rm = TRUE)
    dwb.fore.LOS[,j] <- apply(dwb.boot[,,j], 1, quantile, 
                              prob = (1 - level.p) / 2, na.rm = TRUE)
  }
  
  dimnames(dwb.fore.MEAN) = dimnames(dwb.fore.UPS) = dimnames(dwb.fore.LOS) <-
    list(xx, y_for)
  ci <- list(upper = dwb.fore.UPS[seq_along(x),], 
             lower = dwb.fore.LOS[seq_along(x),])
  out <- list(predicted.value = dwb.fore.MEAN[seq_along(x),], 
              conf.intervals = ci)
  return(out)
}


PBC <- predict.fit_BaselliniCamarda(BC, h = 26, level = 80)










