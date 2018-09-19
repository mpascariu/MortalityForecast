## ------------------------------------------------------------------- ##
##  R code to estimate and forecast adult mortality using the 
##  STAD model described in: Basellini U. and Camarda C.G. (2018), 
##  "Modelling and Forecasting Adult Age-at-death Distributions", 
##  Population Studies
##  
##  Authors: Ugofilippo Basellini & Carlo Giovanni Camarda
##  
##  sessionInfo() details:
##  
##  R version 3.4.2 (2017-09-28)
##  Platform: x86_64-apple-darwin15.6.0 (64-bit)
##  Running under: macOS High Sierra 10.13.5
##  
##  locale: en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
##
##  attached base packages:
##  splines  stats  graphics  grDevices  utils  datasets 
##  methods  base     
## 
##  other attached packages:
##  colorspace_1.3-2  vars_1.5-2  lmtest_0.9-35  urca_1.3-0  strucchange_1.5-1
##  sandwich_2.4-0  zoo_1.8-0  MASS_7.3-47  demography_1.20  forecast_8.3  
##  MortalitySmooth_2.3.4  lattice_0.20-35  svcm_0.1.2  Matrix_1.2-11
##
## ------------------------------------------------------------------- ##

## clean the workspace
rm(list = ls())

## load useful packages
library(MortalitySmooth)
library(demography)
library(vars)
library(colorspace)

## Function for smoothing mortality rates from deaths and exposures with 
## monotonic constraint. Age range can be greater than data, and monotonicity
## is needed for extrapolation
lmx_smooth <- function(x, y, e, w, ncoef, B){
  
  mx <- length(x) ## length of x range
  my <- length(y) ## length of data
  ## extend data if needed
  if (mx > my) {
    yA <- c(y, rep(99, mx - my))
    eA <- c(e, rep(99, mx - my))
    wA <- c(w, rep(0, mx - my))
    
  }else{
    yA <- y
    eA <- e
    wA <- w
  }
  ## starting value: glm with weights
  options(warn = -1)
  fit0 <- glm(round(yA) ~ x, offset = log(eA), weights = wA, family = poisson()) 
  options(warn = 0)
  # including monotonicity
  Dmon <- diff(diag(ncoef), diff = 1)
  wmon <- rep(0, nrow(Dmon))
  Wmon <- diag(wmon)
  kappa <- 10^4
  eps <- 0.1
  epsvec <- rep(eps, ncoef - 1)
  ## selecting lambda
  lambdas <- 10^seq(-1, 6, 0.5)
  nl <- length(lambdas)
  BICs <- numeric(nl)
  COEFs <- matrix(0, ncoef, nl)
  ## difference matrix
  D <- diff(diag(ncoef), diff = 2)
  tDD <- t(D) %*% D
  ## loop over the lambdas
  for (l in 1:nl) {
    P <- lambdas[l] * tDD
    ## starting linear predictor
    eta <- fit0$coef[1] + fit0$coef[2] * x
    
    for (it in 1:100) {
      ## monotonicity penalty
      Pmon0 <- t(Dmon) %*% Wmon %*% Dmon
      Pmon <- kappa * Pmon0
      vmon <- kappa * t(Dmon) %*% (wmon * epsvec)
      ## regression part
      mu <- exp(eta) * eA
      z <- wA * ((yA - mu)/mu + eta)
      z[is.na(z)] <- 0
      w <- c(wA * mu)
      BtWB <- t(B) %*% (w * B)
      BtWBpP <- BtWB + P + Pmon
      BtWz <- t(B) %*% (w * z)
      a <- solve(BtWBpP, BtWz + vmon)
      eta.old <- eta
      eta <- B %*% a
      ## update monotonicity weights
      diff.a <- diff(a) >= eps
      wmon <- rep(1, nrow(Dmon))
      wmon[diff.a] <- 0
      Wmon <- diag(wmon)
      ## convergence critera  
      deta <- max(abs((eta.old - eta)/eta.old))
      # cat(deta, "\n")
      if (it >= 4 & deta <= 10^-5) break   
    }
    
    ## deviance
    y1 <- yA
    mu1 <- mu
    y1[y1 == 0] <- 10^(-5)
    mu1[mu1 == 0] <- 10^(-5)
    dev <- 2 * sum(wA*(y1 * log(y1/mu1)), na.rm = TRUE)
    ## ED
    H <- solve(BtWBpP, BtWB)
    h <- diag(H)
    ed <- sum(h)
    ## BIC
    BICs[l] <- dev + log(sum(wA)) * ed
    ## fitted coef
    COEFs[, l] <- a
  }
  ## evaluate mortality at finer grid ~ hazard function
  pmin <- which.min(BICs)
  coef.hat <- COEFs[, pmin]
  out <- list(coef = coef.hat)
  return(out)
}

## function for aligning fx given a shifting parameter
fx_shift <- function(x, fx, shift, ndx = 25, deg = 3){
  ## length of x range
  xs <- x
  ms <- length(xs)
  delta <- round(diff(x)[1],4)
  ages.add <- 0
  ## augment basis
  if (shift > 0) ages.add = ceiling(shift + 5)
  if (shift < 0) ages.add = ceiling(abs(shift - 5))
  ## weights
  w <- c(rep(0,ages.add/delta),rep(1,ms),rep(0,ages.add/delta))
  ## augmented data
  xA <- c(rev(seq(from = xs[1] - delta, by = -delta, length = ages.add / delta)), 
          xs, 
          seq(from = xs[ms] + delta, by = delta, length = ages.add / delta))
  fA <- c(rep(999, ages.add / delta), fx, rep(999, ages.add / delta))
  ## new basis & B-splines parameters
  xl <- min(xA)
  xr <- max(xA)
  xmin <- round(xl - 0.01 * (xr - xl), 3)
  xmax <- round(xr + 0.01 * (xr - xl), 3)
  nb <- ndx + deg
  ## new B-splines bases
  BA <- MortSmooth_bbase(xA, xmin, xmax, ndx, deg)
  ## smoothing + extrapolation
  D <- diff(diag(nb), diff = 2)
  tDD <- t(D) %*% D
  P <- 10^-5 * tDD
  betasA <- solve(t(BA) %*% (w*BA) + P, t(BA) %*% (w * log(fA)))
  ## shifting
  ws.i <- xs - shift
  ## evaluating B on shifted x
  Bshift <- MortSmooth_bbase(ws.i, xl = xmin, xr = xmax, ndx, deg)
  ## shifted density
  fshift <- exp(Bshift %*% betasA)
  out <- fshift
  if (shift == 0) out <- fx
  return(out)
}

## function for obtaining coefficients of standard
coeff_stand <- function(x, fx, ages.add.l = 30, ages.add.r = 20, ndx = 30, deg = 3){
  ## Extrapolate age range to left and right + derive coefficients
  ## length of age range
  xs <- x
  delta <- round(diff(xs)[1],4)
  ms <- length(xs)
  ## weights
  w <- c(rep(0, ages.add.l / delta), rep(1, ms), rep(0, ages.add.r / delta))
  ## augmented data
  xA <- c(rev(seq(from = xs[1] - delta, by = -delta, length = ages.add.l / delta)), 
          xs, 
          seq(from = xs[ms] + delta, by = delta, length = ages.add.r / delta))
  fA <- c(rep(999, ages.add.l/delta), fx, 
          rep(999, ages.add.r/delta))
  ## new basis & B-splines parameters
  xl <- min(xA)
  xr <- max(xA)
  xmin <- round(xl - 0.01 * (xr - xl),3)
  xmax <- round(xr + 0.01 * (xr - xl),3)
  nb <- ndx + deg
  ## B-splines bases standard expanded
  BA <- MortSmooth_bbase(xA, xmin, xmax, ndx, deg)
  ## simple smoothing + extrapolation
  D <- diff(diag(nb), diff = 2)
  tDD <- t(D) %*% D
  P <- 10^-5 * tDD
  ## finding betas for the standard distribution expanded
  betasA <- solve(t(BA) %*% (w * BA) + P, t(BA) %*% (w * log(fA)))
  ## effective dimension
  tBWB <- t(BA) %*% (w*BA)
  tBWBpP <- tBWB + P 
  H <- solve(tBWBpP, tBWB)
  h <- diag(H)
  ed <- sum(h)
  ## return
  out <- list(betasA = betasA, ed = ed)
  return(out)
}

## function to compute mx from dx
mx_from_dx <- function(dx, ax = NULL){
  ## dimension of dx
  m <- length(dx)
  ## template vectors
  lx <- Lx <- mx <- rep(NA,m)
  ## set the radix of life table
  lx[1] <- sum(dx,na.rm = T)
  ## compute l(x+1)=l(x)-d(x) 
  for (i in 2:m) {
    lx[i] <- lx[i - 1] - dx[i - 1]
  }
  ## set ax = 1/2
  if (is.null(ax)) ax <- rep(1/2,m)
  ## compute Lx = l(x+1) + ax*dx
  Lx[-m] <- lx[-1] + ax[-m] * dx[-m]
  ## compute mx
  mx <- dx/Lx
  ## return mx value
  return(mx)
}

## Function to compute dx from mx 
dx_from_mx <- function(x,mx){
  ## length of age range
  xs <- x
  ms <- length(xs)
  delta <- round(diff(xs)[1],4)
  ## identity matrix and C matrix
  I <- diag(ms)
  C <- lower.tri(I, diag = TRUE)
  C[C == 1] <- -delta 
  ## compute dx
  dx <- mx * exp(C %*% mx)
  return(dx)
}


## optimization function for bL and bU
MLE_obj_FUN <- function(par, x, xA,
                        Mstand, shat, xlo, xup, coeff.stand, Dx, Ex, 
                        xmin, xmax, ndx, deg){
  ## starting b
  bL <- par[1]
  bU <- par[2]
  ## segment a linear transformation function
  ## below the mode
  aL <- Mstand - bL * (shat + Mstand)
  wL <- aL + bL * xlo
  ## above the mode
  aU <- Mstand - bU * (shat + Mstand)
  wU <- aU + bU * xup
  ## unique transformation function
  wb <- c(wL, wU)
  wb <- sort(wb)
  ## B-splines on transformed ages
  Bwb <- MortSmooth_bbase(x = c(wb), xmin, xmax, ndx = ndx, deg = deg)
  ## transformed density
  dwb <- as.vector(exp(Bwb %*% coeff.stand))
  dwb <- dwb[xA %in% x]
  if (sum(dwb) > 1) dwb <- dwb / sum(dwb)
  ## hazard
  eta <- log(mx_from_dx(dx = dwb))
  mu <- exp(eta)
  ## minimise minus the Log-Likelihood (maximise the LL)
  Lij <- -sum(Dx * eta[1:length(Dx)] - Ex * mu[1:length(Ex)], 
              na.rm = TRUE)
  return(Lij)
}

# ----------------------------------------------
# ----------------------------------------------
# ----------------------------------------------


fit_BaselliniCamarda <- function(FittingData, x, y, deg = 3) {
  
  xx <- min(x):120 ## expanded ages
  xs <- seq(min(xx), max(xx), by = 0.1)
  m  <- length(xx)
  ms <- length(xs)
  n  <- length(y)
  
  ## final starting data: actual death counts and exposures
  E <- FittingData$pop$female
  Z <- E * FittingData$rate$female
  
  # Mx = FittingData$rate$female
  # E <- Mx * 0 + 1
  # Z <- Mx
  
  ## B-splines bases
  xl <- min(xx)
  xr <- max(xx)
  xmin <- round(xl - 0.01 * (xr - xl), 3)
  xmax <- round(xr + 0.01 * (xr - xl), 3)
  ndx <- floor(m/3)
  nbx <- ndx + deg
  
  B <- MortSmooth_bbase(xx, xmin, xmax, ndx, deg)
  Bs <- MortSmooth_bbase(xs, xmin, xmax, ndx, deg)
  
  ## -- FITTING -------------
  ## smooth mortality for each y, forcing monotonicity at all ages
  w <- matrix(1, length(x), n) # weigths to avoid issue with zero
  w[Z < 0.00001] <- 0
  w[is.na(Z)] <- 0
  lMX.smooth <- matrix(NA, ms, n)
  for (i in 1:n) {
    lmx.smooth.coeff <- lmx_smooth(x = xx, y = Z[, i], e = E[, i], 
                                   w = w[, i], ncoef = nbx, B = B)$coef
    lMX.smooth[,i] <- Bs %*% lmx.smooth.coeff
  }
  
  ## compute density from smooth mortality rates
  FXs <- matrix(0, nrow = ms, ncol = n)
  for (i in 1:n) FXs[, i] = dx_from_mx(x = xs, mx = exp(lMX.smooth[, i]))
  
  ## compute modal age at death
  M <- xs[apply(FXs, 2, which.max)]
  
  ## ESTIMATE THE STANDARD DISTRIBUTION -----------------------------
  ## 1. Alignment procedure: get shifting parameter
  s <- M - M[1]
  ## derive aligned distributions
  FXs.align <- matrix(0, nrow = ms, ncol = n)
  for (i in 1:n) {
    FXs.align[, i] <- fx_shift(x = xs, fx = FXs[, i], shift = -s[i], 
                               ndx = ndx, deg = deg)
  }
  ## 2. Standard = mean of the aligned density
  FXallmean <- apply(FXs.align, 1, mean)
  FXstand <- FXallmean
  ## 3. Mode of the standard (extrapolate to left)
  Mstand <- M[1]
  
  ## 4. Coefficients of the standard: need to augment x-axis
  ages.add.l <- 40
  ages.add.r <- 30
  delta1 <- 1
  
  ## define new augmented x-axis 
  xA <- c(rev(seq(from = xx[1] - delta1, by = -delta1, length = ages.add.l / delta1)), 
          xx, 
          seq(from = xx[m] + delta1, by = delta1, length = ages.add.r / delta1))
  mA <- length(xA)
  
  ## new B-splines parameters on augmented axis
  xlA <- min(xA)
  xrA <- max(xA)
  xminA <- round(xlA - 0.01 * (xrA - xlA),3)
  xmaxA <- round(xrA + 0.01 * (xrA - xlA),3)
  ndxA <- floor((xA[mA] - xA[1])/3)
  
  ## augmented B-splines
  # BA <- MortSmooth_bbase(xA, xminA, xmaxA, ndx = ndxA, deg = deg)
  
  ## Derive coefficients of the standard
  Standard <- coeff_stand(x = xs, fx = FXstand, ndx = ndxA, deg = deg,
                          ages.add.l = ages.add.l, ages.add.r = ages.add.r)
  coeff_Stand <- Standard$betasA
  
  ## MLE ESTIMATION -------------------------------------------------
  ## Break point of the x axis for each year 
  PLO <- PUP <- XLO <- XUP <- list()
  for (i in 1:n) {
    PLO[[i]] <- which(xA <= floor(M[i]))
    PUP[[i]] <- which(xA > floor(M[i]))
    XLO[[i]] <- xA[PLO[[i]]]
    XUP[[i]] <- xA[PUP[[i]]]
  }
  
  ## empty vectors and matrices to store results 
  ss <- bsLO <- bsUP <- numeric(n)
  DXstad <- matrix(NA, m, n)
  
  ## Estimation
  for (i in 1:n) {
    # conv.stad <- FALSE
    ## starting values
    if (i == 1) {
      ## start from 1 if first year  
      start.value <- c(1,1)
    }else{
      ## start from previously estimated pars
      start.value <- c(bsLO[i - 1], bsUP[i - 1])
    }
    ## MLE
    opt <- optim(par = start.value, fn = MLE_obj_FUN, x = xx, xA = xA,
                 Mstand = Mstand, shat = s[i], xlo = XLO[[i]], xup = XUP[[i]],
                 coeff.stand = coeff_Stand, Dx = Z[,i], Ex = E[,i],
                 xmin = xminA, xmax = xmaxA, ndx = ndxA, deg = deg)
    if (opt$convergence != 0) break
    
    ## assign 
    shat <- s[i]
    bLhat <- opt$par[1]
    bUhat <- opt$par[2]
    
    ## compute dx:
    ## segment a linear transformation function
    ## below the mode
    aL <- Mstand - bLhat * (shat + Mstand)
    wL <- aL + bLhat * XLO[[i]]
    ## above the mode
    aU <- Mstand - bUhat * (shat + Mstand)
    wU <- aU + bUhat * XUP[[i]]
    ## unique transformation function
    wb <- c(wL, wU)
    wb <- sort(wb)
    ## B-splines on transformed ages
    Bwb <- MortSmooth_bbase(x = c(wb), xminA, xmaxA, ndx = ndxA, deg = deg)
    ## transformed density
    fwb <- as.vector(exp(Bwb %*% coeff_Stand))
    fwb <- fwb[xA %in% xx]
    ## check if density is greater than one
    if (sum(fwb) > 1) fwb <- (fwb)/sum(fwb)
    ## save parameters, density and mx
    ss[i] <- shat
    bsLO[i] <- bLhat
    bsUP[i] <- bUhat
    DXstad[,i] <- fwb
  }
  
  cf <- data.frame(s = ss, bl = bsLO, bu = bsUP)
  dimnames(DXstad) <- list(xx, y)
  
  out <- list(x = x, y = y,
              coefficients = cf, fitted.values = DXstad[seq_along(x), ], 
              Mstand = Mstand, xA = xA,
              PLO = PLO, PUP = PUP, XLO = XLO, XUP = XUP,
              coeff_Stand = coeff_Stand)
  return(out)
}

## -- READING THE DATA -------------

## Read mortality data   
cou <- "SWE"   ## choose any country from HMD
user <- "zpascariu@outlook.com"   ## set your HMD credentials
pass <- "1513208587"   ## set your HMD credentials
FullData <- hmd.mx(cou, username = user, password = pass)


## original ages
x <- 30:110
## years
y <- 1980:2014

## subset age and year-specific data
FittingData <- extract.years(FullData, years = y)
FittingData <- extract.ages(FittingData, ages = x)

BC <- fit_BaselliniCamarda(FittingData, x, y)


# apply(BC$fitted, 2, sum)










