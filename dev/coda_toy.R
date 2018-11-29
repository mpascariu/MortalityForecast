# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Thu Nov 22 15:48:29 2018
# --------------------------------------------------- #
remove(list = ls())
library(MortalityForecast)
library(MortalityLaws)
library(compositions)
library(forecast)


coda <- function(data){
  
  dx  <- t(data)
  ax  <- apply(dx, 2, mean)
  ax  <- ax/sum(ax)
  tdx <- sweep(acomp(dx), 2, ax, FUN = "-")
  cdx <- clr(acomp(tdx))
  
  S  <- svd(cdx)
  kt <- S$d[1] * S$u[, 1]
  bx <- S$v[,1]
  cf <- list(ax = as.numeric(ax), bx = as.numeric(bx), kt = as.numeric(kt))
  
  fv  <- clrInv(c(kt) %*% t(bx))
  fv  <- sweep(fv, 2, ax, FUN = "+")
  fdx <- unclass(t(fv/rowSums(fv)))
  dimnames(fdx) <- dimnames(data)
  
  out <- list(coefficients = cf, fitted.values = fdx)
  return(out)
}


pred_coda <- function(object, h) {
  C <- coef(object)
  M <- Arima(C$kt, order = c(0,1,0), include.drift = TRUE)
  kt <- forecast(M, h = h)
  
  # mean
  tdx <- clrInv(c(kt$mean) %*% t(C$bx))
  pv <- sweep(tdx, 2, C$ax, FUN = "+")
  pv <- unclass(t(pv))
  colnames(pv) <- 1:h
  
  out <- list(predicted.values = pv, kt = kt, kt.arima = M)
  return(out)
}


codaC <- function(data, B.data){
  
  # Average
  B <- coda(B.data)
  B.tdx <- with(coef(B), clrInv(c(kt) %*% t(bx)))
  
  # Deviation
  dx  <- t(data)
  ax  <- apply(dx, 2, mean)
  ax  <- ax/sum(ax)
  A.tdx <- sweep(acomp(dx), 2, ax, FUN = "-")
  tdx <- A.tdx - B.tdx
  cdx <- clr(acomp(tdx))
  
  S  <- svd(cdx)
  kt <- S$d[1] * S$u[, 1]
  bx <- S$v[,1]
  cf <- list(ax = as.numeric(ax), bx = as.numeric(bx), kt = as.numeric(kt))
  
  fv  <- clrInv(c(kt) %*% t(bx))
  fv  <- fv + B.tdx 
  fv  <- sweep(fv, 2, ax, FUN = "+")
  fdx <- unclass(t(fv/rowSums(fv)))
  dimnames(fdx) <- dimnames(data)
  
  out <- list(coefficients = cf, fitted.values = fdx, benchmark = B)
  return(out)
}


pred_codaC <- function(object, h) {
  B <- object$benchmark
  B.pred <- pred_coda(B, h = h)
  
  C <- coef(object)
  M <- Arima(C$kt, order = c(0,1,0), include.drift = TRUE)
  kt <- forecast(M, h = h)
  
  #mean
  B.tdx <- clrInv(c(B.pred$kt$mean) %*% t(coef(B)$bx))
  A.tdx <- clrInv(c(kt$mean) %*% t(C$bx))
  tdx <- A.tdx + B.tdx
  pv <- sweep(tdx, 2, C$ax, FUN = "+")
  pv <- unclass(t(pv))
  colnames(pv) <- 1:h
  
  out <- list(predicted.values = pv, kt = kt, kt.arima = M, benchmark = B.pred)
  return(out)
}




# ------------------------------------------
# Data

h = 5
x <- 0:100
y <- 1985:2016
dx.k <- HMD_male$dx[["USA"]]
dx.k <- dx.k[paste(x), paste(y)]
dx.k <- convertFx(x = x, data = dx.k, "dx", "dx", lx0 = 1)

dx.B <- HMD_female$dx[["USA"]]
dx.B <- dx.B[paste(x), paste(y)]
dx.B <- convertFx(x = x, data = dx.B, "dx", "dx", lx0 = 1)

#  Check Oeppen
M1 <- coda(dx.k)
M2 <- model_Oeppen(dx.k)
P1 <- pred_coda(M1, h)
P2 <- predict(M2, h, order = c(0, 1, 0), include.drift = T, jumpchoice = "fit")

plot(P1$predicted.values[, 5], log = "y", type = "l")
lines(P2$predicted.values[, 5], col = 2)

ee0 = convertFx(x, dx.k, from = "dx", to = "ex")
ee1 = convertFx(x, as.data.frame(P1$predicted.values), from = "dx", to = "ex")
ee2 = convertFx(x, as.data.frame(P2$predicted.values), from = "dx", to = "ex")

ee3 = convertFx(x, as.data.frame(P2$conf.intervals[[1]]), from = "dx", to = "ex")
ee4 = convertFx(x, as.data.frame(P2$conf.intervals[[4]]), from = "dx", to = "ex")

ee3[1,] - ee4[1,]

ee0[1, ]
ee1[1, ]
ee2[1, ]




# Check OeppenC
N1 <- codaC(data = dx.k, B.data = dx.B)
N2 <- model_OeppenC(data = dx.k, B.data = dx.B)
R1 <- pred_codaC(N1, h)
R2 <- predict(N2, h, order = c(0, 1, 0), include.drift = TRUE, jumpchoice = "fit")


e0 = convertFx(x, dx.k, from = "dx", to = "ex")
e1 = convertFx(x, as.data.frame(R1$predicted.values), from = "dx", to = "ex")
e2 = convertFx(x, as.data.frame(R2$predicted.values), from = "dx", to = "ex")
e3 = convertFx(x, as.data.frame(R2$conf.intervals[[1]]), from = "dx", to = "ex")
e4 = convertFx(x, as.data.frame(R2$conf.intervals[[4]]), from = "dx", to = "ex")

e3[1, ] - e4[1, ]

e0[1, ]
e1[1, ]
e2[1, ]

convertFx(x, as.data.frame(R2$conf.intervals[[4]]), from = "dx", to = "ex")[1,]

R2$kt

n = 65
plot(as.numeric(dx.k[n, ]), pch = 16)
lines(M1$fitted.values[n, ], col = 2)
lines(M2$fitted.values[n, ], col = 4)
lines(N1$fitted.values[n, ], col = 3, lwd = 2)
lines(N2$fitted.values[n, ], col = 5, lwd = 2)



C1 <- N1$coefficients
C2 <- N2$coefficients


par(mfrow = c(1, 3))
plot( C1$ax, type = "l", col = 2)
lines(C2$ax, col = 4)
plot( C1$bx, type = "l", col = 2)
lines(C2$bx, col = 4)
plot( C1$kt, type = "l", col = 2)
lines(C2$kt, col = 4)

# ------------------------------------------

X0 <- model_Oeppen(data = dx.k)
X1 <- model_OeppenC(data = dx.k, B.data = dx.k)
X2 <- model_OeppenC(data = dx.k, B.data = dx.B)

Y0 <- predict(X0, h = 10, jumpchoice = "actual")
Y1 <- predict(X1, h = 10, jumpchoice = "actual")
Y2 <- predict(X2, h = 10, jumpchoice = "actual")


x0 <- convertFx(x, data = Y0$predicted.values, from = "dx", to = "ex")
x1 <- convertFx(x, data = Y1$predicted.values, from = "dx", to = "ex")
x2 <- convertFx(x, data = Y2$predicted.values, from = "dx", to = "ex")

x0[1, ]
x1[1, ]
x2[1, ]

plot(P3$predicted.values[1,])
lines(X3$predicted.values[1,])





















































