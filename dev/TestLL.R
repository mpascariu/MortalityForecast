# --------------------------------------------------- #
# Author: Marius D. Pascariu
# Last update: Sat Nov 17 10:32:09 2018
# --------------------------------------------------- #
remove(list = ls())

library(MortalityForecast)
library(MortalityLaws)
library(forecast)

x <- 0:89
y <- 1985:2014
h <- 20
k <- c("GBRTENW","DNK","NLD","SWE","USA")

mx0 <- HMD_male$mx[k]
mx1 <- lapply(mx0, as.matrix)
mx2 <- lapply(mx1, replace_value_zero)
lapply(mx2, dim)
mx3 <- apply(simplify2array(mx2), 1:2, FUN = prod)^(1/(length(k)))

mx.B <- mx3[paste(x), paste(y)]
mx.k <- HMD_male$mx$USA[paste(x), paste(y)]


dim(mx.k)
dim(mx.B)

r <- 1
nx <- 1
LC <- LC.fun(mx = mx.k, h = h)
LL <- LL.fun(mx = mx.k, mx.av = mx.B, h = h)

LC2 <- fit_LeeCarter2(mx.k, x = x, y = y)
LL2 <- fit_LiLee(data = mx.k, benchmark = mx.B, x = x, y = y)

par(mfrow = c(1, 3))
plot(LC$ax, type = "l", lwd = 3)
lines(LL$ax, col = 2, lwd = 2)
lines(coef(LC2)$ax, col = 4, lwd = 2)
lines(coef(LL2)$ax, col = 3, lwd = 2)

plot(LC$bx, type = "l", ylim = range(LC$bx, LL$bx), lwd = 2)
lines(LL$bx, col = 2, lwd = 2)
lines(LL$bx, col = 2, lwd = 2)
lines(coef(LC2)$bx, col = 4, lwd = 2)
lines(coef(LL2)$bx, col = 3, lwd = 2)

plot(LC$kt, type = "l", ylim = range(LC$kt, LL$kt), lwd = 2)
lines(LL$kt, col = 2, lwd = 2)
lines(coef(LC2)$kt, col = 4, lwd = 2)
lines(coef(LL2)$kt, col = 3, lwd = 2)



plot(LC$ex[1, ])
lines(LL$ex[1, ])

LL$ex

LC$ex[, 1]
