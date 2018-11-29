# --------------------------------------------------- #
# Author: Marius D. Pascariu
# Last update: Fri Nov 16 13:09:05 2018
# --------------------------------------------------- #
remove(list = ls())
library(MortalityForecast)
library(MortalityLaws)

x <- 0:89
y <- 1985:2016
h <- 20
k <- c("GBRTENW","DNK","NLD", "SWE","USA")

mx0 <- HMD_male$mx[k]
mx1 <- lapply(mx0, as.matrix)
mx2 <- lapply(mx1, replace_value_zero)
all(sapply(mx2, dim) == dim(mx2[[1]]))
mx3 <- apply(simplify2array(mx2), 1:2, FUN = function(w) prod(w, na.rm = TRUE)^(1/(length(w))))

mx.B <- mx3[paste(x), paste(y)]

mx.k <- HMD_male$mx[["DNK"]]
mx.k <- replace_value_zero(mx.k)
mx.k <- mx.k[paste(x), paste(y)]



LC <- fit_LeeCarter2(data = mx.k, x = x, y = y)
LL <- fit_LiLee(data = mx.k, benchmark = mx.B, x = x, y = y)



J = "actual"
P1 <- predict(LC, h = h, order = c(0,1,0), include.drift = TRUE, jumpchoice = J)
P2 <- predict(LL, h = h, order = c(0,1,0), include.drift = TRUE, jumpchoice = J)
P1$kt

observed.ex <- convertFx(x, data = mx.k, from = "mx", to = "ex")
observed.ex2 <- convertFx(x, data = LC$observed.values, from = "mx", to = "ex")
fx1 <- convertFx(x, data = fitted(LC), from = "mx", to = "ex")
fx2 <- convertFx(x, data = fitted(LL), from = "mx", to = "ex")

ex1 <- convertFx(x, data = P1$predicted.values, from = "mx", to = "ex")
ex2 <- convertFx(x, data = P2$predicted.values, from = "mx", to = "ex")



c1 = convertFx(x, data = P1$predicted.values[, 1], from = "mx", to = "ex")
c2 = convertFx(x, data = LC$observed.values[, "2016"], from = "mx", to = "ex")
c1 - c2

age <- "0"
e0 <- as.numeric(observed.ex[age, ])
e1 <- as.numeric(ex1[age, ])
e2 <- as.numeric(ex2[age, ])
f1 <- as.numeric(fx1[age, ])
f2 <- as.numeric(fx2[age, ])

plot(y, e0, pch = 16, ylim = range(e0, e1, e2), xlim = range(y, P1$y, P2$y))
points(y, observed.ex2[age, ], pch = 16, col = 4)
lines(P1$y, e1, col = 3)
lines(P2$y, e2, col = 2)
lines(y, f1, col = 3)
lines(y, f2, col = 2)


