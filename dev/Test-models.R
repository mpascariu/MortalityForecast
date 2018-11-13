# Marius D. Pascariu --- Tue Nov 13 14:59:44 2018 ------------------------------
remove(list = ls())
library(MortalityForecast)
library(MortalityLaws)


x  <- 0:100
y  <- 1985:2014
mx <- MortalityForecast.data$mx[paste(x), paste(y)]
dx <- MortalityForecast.data$dx[paste(x), paste(y)]

M0 <- fit_MRW(data = log(mx), x = x, y = y, include.drift = TRUE)
M1 <- fit_LeeCarter2(data = mx, x = x, y = y)
M2 <- fit_HyndmanUllah(data = mx, x = x, y = y)
M3 <- fit_Oeppen(data = dx, x = x, y = y)
M4 <- fit_MEM(data = dx, x = x, y = y, n = 5)

M0
M1
M2
M3
M4


LifeTable(x, mx = exp(fitted(M0)))$lt[102:110, ]
LifeTable(x, mx = fitted(M1))
LifeTable(x, mx = fitted(M2))
LifeTable(x, dx = fitted(M3))
LifeTable(x, dx = fitted(M4))$lt[102:110, ]

R0 <- residuals(M0)
R1 <- residuals(M1)
R2 <- residuals(M2)
R3 <- residuals(M3)
R4 <- residuals(M4)

R0
R1
R2
R3
R4

plot(R0)
plot(R1)
plot(R2)
plot(R3)
plot(R4)

P0 <- predict(M0, h = 20)
P1 <- predict(M1, h = 20)
P2 <- predict(M2, h = 20)
P3 <- predict(M3, h = 20)
P4 <- predict(M4, h = 20)

P0
P1
P2
P3
P4

ls(M0)
ls(M1)
ls(M2)
ls(M3)
ls(M4)


LifeTable(x, mx = exp(P0$predicted.values)) 
LifeTable(x, mx = P1$predicted.values)
LifeTable(x, mx = P2$predicted.values)
LifeTable(x, dx = P3$predicted.values)
LifeTable(x, dx = P4$predicted.values)









































