# Tue Aug 14 13:51:37 2018 --------- Marius D. Pascariu ---
remove(list = ls())

library(MortalityAccuracy)
library(grid)
library(gridExtra)

x = 0:98
y1 = 1980:1999
y2 = 2000:2016
y  = c(y1, y2)

D <- dxForecast::dxForecast.data$dx$male[paste(x), paste(y)]
ex <- dxForecast::dxForecast.data$ex$male
exogen <- ex[paste(y1)]

MM <- c("LC", "FDM", "CoDa", "M6")
B.ex <- doBackTesting(data = D, x = x, 
                      y.fit = y1, y.for = y2, 
                      data.type = "dx", 
                      what = "ex", 
                      exogen = exogen,
                      models = MM)

A <- B.ex$accuracy
A
plot(A)








