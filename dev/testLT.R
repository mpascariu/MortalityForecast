
# Tue Aug  7 12:55:21 2018 --------- Marius D. Pascariu ---
remove(list = ls())
library(MortalityAccuracy)
library(MortalityLaws)
library(gnm)
library(tidyverse)
# library(StMoMo)

x = 0:99
y1 = 1980:1999
y2 = 2000:2016
y  = c(y1, y2)
h = max(y2) - max(y1)

D <- dxForecast::dxForecast.data$dx$male[paste(x), paste(y)]
D1 <- dxForecast::dxForecast.data$dx$male[paste(x), paste(y1)]
D2 <- dxForecast::dxForecast.data$dx$male[paste(x), paste(y2)]
ex <- dxForecast::dxForecast.data$ex$male
exogen <- ex[paste(y1)]

M <- doMortalityModels(data = D1, x, y1, data.type = "dx", exogen = exogen)
P <- doForecasts(M, h)
A <- getAccuracy(P, D2, x = 0:100, y2, data.type = "dx", what = "ex")
A

oex <- getObserved(M, what = "ex")
fex <- getFitted(M, what = "ex")
rex <- getResiduals(M, what = "ex")
pex <- getForecasts(P, what = "ex")

# ----------------------------------------------
# BackTesting

bt <- doBackTesting(data = D, 
                    x = x, 
                    y.fit = y1, 
                    y.for = y2, 
                    data.type = "dx", 
                    what = "ex", 
                    exogen = exogen)

plot(bt)

# ----------------------------------------------
# Mon Aug 13 13:35:59 2018 --------- Marius D. Pascariu ---
