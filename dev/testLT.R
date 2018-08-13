
# Tue Aug  7 12:55:21 2018 --------- Marius D. Pascariu ---
remove(list = ls())
library(MortalityAccuracy)
library(MortalityLaws)
library(gnm)
library(tidyverse)
# library(StMoMo)

x = 0:98
y1 = 1980:1999
y2 = 2000:2016
y  = c(y1, y2)
h = max(y2) - max(y1)
models = c("LC", "CoDa", "M6")

D <- dxForecast::dxForecast.data$dx$male[paste(x), paste(y)]
D1 <- dxForecast::dxForecast.data$dx$male[paste(x), paste(y1)]
D2 <- dxForecast::dxForecast.data$dx$male[paste(x), paste(y2)]
ex <- dxForecast::dxForecast.data$ex$male
exogen <- ex[paste(y1)]

M <- doMortalityModels(data = D1, x, y1, data.type = "dx", exogen = exogen,
                       models)
ls(M)

P <- doForecasts(M, h)
ls(P)
A <- getAccuracy(P, D2, x = 0:100, y2, data.type = "dx", what = "ex")
A

oex <- getObserved(M, what = "ex")
fex <- getFitted(M, what = "ex")
rex <- getResiduals(M, what = "ex")
pex <- getForecasts(P, what = "ex")



# ----------------------------------------------
# BackTesting

B.ex <- doBackTesting(data = D, x = x, 
                      y.fit = y1, y.for = y2, 
                      data.type = "dx", 
                      what = "ex", 
                      exogen = exogen,
                      models = models)

B.qx <- doBackTesting(data = D, x = x, 
                      y.fit = y1, y.for = y2, 
                      data.type = "dx", 
                      what = "qx", 
                      exogen = exogen,
                      models = models)

B.mx <- doBackTesting(data = D, x = x, 
                      y.fit = y1, y.for = y2, 
                      data.type = "dx", 
                      what = "mx", 
                      exogen = exogen,
                      models = models)

B.lx <- doBackTesting(data = D, x = x, 
                      y.fit = y1, y.for = y2, 
                      data.type = "dx", 
                      what = "lx", 
                      exogen = exogen,
                      models = models)


x2 <- c(0, 25, 65, 75, 85, 95)
wi = 12
he = 6
pdf("Forecast_mx.pdf", width = wi, height = he)
plot(B.mx, x2); dev.off()
pdf("Forecast_qx.pdf", width = wi, height = he)
plot(B.qx, x2); dev.off()
pdf("Forecast_ex.pdf", width = wi, height = he)
plot(B.ex, x2); dev.off()
pdf("Forecast_lx.pdf", width = wi, height = he)
plot(B.lx, x2); dev.off()

# ----------------------------------------------
# Mon Aug 13 13:35:59 2018 --------- Marius D. Pascariu ---
