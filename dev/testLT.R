
# Tue Aug  7 12:55:21 2018 --------- Marius D. Pascariu ---
remove(list = ls())
library(MortalityAccuracy)
library(MortalityLaws)
library(gnm)
# library(StMoMo)

x = 0:110
y = 1980:1999
h = 17
dxm <- dxForecast::dxForecast.data$dx$male[paste(x), paste(y)]
ex <- dxForecast::dxForecast.data$ex$male
exogen <- ex[paste(y)]

M <- doMortalityModels(data = dxm, x, y, data.type = "dx", exogen = exogen)
P <- doForecasts(M, h, ci = 95, jumpchoice = "actual")


oex <- getObserved(M, type = "ex")
fex <- getFitted(M, type = "ex")
rex <- getResiduals(M, type = "ex")
pex <- getForecasts(P, type = "ex")


y2 <- max(y) + 1:h
Tdata <- dxForecast::dxForecast.data$dx$male[paste(x), paste(y2)]
doBackTesting(Tdata, P, data.type = "dx", type = "ex")
doBackTesting(Tdata, P, data.type = "dx", type = "mx")
doBackTesting(Tdata, P, data.type = "dx", type = "qx")


# plot(M2)
# plot(M3)

# for (i in 1:length(y)) {
#   plot(dx[, i], pch = 16)
#   lines(cbind(fit_dx1[,1], fit_dx1)[, i], col = 2)
#   lines(fit_dx2[, i], col = 3)
#   lines(fit_dx3[, i], col = 4)
#   Sys.sleep(1)
# }
# 
# for (k in 1:h) {
#   yr = max(y) + k
#   plot(prd_dx1[, k], type = "l", col = 2, main = yr)
#   lines(prd_dx2[, k], col = 3)
#   lines(prd_dx3[, k], col = 4)
#   legend("topleft", legend = c("Mom", "CoDa-LC", "LC"), col = 2:4, lwd = 2)
#   Sys.sleep(1)
# }

