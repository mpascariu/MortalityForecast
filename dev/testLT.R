
# Tue Aug  7 12:55:21 2018 --------- Marius D. Pascariu ---
remove(list = ls())
library(MortalityAccuracy)
library(gnm)
# library(StMoMo)

x = 0:110
y = 1960:1980
dxm <- dxForecast::dxForecast.data$dx$male[paste(x), paste(y)]

M <- doMortalityModels(data = dxm, x, y, data.type = "dx")
P <- doForecasts(M, h = 16, ci = 95, jumpchoice = "actual")



oex <- getObserved(M, type = "ex")
fex <- getFitted(M, type = "ex")
rex <- getResiduals(M, type = "ex")
pex <- getForecasts(P, type = "ex")


aex <- getAccuracy(P, type = "ex")
aex


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


fit1 <- rwf(EuStockMarkets[1:200,1],h=100)
fit2 <- meanf(EuStockMarkets[1:200,1],h=100)
accuracy(fit1)
accuracy(fit2)
accuracy(fit1,EuStockMarkets[201:300,1])
accuracy(fit2,EuStockMarkets[201:300,1])
plot(fit1)
lines(EuStockMarkets[1:300,1])


