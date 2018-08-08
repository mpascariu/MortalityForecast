
# Tue Aug  7 12:55:21 2018 --------- Marius D. Pascariu ---
remove(list = ls())
library(MortalityAccuracy)
library(gnm)
# library(StMoMo)

x = 0:110
y = 1960:2016
dxm <- dxForecast::dxForecast.data$dx$male[paste(x), paste(y)]
dx  <- apply(dxm, 2, function(x) x/sum(x))

M <- getMortalityModels(data = dx, x, y, data.type = "dx")
P <- predict(M, h = 16, ci = 95, jumpchoice = "actual")


fmx <- fitted(M, type = "dx")
fex <- fitted(M, type = "ex")

pex <- predicted(P, type = "ex")



# plot(M2)
# plot(M3)

for (i in 1:length(y)) {
  plot(dx[, i], pch = 16)
  lines(cbind(fit_dx1[,1], fit_dx1)[, i], col = 2)
  lines(fit_dx2[, i], col = 3)
  lines(fit_dx3[, i], col = 4)
  Sys.sleep(1)
}

for (k in 1:h) {
  yr = max(y) + k
  plot(prd_dx1[, k], type = "l", col = 2, main = yr)
  lines(prd_dx2[, k], col = 3)
  lines(prd_dx3[, k], col = 4)
  legend("topleft", legend = c("Mom", "CoDa-LC", "LC"), col = 2:4, lwd = 2)
  Sys.sleep(1)
}






