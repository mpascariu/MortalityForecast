
# Tue Aug  7 12:55:21 2018 --------- Marius D. Pascariu ---
remove(list = ls())
library(MortalityAccuracy)
library(MortalityLaws)
library(gnm)
library(tidyverse)
# library(StMoMo)

x = 0:95
y1 = 1980:1999
y2 = 2000:2010
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


filter.x = c(0, 60, 70, 80, 90, 100)
O2 <- convertFx(x, D, In = "dx", Out = "ex", lx0 = 1)
O3 <- wide2long(data = O2, x, filter.x)
O3$DATA <- NA
O3[O3$y %in% y1, "DATA"] <- "Training Set"
O3[O3$y %in% y2, "DATA"] <- "Validation Set"


H <- wide.list.2.long.df(data = A$forecasts, x, filter.x)

ggplot(H) + 
  facet_wrap(~x, scales = "free") +
  geom_line(aes(x = y, y = value, color = Name, linetype = Name)) +
  geom_point(data = O3, aes(x = y, y = value, fill = DATA), shape = 21) +
  scale_fill_manual(values = 1:2) +
  guides(colour = guide_legend("Mortality Forecast\nModel"), 
         linetype = guide_legend("Mortality Forecast\nModel"), 
         fill = guide_legend("Demographic Data"))


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

