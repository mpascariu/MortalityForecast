
# Tue Aug  7 12:55:21 2018 --------- Marius D. Pascariu ---
remove(list = ls())
library(MortalityAccuracy)
library(MortalityLaws)
library(gnm)
library(tidyverse)
# library(StMoMo)

x = 0:100
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
y3 <- c(y, y2)
Tdata <- dxForecast::dxForecast.data$dx$male[paste(x), paste(y2)]

# ----------------------------------------------
# BackTesting
what = "ex"
filter.x = c(0, 10, 40, 60, 70, 80, 90, 100)

bt <- doBackTesting(Tdata, P, data.type = "dx", type = what)
O1 <- dxForecast::dxForecast.data$dx$male[paste(x), paste(y3)]
O2 <- convertFx(x, O1, In = "dx", Out = what, lx0 = 1)
O3 <- wide2long(data = O2, x, filter.x)
O3$DATA <- NA
O3[O3$y %in% y, "DATA"] <- "Training Set"
O3[O3$y %in% y2, "DATA"] <- "Validation Set"

head(O3)
tail(O3)

H <- wide.list.2.long.df(data = bt$forecasts, x, filter.x)

ggplot(H) + 
  geom_line(aes(x = y, y = value, color = Name, linetype = Name)) +
  facet_wrap(~x, scales = "free") +
  geom_point(data = O3, aes(x = y, y = value, fill = DATA), color = c(2))


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

