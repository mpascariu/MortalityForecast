# Mon Aug 20 17:36:19 2018 --------- Marius D. Pascariu ---
remove(list = ls())
library(MortalityForecast)


x = 0:98
y1 = 1980:1999
y1 = 1990:2005
y2 = 2006:2016
y  = c(y1, y2)
h = max(y2) - max(y1)

D <- dxForecast::dxForecast.data$dx$male[paste(x), paste(y)]
ex <- dxForecast::dxForecast.data$mode$male
exogen <- ex[paste(y1)]

MM <- c("LC", "FDM", "CoDa", "M5", "M5X", "M6")
B.ex <- doBackTesting(data = D, x = x, 
                      y.fit = y1, y.for = y2, 
                      data.type = "dx", 
                      what = "ex", 
                      exogen = exogen,
                      models = MM)

B.mx <- doBackTesting(data = D, x = x, 
                      y.fit = y1, y.for = y2, 
                      data.type = "dx", 
                      what = "mx", 
                      exogen = exogen,
                      models = MM)

B.dx <- doBackTesting(data = D, x = x, 
                      y.fit = y1, y.for = y2, 
                      data.type = "dx", 
                      what = "dx", 
                      exogen = exogen,
                      models = MM)

B.ex$accuracy
B.mx$accuracy
B.dx$accuracy

plot(B.ex, facet = "x")
plot(B.ex, facet = "y")

plot(B.mx, facet = "x")
plot(B.mx, facet = "y")

plot(B.dx, facet = "x")
plot(B.dx, facet = "y")
