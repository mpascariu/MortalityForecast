
# Wed Aug 22 16:58:09 2018 --------- Marius D. Pascariu ---
remove(list = ls())
library(MortalityForecast)

x = 0:98              # Ages
y1 = 1965:1990        # Training period
y2 = 1991:2016        # Validation period
y  = c(y1, y2)
h = max(y2) - max(y1) # Forecasting horizon

D <- MortalityForecast.data$dx[paste(x), paste(y)] # DATA

MM <- c("MRWD", "LC", "CoDa", "M6")
B.ex <- doBackTesting(data = D, x = x, 
                      y.fit = y1, y.for = y2, 
                      data.type = "dx", 
                      what = "ex", 
                      models = MM)
B.mx <- doBackTesting(data = D, x = x, 
                      y.fit = y1, y.for = y2, 
                      data.type = "dx", 
                      what = "mx", 
                      models = MM)
B.ex$accuracy
B.mx$accuracy

plot(B.ex, facet = "x")
plot(B.qx, facet = "x")
plot(B.ex, facet = "y")  
