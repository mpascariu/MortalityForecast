# Tue Aug 14 21:04:37 2018 --------- Marius D. Pascariu ---
remove(list = ls())
library(MortalityForecast)

# Input data
x  <- 0:95
y  <- 1960:2016
dx <- MortalityForecast.data$dx[paste(x), paste(y)]

BB <- doBBackTesting(data = dx, x, y,
                     data.in = "dx", 
                     data.out = "ex", 
                     models = c("MRWD", "FDM", "CoDa", "MEM5", "MEM6"),
                     strategy = c(20, 20, 2),
                     xa = 1:94)

A <- BB$accuracy
A
doRanking(A)

BB$results$S1$MortalityModels$


# ----------------------------------------------

xf = c(0, 25, 45, 65, 80, 90)
plot(B[[1]], facet = "x", which = xf)
plot(B[[2]], facet = "x", which = xf)
plot(B[[3]], facet = "x", which = xf)
plot(B[[4]], facet = "x", which = xf)

B[[1]]$accuracy
B[[2]]$accuracy
B[[3]]$accuracy
B[[4]]$accuracy
# ----------------------------------------------




