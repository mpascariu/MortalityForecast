# Tue Aug 14 21:04:37 2018 --------- Marius D. Pascariu ---
remove(list = ls())
library(MortalityForecast)
library(tibble)

# Input data
x  <- 0:95
y  <- 1970:2016
dx <- MortalityForecast.data$dx[paste(x), paste(y)]

BB <- doBBackTesting(data = dx, x, y,
                     data.in = "dx", 
                     data.out = "ex", 
                     models = c("MRWD", "LC"),
                     strategy = c(20, 20, 5))

ls(BB)

BB$scenarios[[2, 1]]
BB$scenarios[[2, 2]]
BB$scenarios[[2, 3]]

BB$
A_ = round(A/nc, 4)
A_
doRanking(A_)

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




