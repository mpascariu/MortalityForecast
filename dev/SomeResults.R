
# Thu Aug 23 15:52:26 2018 --------- Marius D. Pascariu ---
remove(list = ls())
library(MortalityForecast)


x0  <- 0:95
y   <- 1960:2016
dx0 <- MortalityForecast.data$dx[paste(x0), paste(y)]
MM <- c("MRWD", "LeeCarter", "HyndmanUllah", "CoDa", "MEM4", "MEM5", "MEM6")

BB0 <- doBBackTesting(data = dx0, x0, y,
                     data.in = "dx",
                     data.out = "ex",
                     models = MM,
                     strategy = c(20, 20, 1))
BB0$accuracy

x65  <- 65:95
dx65 <- MortalityForecast.data$dx[paste(x65), paste(y)]
BB65 <- doBBackTesting(data = dx65, x65, y,
                     data.in = "dx",
                     data.out = "ex",
                     models = MM,
                     strategy = c(20, 20, 1))
BB65$accuracy



BB65$scenarios[[1, 3]]

plot(BB65$results$S1, facet = "y")
plot(BB0$results$S1, facet = "x")
