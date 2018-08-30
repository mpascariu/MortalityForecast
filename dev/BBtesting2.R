
# Fri Aug 24 21:50:46 2018 --------- Marius D. Pascariu ---
remove(list = ls())
library(MortalityForecast)


x  <- 0:95
y  <- 1960:2016
dx <- MortalityForecast.data$dx[paste(x), paste(y)]

MM <- c("MRWD", "LeeCarter", "HyndmanUllah", "CoDa", "MEM5", "MEM6")
BB <- doBBackTesting(data = dx, x, y,
                     data.in = "dx", 
                     models = MM,
                     strategy = c(20, 20, 1))




A <- evalAccuracy(BB, data.out = "ex")
A[A$Scenario == "Total", ]

rankA <- doRanking(A)
rankA[rankA$Scenario == "Total", ]



R <- evalRobustness(BB, data.out = "ex")
R[R$Scenario == "Total", ]

rankR <- doRanking(R)
rankR[rankR$Scenario == "Total", ]

# ----------------------------------------------

x = 0:98              # Ages
y1 = 1980:2000        # Training period
y2 = 2001:2016        # Validation period
y  = c(y1, y2)
h = max(y2) - max(y1) # Forecasting horizon

D <- MortalityForecast.data$dx[paste(x), paste(y)] # DATA

# Select various mortality models
MM <- c("MRWD", "LeeCarter", "CoDa", "MEM6")
# Fit & Forecast the models 
B <- doBackTesting(data = D, x = x,
                   y.fit = y1, y.for = y2,
                   data.in = "dx",
                   models = MM)

# Compute accuracy measures
# The measures can be computed for different indicators. Even if it is not 
# impossible to get a diffrent classification and ranking the 
# outcome should be in general the same.
evalAccuracy(B, data.out = "mx")
evalAccuracy(B, data.out = "qx")
evalAccuracy(B, data.out = "dx")

# Rank the model's performance.
A <- evalAccuracy(B, data.out = "ex")
A
R <- doRanking(A)
R

# Visualize the results
plot(B, data.out = "mx", facet = "x")
plot(B, data.out = "mx", facet = "y") 



ex = getForecasts(B$Forecast, "ex")


D2 <- MortalityForecast.data$dx[paste(x), paste(y2)] # DATA

u <- as.matrix(convertFx(x, data = D2, from = "dx", to = "ex"))
u.hat = as.matrix(ex$MEM6)
b <- as.matrix(ex$MRWD)

E  <- u - u.hat # errors
AE <- abs(E)    # absolute errors

f   <- mean(abs(diff(t(u))), na.rm = T)
SE  <- E / f # scaled errors
ASE <- abs(SE)
# 16.Mean Absolute Scaled Error
MASE <- mean(ASE, na.rm = T)
MASE


meanT <- function(z) mean(z, na.rm = TRUE)
f2 = apply(abs(diff(t(u))), 2, meanT)  # compute a scale factor for each time series
SE2 = sweep(E, 1, f2, FUN = "/")
ASE2 = abs(SE2)
MASE2 <- mean(ASE2, na.rm = T)
MASE2

median(ASE, na.rm = T)
median(ASE2, na.rm = T)


