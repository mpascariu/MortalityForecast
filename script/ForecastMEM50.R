# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Mon Dec  3 15:05:08 2018
# --------------------------------------------------- #

# Goal: Forecast mortality for Danish population 50 years into the future
# Method: Using MEM models. Use age 0-110 to fit model. 
# How: Using the MortalityForecast R package. To be installed from GitHub:

# devtools::install_github("mpascariu/MortalityForecast")

# ------------------------------------------
remove(list = ls())
library(MortalityForecast)
library(MortalityLaws)
library(tidyverse)

# Input
x  <- 0:110          # Ages
y  <- 1960:2016      # Years
K  <- "DNK"          # Country to be analysed
h  <- 50             # Forecasting horizon

MM <- c("MEM4", "MEM5")   

# Mortality data for country: K
# mx.country <- HMD_female$mx[[K]][paste(x), paste(y)] %>% replace.zeros
mx.country <- HMD_male$mx[[K]][paste(x), paste(y)] %>% replace.zeros


M <- do.MortalityModels(data = mx.country, 
                        x = x, 
                        y = y, 
                        data.in = "mx", 
                        models = MM)

P <- do.MortalityForecasts(object = M, 
                           h = h, 
                           level = 95, 
                           jumpchoice = "actual")

omx <- get.Observed(M, data.out = "mx")
oex <- convertFx(x = x, data = omx, from = "mx", to = "ex")

fmx <- get.Forecasts(P, data.out = "mx")
fex <- lapply(fmx, FUN = function(w) convertFx(x = x, data = w, from = "mx", to = "ex"))


# Check forecast e[x] by model
Fe0 <- fex %>% lapply(as.matrix) %>% simplify2array()
Fe0[1,,] # e[0]
Fe0[66,,] # e[65]

# Plots
Fmx <- fmx %>% lapply(as.matrix) %>% simplify2array()
par(mfrow = c(1, 2))
matplot(x = x, y = log(Fmx[,, "MEM4"]), type = "l", main = "MEM4")
matplot(x = x, y = log(Fmx[,, "MEM5"]), type = "l", main = "MEM5")

plot(P$MEM4)
plot(P$MEM5)









