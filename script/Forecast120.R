# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Sat Dec  1 17:33:19 2018
# --------------------------------------------------- #

# Goal: Forecast mortality for Danish population 120 years in to the future
# Method: Using 6 extrapolative models. Use age 0-95 to fit model. Get old-age
# mortality using the Kannisto model (up to 120).
# How: Using the MortalityForecast R package. To be installed from GitHub:

# devtools::install_github("mpascariu/MortalityForecast")

# ------------------------------------------
remove(list = ls())
library(MortalityForecast)
library(MortalityLaws)
library(tidyverse)

# Input
x  <- 0:95            # Ages
y  <- 1960:2016       # Years
K  <- "DNK"           # Country to be analysed
h  <- 120             # Forecasting horizon

MM <- c("MRWD",      # Multivariate Random Walk w Drift model
        "LeeCarter", # (0,1,0)
        "LiLee",     # (0,1,0) - (0,1,0)
        "HyndmanUllah", 
        "Oeppen",    # (0,1,0)
        "OeppenC")   # (0,1,0) - (0,1,0)

# Mortality data for country: K
mx.country <- HMD_female$mx[[K]][paste(x), paste(y)] %>% replace.zeros
# mx.country <- HMD_male$mx[[K]][paste(x), paste(y)] %>% replace.zeros

# Create a mortality data for a benchmark population
# to be used in "LiLee" and "OeppenC" models
average_mx <- function(w) prod(w, na.rm = TRUE)^(1/(length(w)))

Mx <- HMD_male$mx %>% 
  lapply(as.matrix) %>% 
  lapply(replace.zeros) %>% 
  simplify2array() %>% 
  apply(., 1:2, FUN = average_mx) %>% 
  as.data.frame()

mx.benchmark <- Mx[paste(x), paste(y)]

# Note: We are fitting various mortality model to males data. Most of them 
# are single population model i.e. the estimates are resulted only from data 
# specific to that population. However, the coherent models like "OeppenC" or 
# "LiLee" require additional data. In these cases, here, female data is used.


M <- do.MortalityModels(data = mx.country, 
                        data.B = mx.benchmark,
                        x = x, 
                        y = y, 
                        data.in = "mx", 
                        models = MM)

P <- do.MortalityForecasts(object = M, 
                           h = h, 
                           level = 95, 
                           jumpchoice = "actual")

omx <- get.Observed(M, data.out = "mx")
fmx <- get.Forecasts(P, data.out = "mx") # Forecast m[x]

# ------------------------------------------

# Extrapolate old-age human mortality curve using the Kannisto model
extra_mortality <- function(mx,              # Matrix of death-rates
                            x,               # ages
                            x_fit = 80:95,   # ages to be considered in calibrarion
                            x_extr = 91:120, # ages for which to extrapolate the death-rates
                            law = "kannisto"){
  # Fit the mortality model
  M <- MortalityLaw(x = x, 
                    mx = mx, 
                    fit.this.x = x_fit, 
                    law = law, 
                    opt.method = "LF2")
  pv <- predict(M, x = x_extr)
  
  # which ages are not to be replaced with fitted values?
  L  <- !(x %in% x_extr)  
  
  # Create the output object
  if (is.vector(mx)) {
    names(mx) <- x
    values <- c(mx[L], pv)
    
  } else {
    rownames(mx) <- x
    values <- rbind(mx[L,], pv)
  }
  
  return(values)
}

# Observed m[x] and e[x]
Omx <- extra_mortality(mx = omx, x = x)
Oex <- convertFx(x = 0:120, data = Omx, from = "mx", to = "ex")

# Forecast m[x] and e[x]
Fmx <- lapply(fmx, FUN = function(w) extra_mortality(mx = w, x = x))
Fex <- lapply(Fmx, FUN = function(w) convertFx(x = 0:120, data = w, from = "mx", to = "ex"))
  
# Check forecast e[0] by model
Fe0 <- Fex %>% lapply(as.matrix) %>% simplify2array()
Fe0[1,,]


# To get the forecasts for male popualtion just run again this script using 
# male data in line 35.

















