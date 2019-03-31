# --------------------------------------------------- #
# Author: Marius D. PASCARIU
# License: GNU General Public License v3.0
# Last update: Tue Mar 19 15:29:37 2019
# --------------------------------------------------- #
remove(list = ls())
set.seed(1234)
library(MortalityForecast)
library(MortalityLaws)
library(tidyverse)


# Download HMD data
cntr <- c("AUS","CAN","FRATNP","GBRTENW","ITA","NLD","ESP","CHE","SWE","USA")
cntrN <- c("Australia", "Canada", "France", "England and Wales", "Italy",
           "The Netherlands", "Spain", "Switzerland", "Sweden", "USA")

load("HMD_LT_m.Rdata")
HMD <- HMD_LT_m

doData <- function(H, x, y, pop, type = "dx") {
  dxm <- H$data %>%
    filter(country %in% pop, Year %in% y, Age %in% x) %>%
    select(Year, Age, type) %>%
    spread(key = Year, value = type) %>%
    select(-Age) %>%
    apply(., MARGIN = 2, FUN = as.numeric) %>%
    as.data.frame() %>%
    replace.zeros()
  rownames(dxm) <- x
  return(dxm)
}

input_data <- "dx"
output_data <- "ex"

# ----------------------------------------------
# 1. Generate data for Figure1
x1 = 0:110              # ages
y1 = 1988:1990          # years
n1 = 2:7
k1 = "USA"              # country

D = doData(HMD, x = x1, y = y1, pop = k1)
D = convertFx(x1, D, input_data, "dx")
D[112:121, ] <- 0

for (n in n1) {
  assign(x = paste0('C', n, '_', k1), model.MEM(data = D, x = 0:120, n = n))
  cat(".")
}


# ----------------------------------------------
# 2. Fit & forecast:
x2 = 0:110
y2 = 1980:2016
k2 = "GBRTENW"                   # country
n2 = 6

D = doData(HMD, x = x2, y = y2, pop = k2, type = input_data)
D = convertFx(x2, D, input_data, "dx")

M6_GBRTENW = model.MEM(data = D, x = x2, n = n2)
P6_GBRTENW = predict(M6_GBRTENW, h = 24)

y2_ = 1960:2016
D = doData(HMD, x = x2, y = y2_, pop = k2)
D = convertFx(x2, D, input_data, "dx")
M6_GBRTENW_0 = model.MEM(data = D, x = x2, n = 3)


# ----------------------------------------------
# ----------------------------------------------
# 3. Back-Testing Age: 0 - 95 - All countries

x <- 0:95
y <- 1960:2016
S <- c(f = 20, h = 20, s = 1)
M <- c("MRWD", "LeeCarter", "HyndmanUllah", "RenshawHaberman", "Oeppen", "MEM6")

B = A = R <- list()
for (i in seq_along(cntr)) {
  Mi <- M
  K  <- cntr[i]
  D  <- doData(HMD, x = x, y = y, pop = K)
  if (K %in% c("NLD", "SWE", "ITA", "ESP", "CHE")) {
    Mi[Mi == "MEM6"] <- "MEM5"
  }
  cat("\n", i, ".", K)

  B[[i]] <- do.BBackTesting(data = D,
                            x = x,
                            y = as.numeric(colnames(D)),
                            data.in = input_data,
                            models = Mi,
                            strategy = S)

  A[[i]] <- evalAccuracy(B[[i]], output_data)
  R[[i]] <- do.Ranking(A[[i]])
}

names(B) = names(A) = names(R) <- cntr


# ----------------------------------------------
# Back-Testing MAxEntMortality only + MRWD

MEM <- c("MRWD", "MEM2", "MEM3", "MEM4", "MEM5", "MEM6")
K  <- "GBRTENW"
D  <- doData(HMD, x = x, y = y, pop = K)
B2 <- do.BBackTesting(data = D,
                     x = x,
                     y = y,
                     data.in = input_data,
                     models = MEM,
                     strategy = S)

A2 <- evalAccuracy(B2, output_data)
R2 <- do.Ranking(A2)


remove(D, HMD, HMD_LT_m, Mi, MEM, n, y1, y2, y2_, x, x1, x2, i,
       input_data, output_data, K, k1, k2, S, doData)
# ----------------------------------------------
# Save Results
save.image(file = "dxForecastResults.Rdata")



