# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Tue Nov 27 15:36:32 2018
# --------------------------------------------------- #
remove(list = ls())
library(MortalityForecast)
library(MortalityLaws)
library(tidyverse)


# Download HMD data
cntr <- c("AUS","CAN","FRATNP","GBRTENW","ITA","NLD","ESP","CHE","SWE","USA", "DNK")
cntrN <- c("Australia", "Canada", "France", "England and Wales", "Italy", 
           "The Netherlands", "Spain", "Switzerland", "Sweden", "USA", "Denmark")

load("dev/HMD_LT_f.Rdata")
HMD <- HMD_LT_f

x  <- 0:100
y  <- 1900:2016
s2 <- 20 # forecasting horizon
s3 <- 3  # step
n_fit_max <- min(60, diff(range(y)) - s2 + 1)
nfit <- 5:n_fit_max # no. of years to be used in fitting 
input_data <- "mx"  # indicator used as input in fitting
output_data <- "ex" # indicator used in accuracy evaluation
M <- c("MRWD", "LeeCarter", "HyndmanUllah", "Oeppen") # Tested models

doData <- function(H, x, y, pop, type = "mx") {
  dxm <- H$data %>% 
    filter(country %in% pop, Year %in% y, Age %in% x) %>% 
    select(Year, Age, type) %>% 
    spread(key = Year, value = type) %>% 
    select(-Age) %>% 
    apply(., MARGIN = 2, FUN = as.numeric) %>% 
    as.data.frame()
  rownames(dxm) <- x
  return(dxm)
}

mx <- doData(HMD, x = x, y = y, pop = "DNK", type = input_data)
mx <- replace.zeros(mx)

BB <- list()
for (i in seq_along(nfit)) {
  s1 <- nfit[i]
  cat("\n\nnfit = ", s1)
  BB[[i]] <- doBBackTesting(data = mx, x, y, data.in = input_data, models = M,
                            strategy = c(s1, s2, s3))
}

A <- NULL
n <- seq_along(M)
for (j in seq_along(nfit)) {
  A <- rbind(A, evalAccuracy(BB[[j]], data.out = output_data)[n, ])
  cat(".")
}

A1 <- A[A$Model == "MRWD", ]
A2 <- A[A$Model == "LeeCarter", ]
A3 <- A[A$Model == "HyndmanUllah", ]
A4 <- A[A$Model == "Oeppen", ]

# print(tbl_df(A), n=48)

plot_results <- function(A, nfit) {
  par(mfrow = c(3, 6), cex.axis = 1)
  for (j in 4:19) {
    cn <- colnames(A)
    plot(nfit, A[[j]], type = "l", main = cn[j], xlim = range(nfit),
         col = 2, axes = F,
         xlab = "Years in fit", ylab = "Accuracy measure")
    axis(1, at = trunc(quantile(nfit)))
    axis(2)
  }
}

# quartz()
plot_results(A1, nfit)
plot_results(A2, nfit)
plot_results(A3, nfit)
plot_results(A4, nfit)



plot(BB[[7]]$results$S10, data.out = "mx", facet = "x")




