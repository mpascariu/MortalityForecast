# Marius D. Pascariu --- Wed Sep 26 13:43:34 2018 ------------------------------
remove(list = ls())
library(MortalityForecast)
library(MortalityLaws)
library(tidyverse)


# Download HMD data
cntr <- c("AUS","CAN","FRATNP","GBRTENW","ITA","NLD","ESP","CHE","SWE","USA", "DNK")
cntrN <- c("Australia", "Canada", "France", "England and Wales", "Italy", 
           "The Netherlands", "Spain", "Switzerland", "Sweden", "USA", "Denmark")
# cntr = "DNK"
# user <- 'zpascariu@outlook.com'
# pwd  <- '1513208587'
# ReadHMD(what = 'LT_f', countries = cntr, username = user, password = pwd, save = TRUE)

load("HMD_LT_f.Rdata")
HMD <- HMD_LT_f



x  <- 0:95
y  <- 1950:2016
s2 <- 5 # forecasting horizon
s3 <- 1  # step
n_fit_max <- min(30, diff(range(y)) - s2 + 1)
nfit <- 3:n_fit_max # no. of years to be used in fitting 
input_data <- "mx"  # indicator used as input in fitting
output_data <- "ex" # indicator used in accuracy evaluation
M <- c("MRWD", "LeeCarter") # Tested models
accm <- c("ME", "MAE", "MAPE", "sMAPE", "MASE") # Considered accuracy measures 

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

mx[mx == 0] <- 1e-10

BB <- list()
A <- NULL
for (i in seq_along(nfit)) {
  s1 <- nfit[i]
  cat("\n\nnfit = ", s1)
  BB[[i]] <- doBBackTesting(data = mx, x, y, data.in = input_data, models = M,
                            strategy = c(s1, s2, 1))
  A <- rbind(A, evalAccuracy(BB[[i]], data.out = output_data, measures = accm)[seq_along(M), ])
}


A1 <- A[A$Model == "MRWD", ]
A2 <- A[A$Model == "LeeCarter", ]


R1 <- data.frame(nfit ,doRanking(A1))
R2 <- data.frame(nfit ,doRanking(A2))

R1
R2


plot_results <- function(A, nfit) {
  R <- data.frame(nfit, doRanking(A))
  par(mfrow = c(1, 1), cex.axis = 2)
  plot(nfit, A$ME, type = "l", lwd = 2, main = "ME")
  # abline(v = R[A$ME == min(abs(A$ME)), "nfit"], lty = 2, col = 2, lwd = 2)
  par(new = TRUE)
  plot(nfit, A$MAE, type = "l", lwd = 2, main = "MAE")
  # plot(nfit, A$MAPE, type = "l", lwd = 2, main = "MAPE")
  # plot(nfit, A$sMAPE, type = "l", lwd = 2, main = "sMAPE")
  # plot(nfit, A$MASE, type = "l", lwd = 2, main = "MASE")
}

plot_results(A1[-c(1:4),], nfit[-c(1:4)])
plot_results(A2[-c(1:4),], nfit[-c(1:4)])







