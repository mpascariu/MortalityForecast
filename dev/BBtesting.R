# Tue Aug 14 21:04:37 2018 --------- Marius D. Pascariu ---
remove(list = ls())
library(MortalityAccuracy)

# Input data
x  = 0:98
xa = 0:98
y = 1980:2016
n = 20
h = 20
step = 3
MM <- c("LC", "M5X", "M6X", "M5", "M6")
D  = dxForecast::dxForecast.data$dx$female[paste(x), paste(y)]
ex <- dxForecast::dxForecast.data$ex$female
exogen <- ex[paste(y)]

# Build scenarios -  method 1
# bop_fit = seq(from = max(y) - (h + n) + 1, to = min(y), by = -1 * step)
# eop_fit = bop_fit + n - 1
# bop_fc = eop_fit + 1
# eop_fc = bop_fc + h - 1
# S = data.frame(bop_fit, eop_fit, bop_fc, eop_fc)
# S

# Build scenarios -  method 2
bop_fit = min(y)
eop_fc = max(y)
eop_fit = seq(bop_fit + n - 1, eop_fc - 6, by = step)
bop_fc = eop_fit + 1
S = data.frame(bop_fit, eop_fit, bop_fc, eop_fc)
S

# Do Back-testing

nc = nrow(S) # no. of cases
B <- list()
A <- 0
for (k in 1:nc) {
  y_fit_k = S[k, 1]:S[k, 2]
  y_for_k = S[k, 3]:S[k, 4]
  y_k = S[k, 1]:S[k,4]
  Dk <- D[, paste(y_k)]
  cat(paste0("\nTest case ", k, "/", nc, ": "))
  Bk <- doBackTesting(data = Dk, x = x, xa = xa, y_fit_k, y_for_k, 
                      data.type = "dx", what = "dx", exogen = exogen,
                      models = MM)
  Ak <- Bk$accuracy$results
  
  B[[k]] <- Bk
  A <- A + Ak
}

A_ = round(A/nc, 4)
A_

doRanking(A_)

# ----------------------------------------------

xf = c(0, 1, 25, 65, 80, 98)
plot(B[[1]], xf)
plot(B[[2]], xf)
plot(B[[3]], xf)
plot(B[[4]], xf)

B[[1]]$accuracy
B[[2]]$accuracy
B[[3]]$accuracy
B[[4]]$accuracy
