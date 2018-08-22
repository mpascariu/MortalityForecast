# Tue Aug 14 21:04:37 2018 --------- Marius D. Pascariu ---
remove(list = ls())
library(MortalityForecast)

# Input data
x  = 0:95
xa = 0:95
y = 1960:2016
strategy = c(f = 20, h = 20, s = 4)
MM <- c("MRWD", "LC", "FDM", "CoDa", "M5", "M6")

D = MortalityForecast.data$dx[paste(x), paste(y)]

buildScenarios <- function(y, S) {
  # Build scenarios -  method 1
  bop_fit = seq(from = max(y) - S[1] - S[2] + 1, to = min(y), by = -1 * S[3])
  eop_fit = bop_fit + S[1] - 1
  bop_fc  = eop_fit + 1
  eop_fc  = bop_fc + S[2] - 1
  
  out <- tibble(
    scenario = paste0("S", 1:length(bop_fit)),
    fit = mapply(":", bop_fit, eop_fit, SIMPLIFY = F),
    forecast = mapply(":", bop_fc, eop_fc, SIMPLIFY = F)
  )
  
  return(out)
}

doBackTESTING <- function(data, x, y, 
                          in.data = c("qx", "mx", "dx", "lx"),
                          out.data = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"),
                          models = c("MRWD", "LC"), 
                          measures = c("ME", "MAE", "MAPE", "sMAPE", "MRAE", "MASE"), 
                          xA = NULL, yA = NULL,
                          jumpchoice = c("actual", "fit"),
                          level = 95,
                          verbose = FALSE, ...) {
  
  in.data <- match.arg(in.data)
  out.data <- match.arg(out.data)
  jumpchoice <- match.arg(jumpchoice)
  
}

S <- buildScenarios(y, strategy)


# Do Back-testing
nc = nrow(S) # no. of cases
B <- list()
A <- 0
for (k in 1:nc) {
  yf = S[[k, "fit"]]
  yh = S[[k, "forecast"]]
  y_ = c(yf, yh)
  Dk <- D[, paste(y_)]
  cat(paste0("\nTest case ", k, "/", nc, ": "))
  Bk <- doBackTesting(data = Dk, x = x, xa = xa, 
                      y.fit = yf, y.for = yh, models = MM,
                      data.type = "dx", what = "mx")
  Ak <- Bk$accuracy$results
  
  B[[k]] <- Bk
  A <- A + Ak
}

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








y  <- 1900
x  <- as.numeric(rownames(ahmd$mx))
Dx <- ahmd$Dx[, paste(y)]
Ex <- ahmd$Ex[, paste(y)]

LT1 <- LifeTable(x, Dx = Dx, Ex = Ex)















