# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Mon Nov 19 15:13:30 2018
# --------------------------------------------------- #
remove(list = ls())
library(MortalityForecast)
library(MortalityLaws)


# Extrapolate old-age human mortality curve using the Kannisto model
extra_mortality <- function(mx,              # Vector or matrix of age specific death-rates
                            x,               # Ages
                            x_fit = 80:90,   # Ages to be considered in calibrarion
                            x_extr = 91:120, # Ages for which to extrapolate the death-rates
                            law = "kannisto"){
  # Fit the mortality model
  M <- MortalityLaw(x = x, mx = mx, 
                    fit.this.x = x_fit, 
                    law = law, 
                    opt.method = "poissonL")
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

# ------------------------------------------


x <- 0:90
y <- 1950:2016
h <- 54
k <- c("GBRTENW","DNK","NLD", "SWE","USA")
J <- "actual" # jumpchoice

mx0 <- HMD_male$mx[k] %>% 
  lapply(as.matrix) %>% 
  lapply(replace.zeros)
# all(sapply(mx0, dim) == dim(mx0[[1]]))

mx3 <- apply(simplify2array(mx0), 1:2, 
             FUN = function(w) prod(w, na.rm = TRUE)^(1/(length(w))))
mx.B <- mx3[paste(x), paste(y)]


mx.k <- mx0[["DNK"]]
mx.k <- mx.k[paste(x), paste(y)]


M1 <- model_LeeCarter(data = mx.k, x = x, y = y)
M2 <- model_LiLee(data = mx.k, data.B = mx.B, x = x, y = y)


J = "fit"
P1 <- predict(M1, h = h, jumpchoice = J)
P2 <- predict(M2, h = h, order = c(0,1,0), include.drift = F, jumpchoice = J)


e0 <- M1$observed.values %>% extra_mortality(mx = ., x = x) %>% 
  convertFx(x = 0:120, data = ., from = "mx", to = "ex")
e1 <- M1$fitted.values %>% extra_mortality(mx = ., x = x) %>% 
  convertFx(x = 0:120, data = ., from = "mx", to = "ex")
e2 <- M2$fitted.values %>% extra_mortality(mx = ., x = x) %>% 
  convertFx(x = 0:120, data = ., from = "mx", to = "ex")
f1 <- P1$predicted.values %>% extra_mortality(mx = ., x = x) %>% 
  convertFx(x = 0:120, data = ., from = "mx", to = "ex")
f2 <- P2$predicted.values %>% extra_mortality(mx = ., x = x) %>% 
  convertFx(x = 0:120, data = ., from = "mx", to = "ex")


f1[1, ]
f2[1, ]

par(mfrow = c(1,1))
plot(y, e0[1,], pch = 16, xlim = range(y, P1$y), ylim = c(69, 88), 
     xlab = "Year", ylab = "Life Expectancy at birth", 
     main = "Without jump-off adjustment")
lines(y, e1[1, ], col = 2)
lines(y, e2[1, ], col = 3)
lines(P1$y, f1[1, ], col = 2, lwd = 2)
lines(P2$y, f2[1, ], col = 3, lwd = 2)
legend("topleft", legend = c("Observed Data", "Lee-Carter", "Li-Lee"),
       pch = c(16, NA, NA), lty = c(NA, 1, 1), col = 1:3, lwd = 2, bty = "n")
















