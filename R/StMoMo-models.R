# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Thu Nov 29 13:10:19 2018
# --------------------------------------------------- #


#' Lee-Carter Mortality Model as implemented in StMoMo package
#' @inheritParams doMortalityModels
#' @inheritParams StMoMo::lc
#' @param lx0 lx0
#' @keywords internal
LC <- function(data, 
               x, 
               y, 
               link = "logit", 
               lx0 = 1e5, 
               verbose = FALSE) {
  
  LCfit <- StMoMo::fit(object = lc(link = link), 
                       Dxt = data * lx0, 
                       Ext = data * 0 + lx0, 
                       ages = x, 
                       years = y, 
                       ages.fit = x, 
                       wxt = NULL, 
                       verbose = verbose)
  return(LCfit)
}


#' @keywords internal
x_mean_ages <- function(x, ages) mean(ages) - x


#' Lee-Carter Mortality Model as implemented in StMoMo package
#' @inheritParams LC
#' @keywords internal
PLAT <- function(data, 
                 x, 
                 y, 
                 link = "log", 
                 lx0 = 1e5, 
                 verbose = FALSE) {
  
  Dx  <- data * lx0
  Ex  <- Dx * 0 + lx0
  wxt <- genWeightMat(ages = x, years = y, clip = 3) # weighting matrix
  
  # Model specification
  constPlat <- function(ax, bx, kt, b0x, gc, wxt, ages){
    nYears <- dim(wxt)[2]
    x <- ages
    t <- 1:nYears
    c <- (1 - tail(ages, 1)):(nYears - ages[1])
    xbar   <- mean(x)
    phiReg <- lm(gc ~ 1 + c + I(c ^ 2), na.action = na.omit)
    phi    <- coef(phiReg)
    
     gc <- gc - phi[1] - phi[2] * c - phi[3] * c^2
     kt[2,] <- kt[2,] + 2 * phi[3] * t
     kt[1,] <- kt[1,] + phi[2] * t + phi[3] * (t^2 - 2 * xbar * t)
     ci <- rowMeans(kt, na.rm = TRUE)
     
     ax <- ax + phi[1] - phi[2] * x + phi[3] * x^2
     ax <- ax + ci[1] + ci[2] * (xbar - x)
     kt[1, ] <- kt[1, ] - ci[1]
     kt[2, ] <- kt[2, ] - ci[2]
     
     out <- list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)
     return(out)
  }
  
  M <- StMoMo(link = link, 
              staticAgeFun = TRUE,
              periodAgeFun = c("1", x_mean_ages), 
              cohortAgeFun = "1", 
              constFun = constPlat)
  
  PLATfit <- StMoMo::fit(M, Dxt = Dx, Ext = Ex, ages = x, years = y,
                         ages.fit = x, wxt = wxt, verbose = verbose)
  return(PLATfit)
}







