# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Fri Nov 30 11:42:38 2018
# --------------------------------------------------- #
remove(list = ls())


# Code used for downloading and creating of the testing data object in the packages
library(MortalityLaws)
library(tidyverse)


# Download HMD data
source("/Users/mpascariu/Desktop/HMD-credentials.r")
cntr <- c("DNK", "GBRTENW", "NLD", "SWE", "USA")
x <- 0:110
y <- 1950:2016

# Male data
HMD_LTm <- ReadHMD(what = 'LT_m', 
                   countries = cntr,
                   username = username, 
                   password = password, 
                   save = FALSE)

M <- D <- list()
for (k in seq_along(cntr)) {
  K <- cntr[k]
  # dx values
  dxm <- HMD_LTm$data %>% filter(country == K, Year %in% y) %>% 
    select(Year, Age, dx) %>% spread(key = Year, value = dx) %>% select(-Age)
  # mx values
  mxm <- convertFx(x, dxm, from = "dx", to = "mx")
  rownames(dxm) <- rownames(mxm) <- x
  
  D[[k]] <- dxm
  M[[k]] <- mxm
  
}
names(D) <- names(M) <- cntr
HMD_male <- list(dx = D, mx = M, x = x, y = y, countries = cntr)
devtools::use_data(HMD_male, overwrite = TRUE)

# Female data -------------------------------------------------------------
HMD_LTf <- ReadHMD(what = 'LT_f', 
                   countries = cntr,
                   username = username, 
                   password = password, 
                   save = FALSE)

M <- D <- list()
for (k in seq_along(cntr)) {
  K <- cntr[k]
  # dx values
  dxf <- HMD_LTf$data %>% filter(country == K, Year %in% y) %>% 
    select(Year, Age, dx) %>% spread(key = Year, value = dx) %>% select(-Age)
  # mx values
  mxf <- convertFx(x, dxf, from = "dx", to = "mx")
  rownames(dxf) <- rownames(mxf) <- x
  
  D[[k]] <- dxf
  M[[k]] <- mxf
  
}
names(D) <- names(M) <- cntr
HMD_female <- list(dx = D, mx = M, x = x, y = y, countries = cntr)
devtools::use_data(HMD_female, overwrite = TRUE)


