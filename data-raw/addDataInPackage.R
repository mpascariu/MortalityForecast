# --------------------------------------------------- #
# Author: Marius D. Pascariu
# Last update: Wed Nov 14 14:27:03 2018
# --------------------------------------------------- #
remove(list = ls())

# Code used for downloading and creating of the testing data object in the packages
library(MortalityLaws)
library(dplyr)
library(tidyr)


# Download HMD data
source("/Users/mpascariu/Desktop/HMD-credentials.r")
cntr <- c("GBRTENW", "DNK", "NLD", "SWE", "NOR", "USA")
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


