# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Sat Dec  1 11:19:52 2018
# --------------------------------------------------- #
remove(list = ls())


# Code used for downloading and creating of the testing data object in the packages
library(MortalityLaws)
library(tidyverse)
library(devtools)
source("/Users/mpascariu/Desktop/HMD-credentials.r")

cntr <- c("DNK", "GBRTENW", "NLD", "SWE", "USA")

format_data <- function(data, 
                        x = 0:110, 
                        y = 1950:2016) {
  
  M <- D <- list()
  for (k in seq_along(cntr)) {
    K <- cntr[k]
    # dx values
    dxm <- data %>% 
      filter(country == K, Year %in% y) %>% 
      select(Year, Age, dx) %>% 
      spread(key = Year, value = dx) %>% 
      select(-Age) 
    
    # mx values
    mxm <- convertFx(x, dxm, from = "dx", to = "mx")
    
    rownames(dxm) <- rownames(mxm) <- x
    
    D[[k]] <- dxm
    M[[k]] <- mxm
  }
  
  names(D) <- names(M) <- cntr
  
  out <- list(dx = D, 
              mx = M, 
              x = x, 
              y = y, 
              countries = cntr)
  return(out)
}

# Male data -------------------------------------------------------------
HMD_LTm <- ReadHMD(what = 'LT_m', 
                   countries = cntr,
                   username = username, 
                   password = password, 
                   save = FALSE)

HMD_male <- format_data(data = HMD_LTm$data)
use_data(HMD_male, overwrite = TRUE)

# Female data -------------------------------------------------------------
HMD_LTf <- ReadHMD(what = 'LT_f', 
                   countries = cntr,
                   username = username, 
                   password = password, 
                   save = FALSE)

HMD_female <- format_data(data = HMD_LTf$data)
use_data(HMD_female, overwrite = TRUE)


