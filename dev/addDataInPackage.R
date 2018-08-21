# Tue Aug 21 12:01:50 2018 --------- Marius D. Pascariu ---
# Code used for downloading and creating of the testing data object in the packages
remove(list = ls())
library(MortalityLaws)
library(dplyr)
library(tidyr)

# Download HMD data
HMDusr = 'zpascariu@outlook.com'
psd = '1513208587'
cntr = 'GBRTENW'


# Dx
HMD_LTm <- ReadHMD(what = 'LT_m', countries = cntr,
                  username = HMDusr, password = psd, save = F)

x <- 0:110
yr <- 1960:2016
LTm <- HMD_LTm$data %>% filter(country == cntr, Year %in% yr) %>% 
  select(Year, Age, dx, ex)

# ----------------------------------------------
# dx values
dxm <- LTm %>% select(-ex) %>% spread(key = Year, value = dx) %>% select(-Age)
rownames(dxm) <- x
# ----------------------------------------------
# mx values
mxm <- convertFx(x, dxm, In = "dx", Out = "mx")
# ----------------------------------------------
# ex values
exm <- LTm %>% filter(Age == 0) %>% select(ex) %>% c
names(exm$ex) <- yr
# ----------------------------------------------
# mode values
findMode <- function(D, years) {
  mode <- c()
  for (i in 1:length(years)) {
    Di <- D[D$Year == years[i] & D$Age >= 30, ]
    mode[i] <- Di[Di$dx == max(Di$dx), 'Age'][1]
  }
  names(mode) <- years
  return(mode)
}
modem <- findMode(LTm, years = yr)
# ----------------------------------------------


MortalityForecast.data <- list(dx = dxm, mx = mxm, ex = exm$ex,mode = modem,
                               x = x, y = yr, country = cntr)

devtools::use_data(MortalityForecast.data, overwrite = TRUE)


