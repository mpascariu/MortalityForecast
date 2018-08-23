
# Thu Aug 23 15:52:26 2018 --------- Marius D. Pascariu ---
remove(list = ls())
library(MortalityForecast)
library(MortalityLaws)
library(dplyr)
library(tidyr)


# Download demographic data for 3 countries in 1x1 format 
age_int  <- 1  # age interval: 1,5
year_int <- 1  # year interval: 1,5,10
interval <- paste0(age_int, "x", year_int)  # --> 1x1
# And the 3 countries: Sweden Denmark and USA. We have to use the HMD codes
cntr  <- c('SWE', 'DNK', 'USA')  

# Download death counts. We don't want to export data outside R.
LT <- ReadHMD(what = "LT_m",
              countries = cntr,
              interval  = "1x1",
              username  = "zpascariu@outlook.com",
              password  = "1513208587",
              save = FALSE)

formatData <- function(hmd, Country, ages, years, data.out = "mx") {
  out <- hmd$data %>% 
    filter(country == Country, Year %in% years, Age %in% ages) %>% 
    select("Year", "Age", data.out) %>% 
    spread(key = Year, value = data.out) %>% select(-Age)
  rownames(out) <- ages
  out
}

# ----------------------------------------------
x <- 0:95
y <- 1960:2016
D1 <- formatData(LT, "USA", ages = x, years = y, "mx")
MM <- c("MRWD", "LeeCarter", "HyndmanUllah", "CoDa", "MEM4", "MEM5", "MEM6")

BB0_usa <- doBBackTesting(data = D1, x, y,
                          data.in = "mx",
                          data.out = "ex",
                          models = MM,
                          strategy = c(20, 20, 1))
BB0_usa$accuracy
# ----------------------------------------------
# ----------------------------------------------
x <- 0:95
y <- 1960:2016
D2 <- formatData(LT, "DNK", ages = x, years = y, "mx")
MM <- c("MRWD", "LeeCarter", "HyndmanUllah", "CoDa", "MEM4", "MEM5")

BB0_dnk <- doBBackTesting(data = D2, x, y,
                          data.in = "mx",
                          data.out = "ex",
                          models = MM,
                          strategy = c(20, 20, 1))
BB0_dnk$accuracy
# ----------------------------------------------

x0  <- 0:95
y   <- 1960:2016
dx0 <- MortalityForecast.data$dx[paste(x0), paste(y)]

BB0 <- doBBackTesting(data = dx0, x0, y,
                     data.in = "dx",
                     data.out = "ex",
                     models = MM,
                     strategy = c(20, 20, 1))
BB0$accuracy

x65  <- 65:95
dx65 <- MortalityForecast.data$dx[paste(x65), paste(y)]
BB65 <- doBBackTesting(data = dx65, x65, y,
                     data.in = "dx",
                     data.out = "ex",
                     models = MM,
                     strategy = c(20, 20, 1))
BB65$accuracy



BB65$scenarios[[1, 3]]

plot(BB65$results$S1, facet = "y")
plot(BB0$results$S1, facet = "x")
