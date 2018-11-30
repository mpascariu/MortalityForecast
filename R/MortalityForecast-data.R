# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Fri Nov 30 11:48:24 2018
# --------------------------------------------------- #


#' DATA - Demographic data for male population
#'
#' Demographic data for populations between 1950 and 2016 in the following 
#' countries: England \& Wales, Denmark, the Netherlands, Sweden, and USA.
#' The data is provided in the package for testing purposes only.
#' By the time you are using it, the data may be outdated. Download actual 
#' demographic data free of charge from \insertCite{hmd2018;textual}{MortalityForecast}. 
#' Once a username and a password is created on the 
#' \href{https://www.mortality.org}{website} the 
#' \href{https://CRAN.R-project.org/package=MortalityLaws}{MortalityLaws} 
#' R package can be used to extract data directly into your R console.
#' @references  \insertAllCited{}
#' @source \href{https://www.mortality.org}{Human Mortality Database}
"HMD_male"


#' DATA - Demographic data for female population
#' @inherit HMD_male description references
#' @source \href{https://www.mortality.org}{Human Mortality Database}
"HMD_female"
