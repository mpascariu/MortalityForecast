# Tue Aug 28 11:48:25 2018 --------- Marius D. Pascariu ---
remove(list = ls())
set.seed(1234)
library(MortalityForecast)
library(tidyverse)

country_tab <- data.frame(code = c("GBRTENW"), name = c("England \\& Wales"))

dxHMD <- MortalityForecast.data$dx
x.fit <- as.numeric(rownames(dxHMD))
dxHMD <- apply(dxHMD, 2, as.numeric) %>% as.data.frame()
rownames(dxHMD) <- x.fit

# ----------------------------------------------
# Generate data fof Figure1
k <- "GBRTENW" # country
Y1 <- 1969:1991 # years
R1 <- 2:8       # ranks

dxk <- dxHMD %>% select(paste(Y1))
for (n in R1) {
  assign(x = paste0('C', n, '_', k), 
         value = fitMaxEntMortality(data = dxk, n = n, exogen = exogen))
  cat(".")
}


theme1 <- function(){
  theme_bw() +
    theme(panel.background = element_rect(fill = 'grey98'),
          legend.position = "top",
          legend.background = element_rect(fill = NA),
          legend.key.size   = unit(0.8, "cm"),
          legend.text.align = 0,
          strip.text = element_text(colour = 1, size = 14),
          legend.text  = element_text(colour = 1, size = 13),
          legend.title = element_text(colour = 1, size = 13),
          axis.text    = element_text(colour = 1, size = 13),
          axis.title   = element_text(colour = 1, size = 14, face = 'bold'),
          axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
          axis.line.x  = element_line(linetype = 1),
          axis.ticks   = element_line(size = 1),
          panel.spacing.x = unit(.5, "cm"),
          panel.spacing.y = unit(.2, "cm"),
          plot.margin   = unit(c(0.1, 1, 0.1, 0.1), "cm")
    )
}

plot_coverage <- function(cntr, yr){
  
  models <- paste0('C', R1,'_', cntr)
  dens   <- paste0('Density (N = ', R1, ')')
  
  DX <- NULL
  for (k in 1:length(models)) {
    Mk <- get(models[k])
    dx.obs <- as.data.frame(Mk$observed.values) %>% mutate(age = Mk$x) %>% 
      gather(., year, fx, -age) %>% filter(year == yr) %>% 
      mutate(country = cntr, type = 'Observed', rank = dens[k])
    dx.hat <- as.data.frame(Mk$fitted.values) %>% mutate(age = Mk$x) %>% 
      gather(., year, fx, -age) %>% filter(year == yr) %>% 
      mutate(country = cntr, type = 'Estimated', rank = dens[k])
    DX <- rbind(DX, dx.obs, dx.hat)
  }
  
  # compute coverage
  DX1 <- DX[DX$type == "Estimated" & DX$age %in% 0:110, ]
  DX2 <- DX[DX$type == "Observed" & DX$age %in% 0:110, ]
  DX1$area <- pmin(DX1$fx, DX2$fx)
  coverage <- DX1 %>% group_by(rank) %>% summarise(covr = round(sum(area)*100, 2)) 
  
  dx.obs$rank <- 'Density'
  DX          <- rbind(DX, dx.obs)
  # ggplot
  ggplot(DX, aes(x = age, y = fx)) +
    geom_ribbon(aes(ymin = 0, ymax = fx, fill = type, alpha = type, color = type)) +
    facet_wrap(~ rank, nrow = 2) +
    scale_alpha_manual(name = "Legend:", values = c(0.12, 0.5)) + 
    scale_fill_manual(name = "Legend:", values = c(4, 1)) +
    scale_color_manual(name = "Legend:", values = c(4, NA)) +
    scale_x_continuous(breaks = c(0, 35, 70, 110)) +
    xlab("Age, x") + ylab(expression(f(x))) + 
    geom_text(data = coverage, inherit.aes = FALSE, size = 5, hjust = "inward",
              aes(x = 5, y = .027, label = paste0('Coverage:\n  ', covr,'%') )) +
    theme1()
}

# Figure 1  
plot_coverage("GBRTENW", 1990) 
  
# ----------------------------------------------
# Fit & forecast: all countries
# 1 x 2 x 6
x = 0:110
Y2 <- 1980:2016
S2 <- 24
R2 <- 6

dxk <- dxHMD[paste(x), paste(Y2)]

for (n in R2) {
  # Fit
  assign(x = paste0('M', n, '_', k), 
         value = fitMaxEntMortality(data = dxk, n = n, verbose = TRUE))
  # Predict
  assign(x = paste0('P', n, '_', k), 
         value = predict(get(paste0('M', n, '_', k)), h = S2, verbose = TRUE))
}

# Figure 2
plot(M6_GBRTENW, "observed", ny = 5)

# Figure 4
plot(P6_GBRTENW, M6_GBRTENW, plotType = "normalised_moments") +
  scale_x_continuous(breaks = c(1980, 2016, 2040)) +
  theme1()

# Figure 5
plot(P6_GBRTENW, "mean", ny = 5)

# ----------------------------------------------
# Back-Testing Age: 0 - 110

x = 0:95
D = dxHMD[paste(x), ]
MM <- c("MRWD", "LeeCarter", "HyndmanUllah", "CoDa", "MEM6")
BB <- doBBackTesting(data = D, x, y = 1960:2016,
                     data.in = "dx", 
                     models = MM,
                     strategy = c(f = 20, h = 20, s = 1), 
                     level = 95,
                     jumpchoice = "actual")

A <- evalAccuracy(BB, "ex")
R <- doRanking(A)

A[A$Scenario == "Total", ]
R[R$Scenario == "Total", ]


# Figure 6
ages = c(0, 20, 45, 65, 80, 95)
plot(BB$results$S18, data.out = "ex", facet = "x", which = ages) 
# Figure 7
plot(BB$results$S18, data.out = "mx", facet = "x", which = ages)


