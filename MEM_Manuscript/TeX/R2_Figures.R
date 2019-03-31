# --------------------------------------------------- #
# Author: Marius D. PASCARIU
# License: GNU General Public License v3.0
# Last update: Tue Mar 19 15:54:31 2019
# --------------------------------------------------- #
remove(list = ls())
library(tidyverse)
library(MortalityForecast)
load("dxForecastResults.Rdata")


plot_coverage <- function(cntr, yr){

  models <- paste0('C', n1,'_', cntr)
  dens   <- paste0('Moments (N = ', n1, ')')

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

  # dx.obs$rank <- 'Density'
  # DX <- rbind(DX, dx.obs)
  # ggplot
  L <- "Age-at-Death Distribution"
  P = ggplot(DX, aes(x = age, y = fx)) +
    geom_ribbon(aes(ymin = 0, ymax = fx, fill = type, alpha = type, color = type)) +
    facet_wrap(~ rank, nrow = 1) +
    scale_alpha_manual(name = L, values = c(0.12, 0.5)) +
    scale_fill_manual(name = L, values = c(4, 1)) +
    scale_color_manual(name = L, values = c(4, NA)) +
    scale_x_continuous(breaks = c(0, 40, 80, 120)) +
    xlab("Age, x") + ylab(expression(f(x))) +
    geom_text(data = coverage, inherit.aes = FALSE, size = 4, hjust = "inward",
              aes(x = 5, y = .027, label = paste0('Coverage:\n  ', covr,'%') ))
  P
}


theme1 <- function(){
  theme_bw() +
    theme(panel.background = element_rect(fill = 'grey98'),
          strip.text = element_text(colour = 1, size = 13),
          legend.background = element_rect(fill = NA),
          legend.text.align = 0,
          legend.key.size   = unit(0.8, "cm"),
          legend.position = "top",
          legend.text  = element_text(colour = 1, size = 13),
          legend.title = element_text(colour = 1, size = 13),
          axis.text.x  = element_text(colour = 1, size = 13),
          axis.text.y  = element_text(colour = 1, size = 10),
          axis.title   = element_text(colour = 1, size = 14, face = 'bold'),
          axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
          axis.line.x  = element_line(linetype = 1),
          axis.ticks   = element_line(size = 1),
          panel.spacing.x = unit(.5, "cm"),
          panel.spacing.y = unit(.2, "cm"),
          plot.margin   = unit(c(.1, .1, .1, .1), "cm")
    )
}

theme2 <- function(){
  theme1() +
    theme(legend.position = "right",
      legend.text  = element_text(colour = 1, size = 10),
      legend.title = element_text(colour = 1, size = 12))
}

# ----------------------------------------------
fh = 6    # height 1
fh2 = 4.5 # height 2
fw = 11   # width

# Figure 2 ----------------------------------------------
pdf("figure/Figure-Convergence.pdf", height = 3.5, width = fw)
print(plot_coverage("USA", 1990) + theme1())
dev.off()

# Figure 1 ----------------------------------------------
fg = plot(M6_GBRTENW_0, "observed", ny = 5)
fg = fg + labs(title = "", subtitle = "", x = "Age", y = "Year")
pdf("figure/Figure-ObservedDx.pdf", height = fh, width = fw)
fg
dev.off()

# Figure 3 ----------------------------------------------
fg = plot(P6_GBRTENW, M6_GBRTENW, plotType = "normalised_moments")
fg = fg + scale_x_continuous(breaks = c(1980, 2016, 2040)) +
  theme1()

pdf("figure/Figure-Moments.pdf", height = fh, width = fw)
fg
dev.off()


# Forecast Figures --- Back-Test ------------------------------


format_forecast_plot <- function(G, breaks) {

  G +
  scale_x_continuous(breaks = breaks) +
    scale_linetype_manual(values = c(1,2,5,1,2,1)) +
    scale_color_manual(values = c(2,1,3,8,6,4)) +
    guides(fill = guide_legend(title = "Demographic Data", keyheight = 1.3),
           color = guide_legend(title = "Mortality Forecasting\nModels"),
           linetype = guide_legend(title = "Mortality Forecasting\nModels",
                                   keyheight = 1.3)) +
    theme2()
}


ages <- c(0, 25, 45, 65, 75, 85)
for (K in cntr) {
  BT <- B[[K]]$results
  Sn <- BT[[length(BT)]]
  b <- c(range(Sn$input$y.fit), max(Sn$input$y.for))

  F_ex <- plot(Sn, data.out = "ex", facet = "x", which = ages)
  F_mx <- plot(Sn, data.out = "mx", facet = "x", which = ages)
  fig_ex <- format_forecast_plot(F_ex, breaks = b)
  fig_mx <- format_forecast_plot(F_mx, breaks = b)

  file1 <- paste0("figure/Figure_", K, "_ex.pdf")
  file2 <- paste0("figure/Figure_", K, "_mx.pdf")

  pdf(file1, height = fh2, width = fw); print(fig_ex); dev.off()
  pdf(file2, height = fh2, width = fw); print(fig_mx); dev.off()
}


