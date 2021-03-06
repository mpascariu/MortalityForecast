# --------------------------------------------------- #
# Author: Marius D. PASCARIU
# License: GNU General Public License v3.0
# Last update: Tue Jan 22 10:07:39 2019 
# --------------------------------------------------- #

# MortalityForecast Package

#' @details \insertNoCite{*}{MortalityForecast}
#' @references \insertAllCited{}
#' @importFrom MortalityLaws convertFx MortalityLaw
#' @importFrom StMoMo StMoMo lc rh genWeightMat
#' @importFrom compositions acomp geometricmeanCol clr clrInv
#' @importFrom magrittr %>%
#' @importFrom moments all.moments raw2central
#' @importFrom demography demogdata fdm lca
#' @importFrom forecast forecast Arima arimaorder auto.arima
#' @importFrom stats predict fitted median coef lm na.omit quantile cov qnorm
#' lowess
#' @importFrom tidyr gather
#' @importFrom tibble as.tibble tibble add_column
#' @importFrom utils head tail sessionInfo
#' @importFrom forcats fct_recode fct_inorder
#' @importFrom ggplot2 ggplot geom_line  geom_point geom_ribbon facet_wrap    
#' theme theme_minimal element_text labs ylab xlab unit
#' aes alpha guide_legend guides margin
#' scale_color_manual scale_fill_manual scale_x_continuous scale_y_continuous
#' @importFrom ggridges stat_density_ridges
#' @importFrom graphics par plot matplot abline image image.default
#' @importFrom grDevices colorRampPalette grey.colors terrain.colors
#' @importFrom RColorBrewer brewer.pal
#' @importFrom fields image.plot
#' @importFrom pbapply startpb closepb setpb
#' @importFrom parallel makeCluster stopCluster clusterExport 
#' detectCores parLapply mclapply 
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar% foreach 
#' @name MortalityForecast
#' @aliases NULL
#' @docType package
"_PACKAGE"



