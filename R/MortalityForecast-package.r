
# MortalityForecast Package

#' @details \insertNoCite{*}{MortalityForecast}
#' @references \insertAllCited{}
#' @import dxForecast
#' @importFrom StMoMo StMoMo lc genWeightMat
#' @importFrom compositions acomp geometricmeanCol clr clrInv
#' @importFrom MortalityLaws convertFx
#' @importFrom demography demogdata fdm
#' @importFrom forecast forecast Arima arimaorder auto.arima
#' @importFrom stats predict fitted median coef lm na.omit quantile cov qnorm
#' @importFrom tidyr gather
#' @importFrom utils head tail
#' @importFrom ggplot2 aes element_text facet_wrap geom_line geom_point ggplot 
#' guide_legend guides margin scale_fill_manual scale_y_continuous theme unit 
#' xlab ylab
#' @importFrom grid textGrob gpar gTree gList grid.newpage rectGrob
#' @importFrom gridExtra tableGrob grid.arrange ttheme_minimal
#' @importFrom gtable gtable_add_grob
#' @importFrom graphics par plot matplot abline image.default
#' @importFrom grDevices colorRampPalette grey.colors terrain.colors
#' @importFrom RColorBrewer brewer.pal
#' @importFrom fields image.plot
#' @name MortalityForecast
#' @aliases NULL
#' @docType package
"_PACKAGE"



