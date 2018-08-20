
# MortalityForecast Package

#' @import dxForecast CoDa
#' @importFrom StMoMo StMoMo lc genWeightMat
#' @importFrom forecast forecast
#' @importFrom tidyr gather
#' @importFrom MortalityLaws convertFx
#' @importFrom demography demogdata fdm
#' @importFrom ggplot2 aes element_text facet_wrap geom_line geom_point ggplot 
#' guide_legend guides margin scale_fill_manual scale_y_continuous theme unit 
#' xlab ylab
#' @importFrom grid textGrob gpar gTree gList grid.newpage rectGrob
#' @importFrom gridExtra tableGrob grid.arrange ttheme_minimal
#' @importFrom gtable gtable_add_grob
#' @importFrom stats predict fitted median coef lm na.omit quantile
#' @importFrom utils tail
#' @name MortalityAccuracy-package
#' @docType package
"_PACKAGE"
