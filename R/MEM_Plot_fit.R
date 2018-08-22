
#' ggplot the observed and fitted values of a Maximum-Entropy Mortality Model 
#' 
#' @param plotType The type of the plot. The alternatives are 
#' \code{"fitted"}, \code{"observed"}. Default: \code{"fitted"}.
#' @inheritParams plot.residuals.MEM
#' @inheritParams fitted2dens
#' @param quant Sets the number of quantiles the data should be broken into. 
#' Default: c(0.1, 0.9).
#' @details Note: the output is a ggplot2 object. Therefore, one can add or 
#' modify the components of the figure.
#' @examples 
#' x  <- 0:110
#' y  <- 1965:2014
#' dx <- MortalityForecast.data$dx[paste(x), paste(y)]
#' M  <- MEM(dx, x, y, n = 5)
#' 
#' plot(M, plotType = "fitted")
#' plot(M, plotType = "observed")
#' @export   
plot.MEM <- function(x, plotType = c("fitted", "observed"), 
                        ny = 7, quant = c(0.1, 0.9), ...) 
{
  plotType <- match.arg(plotType)
  if (plotType == "fitted") {
    mat = x$fitted.values
    P <- ggplotDistribConvergence(mat, x = x$x, ny, quant) + 
      labs(subtitle = "Fitted Values")
    
  } else if (plotType == "observed") {
    mat = x$observed.values
    P <- ggplotDistribConvergence(mat, x = x$x, ny, quant) + 
      labs(subtitle = "Observed Values")
  } 
  suppressMessages(print(P))
}


#' Plot Convergence in Age at Death Distribution
#' 
#' @inheritParams fitted2dens
#' @inheritParams MEM
#' @keywords internal
ggplotDistribConvergence <- function(mat, x, ny, quant) {
  dx = y = ..quantile.. <- NULL # hack CRAN note
  Z <- mat %>% fitted2dens(ny = ny)
  rx <- range(x)
  
  P <- ggplot(Z, aes(x = dx, y = y, fill = factor(..quantile..))) +
    stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = T, 
                        quantiles = quant, quantile_lines = F) +
    scale_fill_manual(name = "Prob.",
                      values = c("#FF0000A0", "#A0A0A0A0", "#0000FFA0"),
                      labels = c(paste0("=< ", quant[1]*100, "%"), "IQR", 
                                 paste0("> ", quant[2]*100, "%"))) +
    scale_x_continuous(limits = rx,
                       breaks = seq(rx[1], rev(rx)[1], by = 10),
                       expand = c(0.05, 0)) +
    labs(x = "Age (x)", y = "Period (y)",
         title = "Convergence in Age at Death Distribution") +
    theme_minimal(base_family = "mono") +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 15))
}


#' Prepare data for ggplots in plot.MEM function
#' 
#' @param mat Matrix containing the observed or fitted value.
#' @param ny Number of years to be selected from input data and to be added in the plot.
#' @param lx0 Adjusting parameter.
#' @examples 
#' x  <- 0:110
#' y  <- 1965:2014
#' dx <- MortalityForecast.data$dx[paste(x), paste(y)]
#' M  <- MEM(dx, x, y, n = 5)
#' fitted2dens(fitted(M))
#' @keywords internal
#' @export
fitted2dens <- function(mat, ny = 7, lx0 = 300) {
  
  x <- as.numeric(rownames(mat))
  y <- as.numeric(colnames(mat))
  years <- rev(round(quantile(y[-1], probs = seq(0, 1, length.out = ny))))
  
  H <- cbind(x, mat) %>% as.data.frame() %>% 
    tidyr::gather(key = "Year", value = "dx", -x)
  H <- H[H$Year %in% years, ]
  H$Dx <- round(H$dx * lx0)
  
  out <- NULL
  for (i in years) {
    X0 <- H[H$Year == i, ]
    X1 <- NULL
    for (j in 1:nrow(X0)) {
      vect <- rep(X0[j, "x"], X0[j, "Dx"])
      X1 <- c(X1, vect)
    }
    out <- rbind(out, data.frame(y = i, dx = X1))
  }
  out$y <- as.factor(out$y)
  return(out)
}
