# Fri Aug 24 12:35:23 2018 --------- Marius D. Pascariu ---


#' ggplot the observed and fitted values of a Maximum-Entropy Mortality Model 
#' 
#' @param plotType The type of the plot. The alternatives are 
#' \code{"fitted"}, \code{"observed"}. Default: \code{"fitted"}.
#' @inheritParams plot.residMF
#' @inheritParams fitted2dens
#' @param level Sets the number of quantiles the data should be broken into. 
#' Default: 80.
#' @details Note: the output is a ggplot2 object. Therefore, one can add or 
#' modify the components of the figure.
#' @seealso \code{\link{fit_MEM}}
#' @examples 
#' # For examples go to ?fit_MEM
#' @export   
plot.MEM <- function(x, plotType = c("fitted", "observed"), 
                                    ny = 7, level = 80, ...) 
{
  plotType <- match.arg(plotType)
  if (plotType == "fitted") {
    mat = x$fitted.values
    P <- ggplotDistribConvergence(mat, x = x$x, y = x$y[-1], ny, level) + 
      labs(subtitle = "Fitted Values")
    
  } else if (plotType == "observed") {
    mat = x$observed.values
    P <- ggplotDistribConvergence(mat, x = x$x, y = x$y, ny, level) + 
      labs(subtitle = "Observed Values")
  } 
  suppressMessages(print(P))
}


#' Plot Convergence in Age at Death Distribution
#' 
#' @inheritParams fitted2dens
#' @inheritParams fit_MEM
#' @keywords internal
ggplotDistribConvergence <- function(mat, x, y, ny, level) {
  dx = ..quantile.. <- NULL # hack CRAN note
  
  Z  <- fitted2dens(mat, x, y, ny = ny)
  rx <- range(x)
  quant <- abs(c(0, -1) + (1 - level/100)/2)
  
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
#' @inheritParams doMortalityModels
#' @examples 
#' x  <- 0:110
#' y  <- 1965:2014
#' dx <- HMD_male$dx$GBRTENW[paste(x), paste(y)]
#' M  <- fit_MEM(dx, x, y, n = 5)
#' fitted2dens(fitted(M), x, y[-1])
#' @keywords internal
#' @export
fitted2dens <- function(mat, x, y, ny = 7, lx0 = 300) {
  
  years <- rev(round(quantile(y, probs = seq(0, 1, length.out = ny))))
  
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




#' ggplot the predicted values of a Maximum-Entropy Mortality Model
#' 
#' @param x An object of the class \code{\link{predict.MEM}}.
#' @param y An object of the class \code{\link{fit_MEM}}. 
#' Needed only for ploting moments. See \code{plotType}.
#' @param plotType The type of the plot. The alternatives are 
#' \code{"mean", "lower", "upper", "raw_moments", "normalised_moments", "scaled_moments"}. 
#' Default: \code{"mean"}.
#' @inheritParams plot.MEM
#' @examples
#' # For examples go to ?predict.MEM
#' @export
plot.predict.MEM <- function(x, y = NULL,
                             plotType = c("mean", "lower", "upper", 
                                          "raw_moments", "normalised_moments", "scaled_moments"), 
                             ny = 7, level = 80, ...) 
{
  plotType <- match.arg(plotType)
  if (plotType == "mean") {
    mat = x$predicted.values
    P <- ggplotDistribConvergence(mat, x = x$x, y = x$y, ny, level) + 
      labs(subtitle = "Forecast Values - Best estimate")
    
  } else if (plotType == "lower") {
    mat = x$conf.intervals$predicted.values[[1]]
    P <- ggplotDistribConvergence(mat, x = x$x, y = x$y, ny, level) + 
      labs(subtitle = "Forecast Values - lower bound")
    
  } else if (plotType == "upper") {
    mat = x$conf.intervals$predicted.values[[2]]
    P <- ggplotDistribConvergence(mat, x = x$x, y = x$y, ny, level) + 
      labs(subtitle = "Forecast Values - upper bound")
    
  } else {
    P <- ggplotPredict(P = x, M = y, plotType = plotType)
  }
  suppressMessages(print(P))
}

# ----------------------------------------------


#' @keywords internal
ggplotPredict <- function(P, M, plotType = c("raw_moments", 
                                             "normalised_moments",
                                             "scaled_moments")) {
  year = value = lw = up = type <- NULL # hack CRAN note
  plotType <- match.arg(plotType)
  data <- prepare_ggdata(P, M, plotType)
  
  G <- ggplot(data, aes(year, value)) +
    geom_ribbon(aes(ymin = lw, ymax = up, fill = type)) +
    geom_line(aes(colour = type)) +
    facet_wrap(~ moment, scales = "free_y", nrow = 2) +
    scale_color_manual(name = "Legend: ", values = c("blue", 1)) +
    scale_fill_manual(name = "Legend: ", values = alpha(c(4, NA), c(0.2, NA))) +
    xlab("\nTime (years)") + ylab("") 
}


#' @keywords internal
prepare_ggdata <- function(P, M, plotType) {
  moment = value <- NULL # hack CRAN note
  
  L <- list(oRM = M$observed.raw.moments, 
            pRM = P$predicted.raw.moments, 
            uRM = P$conf.intervals$predicted.raw.moments[[1]],
            lRM = P$conf.intervals$predicted.raw.moments[[2]])
  
  fn <- switch(plotType,
               raw_moments = function(Z) Z,
               normalised_moments = function(Z) convertMoments(Z, "raw", "normalized"),
               scaled_moments = function(Z) log(abs(convertMoments(Z, "raw", "normalized"))))
  
  L1 <- lapply(L, as.data.frame)
  L2 <- lapply(L1, fn)
  Mnames <- colnames(L2[[1]])
  
  fn2 = function(Z) {
    D <- data.frame(year = as.numeric(rownames(Z)), Z)
    colnames(D) <- c("year", Mnames)
    D
  }
  L3 <- lapply(L2, fn2)
  L3[[1]]$type <- "Observed"
  L3[[2]]$type <- "Forecast"
  
  H <- rbind(L3[[1]], L3[[2]])
  H <- tidyr::gather(H, moment, value, -c("year", "type"))
  H$moment <- factor(H$moment, levels = Mnames)
  H$lw = H$up <- NA
  H[H$type == 'Forecast', "lw"] <- unlist(L3[[3]][, -1])
  H[H$type == 'Forecast', 'up'] <- unlist(L3[[4]][, -1])
  H$year = as.numeric(H$year)
  return(H)
}





