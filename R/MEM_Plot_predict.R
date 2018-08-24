
#' ggplot the predicted values of a Maximum-Entropy Mortality Model
#' 
#' @param x An object of the class \code{\link{predict.fitMaxEntMortality}}.
#' @param y An object of the class \code{\link{fitMaxEntMortality}}. 
#' Needed only for ploting moments. See \code{plotType}.
#' @param plotType The type of the plot. The alternatives are 
#' \code{"mean", "lower", "upper", "raw_moments", "normalised_moments"}. 
#' Default: \code{"mean"}.
#' @inheritParams plot.fitMaxEntMortality
#' @examples
#' x  <- 0:110
#' y  <- 1965:2014
#' dx <- MortalityForecast.data$dx[paste(x), paste(y)]
#' M  <- fitMaxEntMortality(dx, x, y, n = 5)  
#' P  <- predict(M, h = 16, x.h = 0:120)
#' 
#' plot(P, plotType = "mean")
#' plot(P, plotType = "lower")
#' plot(P, plotType = "upper")
#' 
#' plot(P, M, plotType = "raw_moments")
#' plot(P, M, plotType = "normalised_moments")
#' @export
plot.predict.fitMaxEntMortality <- function(x, y = NULL,
                                plotType = c("mean", "lower", "upper", 
                                             "raw_moments", "normalised_moments"), 
                                ny = 7, level = 80, ...) 
{
  plotType <- match.arg(plotType)
  if (plotType == "mean") {
    mat = x$predicted.values
    P <- ggplotDistribConvergence(mat, x = x$x, ny, level) + 
      labs(subtitle = "Forecast Values - Best estimate")
    
  } else if (plotType == "lower") {
    mat = x$conf.intervals$predicted.values[[1]]
    P <- ggplotDistribConvergence(mat, x = x$x, ny, level) + 
      labs(subtitle = "Forecast Values - lower bound")
    
  } else if (plotType == "upper") {
    mat = x$conf.intervals$predicted.values[[2]]
    P <- ggplotDistribConvergence(mat, x = x$x, ny, level) + 
      labs(subtitle = "Forecast Values - upper bound")
    
  } else {
    P <- ggplotPredict(P = x, M = y, plotType = plotType)
  }
  suppressMessages(print(P))
}

# ----------------------------------------------


#' @keywords internal
ggplotPredict <- function(P, M, plotType = c("raw_moments", 
                                             "normalised_moments")) {
  year = value = lw = up = type <- NULL # hack CRAN note
  plotType <- match.arg(plotType)
  data <- prepare_ggdata(P, M, plotType)
  
  G <- ggplot(data, aes(year, value)) +
    geom_ribbon(aes(ymin = lw, ymax = up, fill = type)) +
    geom_line(aes(colour = type)) +
    facet_wrap(~ moment, scales = "free_y") +
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
  
  if (plotType == "raw_moments") fn = function(Z) Z
  if (plotType == "normalised_moments") {
    fn = function(Z) convertMoments(Z, "raw", "normalized")
  }
  
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




