# --------------------------------------------------- #
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Tue Nov 27 16:09:51 2018
# --------------------------------------------------- #

#' Plot method for objects of the class \code{BackTesting}
#' @param x An object of the class \code{do.BackTesting}.
#' @param facet What facets to include? Options: \code{"x", "y"}.
#' @param which Which \code{x} or which \code{y} to be plotted. Numerical vector.
#' @inheritParams evalAccuracy.BackTesting
#' @inheritParams wide2long
#' @export
plot.BackTesting <- function(x, 
                             data.out, 
                             facet = c("x", "y"),
                             which = NULL, 
                             ...) {
  B  <- x
  x  <- B$input$x
  y1 <- B$input$y.fit
  y2 <- B$input$y.for
  facet <- match.arg(facet)
  
  y_lab <- switch(data.out,
                  mx = "Central death rate at age x, \nm[x]",
                  qx = "Probability of dying between age x and x+1, \nq[x]",
                  dx = "Life Table d[x]",
                  lx = "Survivorship, \nl[x]",
                  Lx = "Lx",
                  Tx = "Tx",
                  ex = "The expectation of life at age x, \ne[x]")
  
  # Observed values
  O  <- get.Observed(x = x, 
                     data.in = B$input$data.in,    
                     data.out = data.out, 
                     data = B$input$data)
  
  # Forecast values
  H  <- get.Forecasts(object = B$Forecast, 
                      data.out = data.out) 
  
  # ggplot method
  if (facet == "y") {
    P <- plot_y_facets(O, H, x, y1, y2, which)
    
  } else {
    P <- plot_x_facets(O, H, x, y1, y2, which)
  }
  
  out = P + ylab(y_lab) +
    guides(colour = guide_legend("Mortality Forecast\nModels"), 
           linetype = guide_legend("Mortality Forecast\nModels"), 
           fill = guide_legend("Demographic Data")) +
    theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))
  out
}


#' @keywords internal
plot_x_facets <- function(O, H, x, y1, y2, which) {
  y = value = Name = DATA <- NULL # hack CRAN note
  
  if (is.null(which)) {
    which <- unique(floor(quantile(x, probs = seq(0, 1, length.out = 6))))
  }
  
  H <- wide.list.2.long.df(data = H, x = x, y = y2, which.x = which)
  H$Name <- change_model_factor_levels(H$Name)
  O <- wide2long(data = O, x = x, y = c(y1, y2), which.x = which)
  O$DATA <- NA
  O[O$y %in% y1, "DATA"] <- "Training Set"
  O[O$y %in% y2, "DATA"] <- "Validation Set"
  
  P <- ggplot(H) + 
    facet_wrap(~ x, scales = "free_y", nrow = 2) +
    geom_point(data = O, aes(x = y, y = value, fill = DATA), shape = 21) +
    geom_line(aes(x = y, y = value, color = Name, linetype = Name), size = .7) +
    scale_fill_manual(values = 1:2) +
    xlab("Time period (Years)") 
  P
}


#' @keywords internal
plot_y_facets <- function(O, H, x, y1, y2, which) {
  y = value = Name = DATA <- NULL # hack CRAN note
  
  if (is.null(which)) {
    which <- unique(floor(quantile(y2, probs = seq(0,1, length.out = 6))))
  }
  
  H <- wide.list.2.long.df(data = H, x = x, y = y2, which.y = which)
  H$Name <- change_model_factor_levels(H$Name)
  O <- wide2long(data = O, x = x, y = c(y1, y2), which.y = which)
  O$DATA <- "Validation Set"
  
  P <- ggplot(H) + 
    facet_wrap(~ y, scales = "free_y", nrow = 2) +
    geom_point(data = O, aes(x = x, y = value, fill = DATA), 
               shape = 21, alpha = 0.3) +
    geom_line(aes(x = x, y = value, color = Name, linetype = Name), size = .7) +
    scale_fill_manual(values = 2) +
    xlab("Age (x)") + 
    scale_y_continuous(trans = 'log10')
  P
}



#' Change factor levels in order to get nice legends in ggplots
#' @param vect A vector of the class \code{factor}.
#' @examples 
#' x <- factor(c("MRWD", "LeeCarter", "HyndmanUllah", "Oeppen", "MEM6"))
#' new.x <- change_model_factor_levels(x)
#' @keywords internal
#' @export
change_model_factor_levels <- function(vect) {
  X <- suppressWarnings(forcats::fct_recode(vect, 
                                            "M.Random-Walk w Drift" = "MRWD", 
                                            "Lee-Carter" = "LeeCarter", 
                                            "Li-Lee" = "LiLee", 
                                            "Hyndman-Ullah" = "HyndmanUllah", 
                                            "Renshaw-Haberman" = "RenshawHaberman", 
                                            "Oeppen" = "Oeppen",
                                            "Oeppen-C" = "OeppenC",
                                            "MEM-2" = "MEM2",
                                            "MEM-3" = "MEM3",
                                            "MEM-4" = "MEM4",
                                            "MEM-5" = "MEM5",
                                            "MEM-6" = "MEM6"))
  out <- fct_inorder(X)
  return(out)
}

