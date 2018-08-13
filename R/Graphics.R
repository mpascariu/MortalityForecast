

#' Plot method for objects of the class doBackTesting
#' @param x An object of the class \code{doBackTesting}.
#' @inheritParams wide2long
#' @export
plot.doBackTesting <- function(x, filter.x = c(0, 60, 70, 80, 90, 100), ...) {
  y = value = Name = DATA <- NULL # hack CRAN note
  
  B  <- x
  x  <- B$input$x
  y1 <- B$input$y.fit
  y2 <- B$input$y.for
  index <- B$input$what
  
  y_lab <- switch(index,
                  mx = "Central death rate at age x, \nm[x]",
                  qx = "Probability of dying between age x and x+1, \nq[x]",
                  dx = "Life Table d[x]",
                  lx = "Survivorship, \nl[x]",
                  Lx = "Lx",
                  Tx = "Tx",
                  ex = "Remaining Life Expectancy at age x, \ne[x]")
  
  # Observed values
  O <- convertFx(x, 
                 data = B$input$data, 
                 In = B$input$data.type, 
                 Out = B$input$what, 
                 lx0 = 1)
  O <- wide2long(data = O, x, filter.x) %>% mutate(DATA = NA)
  O[O$y %in% y1, "DATA"] <- "Training Set"
  O[O$y %in% y2, "DATA"] <- "Validation Set"
  
  # Forecast values
  H <- B$datasets$forecasts
  H <- wide.list.2.long.df(data = H, x, filter.x)
  
  # ggplot method
  P <- ggplot(H) + 
    facet_wrap(~x, scales = "free") +
    geom_line(aes(x = y, y = value, color = Name, linetype = Name)) +
    geom_point(data = O, aes(x = y, y = value, fill = DATA), shape = 21) +
    scale_fill_manual(values = 1:2) +
    guides(colour = guide_legend("Mortality Forecast\nModels"), 
           linetype = guide_legend("Mortality Forecast\nModels"), 
           fill = guide_legend("Demographic Data")) +
    xlab("Time period (Years)") + ylab(y_lab) +
    theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))
  P
}

