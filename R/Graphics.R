

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
    guides(colour = guide_legend("Mortality Forecast\nModel"), 
           linetype = guide_legend("Mortality Forecast\nModel"), 
           fill = guide_legend("Demographic Data"))
  P
}

