

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
                  ex = "The expectation of life at age x, \ne[x]")
  
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


# ----------------------------------------------


#' Plot method for objects of the class \code{getAccuracy}
#' @param x An object of the class \code{getAccuracy}.
#' @param fill Colour for filling the box of the best-performing method.
#' @param lwd Line width for marking the box of the best-performing method.
#' @param ... Ignored.
#' @seealso \code{\link{getAccuracy}}
#' @export
plot.getAccuracy <- function(x, fill = "#addd8e", lwd = 0.5, ...) {
  A  <- x
  my_theme <- ttheme_minimal(base_family = "serif")
  G1 <- tableGrob(round(A$results, digits = 4), theme = my_theme)
  G2 <- tableGrob(A$rank, theme = my_theme)
  G3 <- tableGrob(as.matrix(sort(A$GC)), theme = my_theme)
  
  t1 <- tableGrob(paste("Forecasting Accuracy Measures\nLife Table Index:", A$index),
                  theme = my_theme)
  t2 <- tableGrob('Ranking: Best performing models in each category',
                  theme = my_theme)
  t3 <- tableGrob("General\nClassification",
                  theme = my_theme)
  
  wth <- rep(unit(.12, "npc"), ncol(A$results) + 1)
  G1$widths <- G2$widths <- wth
  
  w <- find_grob_cells(G1, A)
  for (i in 1:length(w)) {
    G1$grobs[w][[i]][["gp"]] <- gpar(fill = fill, lwd = lwd)
    G2$grobs[w][[i]][["gp"]] <- gpar(fill = fill, lwd = lwd)
  }
  
  # very inefficient code here --- clearly I don't know what I am doing with 
  # these grobs
  d = dim(G1)
  G1 <- gtable_add_grob(G1, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                        t = 2, b = d[1], l = 1, r = d[2])
  G1 <- gtable_add_grob(G1, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                        t = 1, b = 1, l = 1, r = d[2])
  G2 <- gtable_add_grob(G2, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                        t = 2, b = d[1], l = 1, r = d[2])
  G2 <- gtable_add_grob(G2, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                        t = 1, b = 1, l = 1, r = d[2])
  
  grid.newpage()
  lay <- rbind(c(1,1,1,3,3,3,5),
               c(2,2,2,4,4,4,6),
               c(2,2,2,4,4,4,6),
               c(2,2,2,4,4,4,6),
               c(2,2,2,4,4,4,6))
  grid.arrange(t1, G1, t2, G2, t3, G3, layout_matrix = lay)
}


#' Utility funtion for plot.getAccuracy
#' @param tg An object of the class \code{\link{tableGrob}}.
#' @param object An object of the class \code{\link{getAccuracy}}.
#' @keywords internal
find_grob_cells <- function(tg, object){
  V  <- object$results
  R  <- object$rank
  N  <- R * 0 + 1:nrow(R) + 1
  t  <- N[R == 1]
  l  <- 1:ncol(V) + 1
  l  <- rep(l, colSums(R == 1))
  tl <- paste(t, l, sep = "-")
  
  L1 <- tg$layout
  L1$tl <- paste(L1$t, L1$l, sep = "-")
  
  w <- which(L1$tl %in% tl & L1$name == "core-bg")
  return(w)
}

