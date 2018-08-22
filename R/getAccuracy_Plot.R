

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
