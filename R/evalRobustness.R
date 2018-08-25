# Fri Aug 24 21:20:42 2018 --------- Marius D. Pascariu ---

#' Generic for evaluating the robustness of the mortality models
#' tested with the \code{doBBackTesting} function.
#'  
#' @param object An object of the class \code{doBBackTesting}.
#' @inheritParams doMortalityModels
#' @keywords internal
#' @author Marius D. Pascariu
#' @examples 
#' # For examples go to ?doBBackTesting
#' @keywords internal
#' @export
evalRobustness = function(object, ...)
  UseMethod("evalRobustness")


#' Evaluating the Robustness of the Forcasting Mortality Models
#' 
#' @inheritParams evalAccuracy.doBBackTesting
#' @inheritParams computeRobustness
#' @seealso 
#' \code{\link{doBBackTesting}}
#' \code{\link{evalAccuracy.doBBackTesting}}
#' @author Marius D. Pascariu
#' @examples 
#' # For examples go to ?doBBackTesting
#' @export  
evalRobustness.doBBackTesting <- function(object,
                           data.out = c("qx", "mx", "dx", "lx", "Lx", "Tx", "ex"),
                           measures = c("MD", "MAD", "sMRAD"),
                           ...) {
  
  data.out <- match.arg(data.out)
  O <- object$results
  
  sn <- nrow(object$scenarios)  # no. of scenarios
  Mn <- object$input$models # Model names
  MN <- length(Mn) # no. of models tested
  
  R  <- tibble()
  RR <- 0
  N  <- -c(1:3)
  for (s in 1:(sn - 1)) {
    F1 <- getForecasts(O[[s]]$Forecast, data.out)
    F2 <- getForecasts(O[[s + 1]]$Forecast, data.out)
    # The benchmark is always the first model specified in Mn (`models` argument)
    b1 <- F1[[1]]
    b2 <- F2[[1]]
    scenario <- paste(s, s + 1, sep = "-")
    
    Rs <- tibble()
    for (m in 1:MN) {
      f1 <- F1[[m]]
      f2 <- F2[[m]]
      r <- computeRobustness(f1, f2, b1, b2, measures = measures)
      r <- add_column(Scenario = scenario, Model = Mn[m], 
                      LifeTableIndex = data.out, r, .before = T)
      Rs <- rbind(Rs, r)
    }
    R <- rbind(R, Rs)
    RR <- RR + Rs[, N]/(sn - 1)  # compute mean values over all scenarios
  }
  Rs[, N] <- RR
  Rs$Scenario <- "Total"
  out <- rbind(Rs, R)
  return(out)
}


#' Compute robustness measures
#' 
#' @param F1 Forecast scenario 1.
#' @param F2 Forecast scenario 2.
#' @param B1 Forecast scenario 1 - benchmark.
#' @param B2 Forecast scenario 2 - benchmark.
#' @param measures What robustness measure to compute? Various alternatives are 
#' available, \itemize{
#'  \item{Mean delta measures: } \code{"MD", "MAD", "sMRAD"};
#'  \item{Median delta measures: } \code{"MdD", "MdAD", "sMdRAD"}.
#'  }
#' @param na.rm A logical value indicating whether NA values should be stripped 
#' before the computation proceeds. Default: \code{TRUE}.
#' @author Marius D. Pascariu
#' @keywords internal
computeRobustness <- function(F1, F2, B1, B2, 
                              measures,
                              na.rm = TRUE) {
  # Identify the common years in the 2 forcasts
  y1 <- colnames(F1)
  y2 <- colnames(F2)
  y  <- y1[y1 %in% y2]
  # Subset the data so that we have the same years in the 2 forecasts
  F1 <- as.matrix(F1[, y])
  F2 <- as.matrix(F2[, y])
  B1 <- as.matrix(B1[, y])
  B2 <- as.matrix(B2[, y])
  
  # The % change in point forecast from A to B
  delta   = (F2/F1 - 1) * 100
  delta.b = (B2/B1 - 1) * 100
  
  N  <- na.rm
  D  <- delta
  bD <- delta.b
  AD <- abs(D)
  bAD  <- abs(bD)
  sRAD <- 2 * AD/(AD + bAD)
  
  MD <- mean(D, na.rm = N)
  MAD <- mean(AD, na.rm = N)
  sMRAD <- mean(sRAD, na.rm = N)
  
  MdD <- median(D, na.rm = N)
  MdAD <- median(AD, na.rm = N)
  sMdRAD <- median(sRAD, na.rm = N)
  
  out <- tibble(MD, MAD, sMRAD,     # mean deltas
                MdD, MdAD, sMdRAD) # median deltas
  out <- out[, measures]
  return(out)
}