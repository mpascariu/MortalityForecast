

#' Rank models based on the accuracy or robustness results
#' @param A Table containing accuracy or robustness measures.
#' @seealso 
#' \code{\link{do.BackTesting}} 
#' \code{\link{do.BBackTesting}}
#' @author Marius D. Pascariu
#' @examples 
#' # For examples go to ?do.BackTesting or ?do.BBackTesting
#' @export
do.Ranking <- function(A) {
  
  S  <- unique(A$Scenario)
  ns <- length(S) # no. of scenarios
  
  if (ns == 1) {
    N  <- -c(1:3)
    if (any(colnames(A) == "ME")) A[, "ME"] <- abs(A[, "ME"])
    if (any(colnames(A) == "MD")) A[, "MD"] <- abs(A[, "MD"])
    
    A[, N] <- as.integer(floor(apply(A[, N], 2, rank)))
    G <- rank(apply(A[, N], 1, median))
    G <- as.integer(floor(G))
    R <- add_column(A, ModelRanking = G)
    
  } else {
    R <- tibble()
    for (s in 1:ns) {
      As <- A[A$Scenario == S[s], ]
      Rs <- do.Ranking(As)
      R  <- rbind(R, Rs)
    }
  }
  R
}




