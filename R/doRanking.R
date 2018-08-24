

#' Rank models based on the accuracy results
#' @param A Table containing accuracy measures.
#' @seealso 
#' \code{\link{doBackTesting}} 
#' \code{\link{doBBackTesting}}
#' @author Marius D. Pascariu
#' @examples 
#' # For examples go to ?doBackTesting or ?doBBackTesting
#' @export
doRanking <- function(A) {
  
  S  <- unique(A$Scenario)
  ns <- length(S) # no. of scenarios
  
  if (ns == 1) {
    N  <- -c(1:3)
    if (any(colnames(A) == "ME")) A[, "ME"] <- abs(A[, "ME"])
    A[, N] <- as.integer(floor(apply(A[, N], 2, rank)))
    G <- rank(apply(A[, N], 1, median))
    G <- as.integer(floor(G))
    R <- add_column(A, ModelRanking = G)
    
  } else {
    R <- tibble()
    for (s in 1:ns) {
      As <- A[A$Scenario == S[s], ]
      Rs <- doRanking(As)
      R  <- rbind(R, Rs)
    }
  }
  R
}




