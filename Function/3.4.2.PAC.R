#################################################
# TODO: Calculate Proportion of Ambiguous Clustering (PAC)
#' @description 
#' The PAC method can be used to determine the optimal number of clusters (K) for consensus clustering. The smaller the PAC, the better the value of K.
#' "Proportion of Ambiguous Clustering" (PAC) is defined as the proportion of sample pairs with the same label that fall within the middle subinterval (x1, x2) ∈ [0, 1] of the CDF curve.
#' A lower PAC value indicates a flatter middle segment, and a flat CDF curve occurs only when the true K is correct.
#' 
#' @param consensus_result: The result of consensus clustering (output from the ConsensusClusterPlus function).
#' @param maxK: The maximum K value used in consensus clustering.
#' @param lower and upper: The lower and upper bounds that determine what is ambiguous.
#'                         In a perfect clustering, the consensus matrix would consist of only 0s and 1s, 
#'                         and the PAC assessed on the (0, 1) interval would have a perfect score of 0. 
#'                         A (0.1, 0.9) interval is commonly used to define ambiguity.
#' 
#' @output A numeric vector storing the PAC values for different K values.
#################################################

computePAC <- function(consensus_result, maxK, lower = 0.1, upper = 0.9) {
	Kvec <- 2:maxK
	PAC <- rep(NA, length(Kvec))
	names(PAC) <- Kvec
	for (i in Kvec) {
		M <- consensus_result[[i]]$consensusMatrix
		Fn <- ecdf(M[lower.tri(M)])
		PAC[i - 1] <- Fn(upper) - Fn(lower)
	}
	return(PAC)
}
