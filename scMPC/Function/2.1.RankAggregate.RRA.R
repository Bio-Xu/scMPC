RankAggregate.RRA <- function(geneList, geneSet, full = FALSE, exact = FALSE){
  library(RobustRankAggreg)
  
  N <- length(unlist(geneSet))
  
  geneRank.RRA <- aggregateRanks(geneList, N = N, full = full, exact = exact)
  
  t.gene <- as.data.frame(table(unlist(geneList)))
  geneRank.RRA$Freq <- t.gene[match(geneRank.RRA$Name, t.gene[,1]), 2]
  
  return(geneRank.RRA)
}