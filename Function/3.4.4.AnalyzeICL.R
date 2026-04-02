###########
# TODO: Analyze Item Consensus and Cluster Consensus Averages for Different Cluster Numbers in Consensus Clustering Results
#
#' @description 
#' Use the `calcICL` function from the ConsensusClusterPlus package to calculate Item Consensus (IC) and Cluster Consensus (CLC).
#' 1) Calculate the average IC for each cluster pair and then calculate the average IC for each number of clusters.
#' 2) Calculate the average within-cluster CLC for each number of clusters.
#' Generate line plots to help determine the optimal number of clusters.
#' 
#' @param res_CCP: The result of consensus clustering (output from the `ConsensusClusterPlus` function).
#' @param outDir: The directory to output the plots.
#' 
#' @return A list containing the average IC: `meanIC` and average CLC: `meanCLC`.
###########
AnalyzeICL <- function(res_CCP, outDir = "./"){
  library(ConsensusClusterPlus)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  #####
  # 1. Calculate IC and CLC
  ####
  icl <- calcICL(res_CCP, title=outDir, plot="png")
  t.ic <- icl$itemConsensus
  t.clc <- icl$clusterConsensus
  
  #####
  # 2. Combine IC with consensus clustering results to calculate the average IC for each cluster
  ####
  t.class <- lapply(res_CCP, function(x){
    enframe(x[[3]], name = "item", value = "cluster")
  })
  t.class <- t.class[-1]
  names(t.class) <- 2:length(res_CCP)

  pairedIC <- list()
  for(i in names(t.class)){
    t.ic_k <- filter(t.ic, k==i)
    pairedIC[[i]] <- inner_join(t.ic_k, t.class[[i]]) #保留都出现的行，即样本与所归属的cluster所对应行
  }

  meanIC <- sapply(pairedIC, function(x){
    mean(x[,"itemConsensus"], na.rm=T)
  })

  meanIC_df <- enframe(meanIC, name = "K", value = "meanIC")
  meanIC_df$K <- as.numeric(meanIC_df$K)
  png(file.path(outDir,"meanIC.png"))
  print(ggplot(meanIC_df, mapping = aes(x = K, y = meanIC, group = 1)) + geom_line(col="steelblue") + geom_point() + theme_bw())
  dev.off()
  
  #####
  # 3. Calculate average CLC within each cluster
  ####
  meanCLC <- tapply(t.clc[,3], t.clc[,1], mean, na.rm = T)
  meanCLC_df <- enframe(meanCLC, name = "K", value = "meanCLC")
  meanCLC_df$K <- as.numeric(meanCLC_df$K)
  png(file.path(outDir,"meanCLC.png"))
  print(ggplot(meanCLC_df, mapping = aes(x = K, y = meanCLC, group = 1)) + geom_line(col="steelblue") + geom_point() + theme_bw())
  dev.off()

  results <- list(meanIC = meanIC, meanCLC = meanCLC)
  return(results)
}