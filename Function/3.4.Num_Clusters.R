###############################################################################
# TODO: Perform Consensus Clustering Based on Signature Score Matrix and Automatically Determine the Optimal Cluster Number
#
#' @param scoreMatrix: Signature score matrix (e.g., ssGSEA), with rows as signatures and columns as samples.
#' @param signatureLabel: A vector of functional labels for the signatures, corresponding to the rows of scoreMatrix.
#' @param maxK: Maximum number of clusters, which can be specified based on annotation categories or the initial consensus clustering results.
#' @param minClusterNum: Minimum number of clusters, which can be specified based on annotation categories or the initial consensus clustering results.
#' @param outDir: Output directory for clustering results and images.
#' 
#' @return The optimal consensus clustering result.
###############################################################################

Num_Clusters <- function(res_CCP,  maxK = 8, minClusterNum = 2,outDir = "GSE158399"){
  #####
  # 2.1.1 PAC
  
  PACs <- computePAC(res_CCP, maxK)
  PAC_df <- data.frame(K = 2:c(length(PACs)+1), PAC = PACs)
  png(file.path(outDir,"PAC.png"))
  print(ggplot(PAC_df, mapping = aes(x = K, y = PAC, group = 1)) + 
          geom_line(col="steelblue") + geom_point() + theme_bw())
  dev.off()
  PACs2 <- PACs[-(1:(minClusterNum-1))]
  # The cluster number corresponding to the minimum PAC is optimal
  candidateK[1] <- as.numeric(names(PACs2)[which.min(PACs2)])
    
  #####
  # 2.1.2 Determine Optimal Cluster Number Based on Item Consensus
  list.IC <- AnalyzeICL(res_CCP, outDir)
  list.IC2 <- lapply(list.IC, function(x){
    x[-(1:(minClusterNum-1))]
  })
  # The cluster number corresponding to the maximum value is optimal
  candidateK[2] <- as.numeric(names(list.IC2$meanIC)[which.max(list.IC2$meanIC)])
  candidateK[3] <- as.numeric(names(list.IC2$meanCLC)[which.max(list.IC2$meanCLC)])
    
  #####
  # 2.1.3 Silhouette
  clusterLabels <- lapply(res_CCP[-1], function(x){
    x[["consensusClass"]]
  })
  distanceMatrix <- lapply(res_CCP[-1], function(x){
    1-x[["consensusMatrix"]]
  })
  silhouetteResults <- list()
  for(i in 1:length(clusterLabels)){
    silhouetteResults[[i]] <- calculateSilhouette(clusterLabels[[i]], distanceMatrix[[i]])
  }
  meanSilhuette <- sapply(silhouetteResults, function(x){
    x$avg.width
  })
  names(meanSilhuette) <- 2:maxK
  meanSilhuette_df <- data.frame(K = 2:c(length(meanSilhuette)+1), Silhuette = meanSilhuette)
  png(file.path(outDir,"Silhuette.png"))
  print(ggplot(meanSilhuette_df, mapping = aes(x = K, y = Silhuette, group = 1)) + 
          geom_line(col="steelblue") + geom_point() + theme_bw())
  dev.off()
  meanSilhuette2 <- meanSilhuette[-(1:(minClusterNum-1))]
   # The cluster number corresponding to the maximum value is optimal
  candidateK[4] <- as.numeric(names(meanSilhuette2)[which.max(meanSilhuette2)])
    
  #####
  # 2.1.4 Cophenetic Correlation Coefficient
  cccResults <- CopheneticCorCoef(res_CCP[-1], linkage = "complete", output = outDir)
  names(cccResults) <- 2:maxK
  cccResults2 <- cccResults[-(1:(minClusterNum-1))]
  #最大值对应的聚类数目为最优
  candidateK[5] <- as.numeric(names(cccResults2)[which.max(cccResults2)])
    
  #####
  # 2.2 Voting to Choose Optimal Cluster Number
  #####
  candidateK
  candidateK_ordered <- sort(table(candidateK), decreasing = T)
  if(any(candidateK_ordered > 1)){
    optimalK <- as.numeric(names(candidateK_ordered)[which(candidateK_ordered == candidateK_ordered[1])])
  }else optimalK <- median(candidateK)
  return(optimalK)
}
