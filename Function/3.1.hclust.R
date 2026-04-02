#####
# TODO: Hierarchical Clustering
#' @description
#' Perform hierarchical clustering based on the feature gene expression profile.
#' 1) Use Nbclust to select the optimal number of clusters.
#' 2) Perform hierarchical clustering using hclust.
#' 
#' @param erMatrix  The feature gene expression matrix, with features as rows and samples as columns.
#' @param nc.distance The distance measure used by NbClust, default is "euclidean".
#' @param nc.min The minimum number of clusters for NbClust.
#' @param nc.max The maximum number of clusters for NbClust.
#' @param nc.method The clustering method used by NbClust, default is "ward.D2".
#' @param seed Random seed, default value is 1234
#' 
#' @param index_list The indices used by NbClust to determine the best number of clusters, excluding those that cannot be computed, graphical methods, and slow computations. 17 indices are retained to jointly determine the optimal number of clusters.
#' 
#' @return list: 1) The optimal number of clusters, `optimalK`; 2) The cluster information for each sample, `i_clus`.
#####


hclustfunc <- function(erMetrix, nc.distance = "euclidean", nc.min = 2, nc.max = 8, nc.method = "ward.D2",seed = 1234){
  library(dplyr)
  library(NbClust)
  library(purrr)
  library(Seurat)
  library(mclust)
  #####
  # 0 Expression profile preprocessing
  ####
  erMetrix2 <- erMetrix[, colSums(abs(erMetrix))!= 0]
  erMetrix2 <- erMetrix2[rowSums(abs(erMetrix2))!= 0, ]
  erMetrix2 <- t(erMetrix2) %>% as.data.frame
  
  #####
  # 1 Perform hierarchical clustering
  ####
  set.seed(seed)
  index_list <- c("kl", "ch", "hartigan", "cindex", "db", "silhouette", "ratkowsky", "ball", "frey", "mcclain", "dunn", "sdindex")
  Clus <- c()
  
  for (i in 1:length(index_list)) {
    nb_clust <- tryCatch({
      NbClust(erMetrix2, distance = nc.distance, min.nc = nc.min, max.nc = nc.max, method = nc.method, index = index_list[i])
    }, error = function(e) {
      message(paste("Error in NbClust for index", index_list[i], ":", e))
      return(NULL)
    })
    
    if (!is.null(nb_clust)) {
      clus <- as.numeric(nb_clust$Best.nc[1])
      Clus <- c(Clus, clus)
    }
  }
  if (length(Clus) == 0) {
    message("No valid clustering results found, using default k=2.")
    optimalK <- 2
  } else {
    t.clust <- table(Clus)
    optimalK <- names(t.clust)[which.max(t.clust)] %>% as.numeric
  }

  baseCluster <-hclust(dist(erMetrix2), method = nc.method)
  
  i_clus<-cutree(baseCluster, k = optimalK)

    
  return(list(optimalK = optimalK,hclust_clus=i_clus))
}
