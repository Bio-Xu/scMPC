#####
# TODO Construct Covariance Matrix Based on Clustering Results
#' @description
#' Perform clustering based on the feature gene expression profile, and generate a covariance matrix between each pair of samples.
#' 1) Use Nbclust to select the optimal number of clusters.
#' 2) Perform K-means clustering.
#' 3) Generate the covariance matrix based on the clustering results.
#' 
#' @param results.CAM  A list that includes the covariance matrix, cluster information, etc.
#' @param nfeatures  The number of highly variable genes selected for dimensionality reduction, default is NULL (no dimensionality reduction).
#' @param erMatrix  The malignant cell expression profile.
#' @param distance The distance measure used when weighting samples, default is "euclidean"; other options include "kendall", "manhattan", "spearman", "person".
#' 
#' @return list: 1) The optimal number of clusters `optimalK`; 2) Covariance matrix `CAM` (a 01 matrix): 1 indicates samples are in the same cluster, 0 indicates they are not in the same cluster.
####


ConstructCAM <- function(i_clus){
    df.clust <- as.data.frame(i_clus) %>%
      tibble::rownames_to_column("cell") #将行名作为第一列并命名为cell
    outer_join <- merge(df.clust, df.clust, by='i_clus', all=TRUE)
    CAM <- table(outer_join$cell.x, outer_join$cell.y) %>%
      as.data.frame.array 
    t.filter <- setdiff(colnames(erMatrix), rownames(CAM)) #删除的样本
    CAM[t.filter, t.filter] <- 0
    CAM[is.na(CAM)] <- 0
    CAM <- as.matrix(CAM)
    diag(CAM) <- 1 
    CAM <- CAM[colnames(erMatrix), colnames(erMatrix)] 
    
    return(list(CAM = CAM))
  
}
