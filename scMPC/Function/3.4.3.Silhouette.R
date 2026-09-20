############################################
# TODO: Calculate Clustering Evaluation Metric - Silhouette Coefficient
#
#' @description 
#' Calculate the silhouette width for each sample and compute the average silhouette coefficient for all samples and each cluster.
#' 
#' @param cluster_labels A vector of cluster assignments (the vector name is the sample, and the value is the cluster label).
#' @param distance_matrix A distance matrix representing pairwise distances between samples.
#' @output list, including the following three elements:
#' 1. `silhouette`: The result of the silhouette function (rows represent samples, the `sil_width` column represents the silhouette width, as shown below):
#         cluster neighbor  sil_width
# 01005       1        3    0.94646509
# 01010       2        3    0.86871323
# 03002       1        3    0.92158098
#' 2. `avg.width`: The average silhouette width for all samples.
#' 3. `clus.avg.widths`: A vector of the average silhouette width values for each cluster (as shown below):
#         1         2         3         4 
# 0.8741171 0.7679973 0.7694033 0.9964693 
###########################################

calculateSilhouette <- function(cluster_labels, distance_matrix, ...){
  library(cluster)
  result <- list()
  # 计算并转换silhouette值
  sil_width <- silhouette(cluster_labels, distance_matrix, ...)
  result$silhouette <- as.data.frame.matrix(sil_width)
  rownames(result$silhouette) <- names(cluster_labels)
  # 获取所有样本的平均silhouette值
  result$avg.width <- summary(sil_width)$avg.width
  # 获取每个cluster的平均silhouette值
  result$clus.avg.widths <- summary(sil_width)$clus.avg.widths
  return(result)
}

