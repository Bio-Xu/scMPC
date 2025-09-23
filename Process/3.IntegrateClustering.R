##################################################################
# TODO: Integrate clustering results from multiple signature gene sets to identify metastasis-related cell subpopulations
#
#' @description
#' Performs k-means clustering on each signature's single-cell expression profile to select optimal cluster number k as base clustering results, generating a co-association matrix.
#' Improves the co-association matrix through two-level refinement at cell and cluster levels. Finally, performs path-based transformation and final clustering.
#' 
#' @param seuratObj: Seurat object storing processed single-cell expression profiles
#' @param signatureList: List of signature sets for base clustering
#' @param nc.idents: Cell subpopulation identifiers for malignant cells
#' @param method: Final clustering method for consensus clustering 
#'                1) hclust: Hierarchical clustering
#'                2) kmeans: K-means clustering
#'                3) PAM: Partitioning Around Medoids clustering
#' @param k: Number of clusters for final clustering
#' @param hc.method: Hierarchical clustering method, default "ward.D2". Options: "complete", "single", "average", "mcquitty", "median", "centroid"
#' @param file: Storage path for base clustering results
#' 
#' @return Group: Matrix with cell IDs as row names and clusters as column names
##################################################################
GroupByBaseCluster <- function(seuratObj, signatureList, nc.idents,
                           method = "kmeans", k = 3, hc.method = "ward.D2", file = "./"){
  library(Seurat)
  library(dplyr)
  library(NbClust)
  library(plyr)
  library(purrr)
  library(cluster)
  library(FCPS)
  
  #####
  # 0 Extract malignant cell expression profiles
  ####
  malignantObj <- subset(seuratObj, idents = nc.idents) # Extract malignant cells
  erMatrix <- as.data.frame(malignantObj@assays$RNA@data) # Extract cell expression matrix
  
  #####
  # I Build co-association matrix for base clustering
  ####
  # Generate co-association matrix for each signature
  results.CAM <- lapply(signatureList, function(x){
    ErMatrix <- erMatrix[x,] %>% na.omit()
    t.result <- ConstructCAM(ErMatrix)
  })
  save(results.CAM, file = file)
  
  #####
  # II Optimize co-association matrix
  ####
  optimize <- TwoLevelRefine(results.CAM, malignantObj, nfeatures = 500, erMatrix, distance = "euclidean")

  #####
  # III Perform clustering
  ####
  if(method == "hclust"){
    # 1 Hierarchical clustering
    fit <- hclust(as.dist(optimize), method = hc.method) # Perform hierarchical clustering
    cut0 <- cutree(fit, k = k)  # Manual cluster specification
    Group <- as.matrix(cut0) # Sample cluster information
  } else if(method == "kmeans"){
    # 2 K-means clustering
    Group <- kmeansClustering(optimize, ClusterNo = k, RandomNo = 10)
    Group <- Group$Cls %>% as.matrix()
  } else if(method == "pam"){
    # 3 PAM clustering
    Group <- pam(optimize, k, diss = TRUE)
    Group <- Group$clustering %>% as.matrix()
  }
  colnames(Group) <- "group"
  return(Group)
}

#####
# Function 1: Build co-association matrix
#'
#' @description
#' Performs base clustering on signature gene expression profiles to generate co-association matrix between samples.
#' 1) Uses NbClust to determine optimal cluster number
#' 2) Performs base clustering using k-means
#' 3) Generates co-association matrix from clustering results
#' 
#' @param erMatrix: Signature gene expression matrix (features as rows, samples as columns)
#' @param nc.distance: Distance measure for NbClust, default "euclidean"
#' @param nc.min: Minimum cluster number for NbClust
#' @param nc.max: Maximum cluster number for NbClust
#' @param nc.method: Clustering method for NbClust, default "kmeans"
#' @param nc.index: Index for determining optimal cluster number in NbClust, default "all" (uses 26 clustering indices)
#' @param cluster: Base clustering method, default "kmeans". Currently also supports "pam"
#' 
#' @return list: 1) optimalK: optimal cluster number; 2) CAM: co-association matrix (binary matrix where 1 = same cluster, 0 = different cluster); 3) i_clus: cluster information for cell samples
####
ConstructCAM <- function(erMetrix, nc.distance = "euclidean", nc.min = 2, nc.max = 8, nc.method = "kmeans", nc.index = "all",
                         cluster = "kmeans"){
  library(dplyr)
  library(NbClust)
  library(purrr)
  library(Seurat)
  
  #####
  # 0 Preprocess expression matrix
  ####
  # 0.1 Filter non-expressing columns
  erMetrix2 <- erMetrix[, colSums(abs(erMetrix)) != 0]
  # 0.2 Filter non-expressing genes
  erMetrix2 <- erMetrix2[rowSums(abs(erMetrix2)) != 0, ]
  # 0.3 Format conversion
  erMetrix2 <- t(erMetrix2) %>% as.data.frame
  
  #####
  # 1 Generate base clustering
  ####
  # 1.1 Calculate clustering indices for different k values
  set.seed(1234) # Seed for reproducibility since kmeans is stochastic
  nb_clust <- NbClust(erMetrix2, distance = nc.distance, min.nc = nc.min, max.nc = nc.max, method = nc.method,
                      index = nc.index)
  # 1.2 Determine optimal k using voting method
  t.clust <- table(nb_clust$Best.nc[1,])
  optimalK <- names(t.clust)[which.max(t.clust)] %>% as.numeric
  # 1.3 Perform base clustering
  baseCluster <- switch(cluster,
                        kmeans = kmeans(erMetrix2, optimalK),
                        pam = pam(erMetrix2, optimalK))
  i_clus <- baseCluster$cluster # Base clustering information
  
  #####
  # 2 Build co-association matrix from base clustering results
  ####
  # 2.1 Construct co-association matrix
  df.clust <- as.data.frame(i_clus) %>%
    tibble::rownames_to_column("cell") # Convert row names to column named "cell"
  outer_join <- merge(df.clust, df.clust, by = 'i_clus', all = TRUE)
  CAM <- table(outer_join$cell.x, outer_join$cell.y) %>%
    as.data.frame.array # Convert to data frame
  # 2.2 Maintain original sample count (pad with 0 for removed samples)
  t.filter <- setdiff(colnames(erMetrix), rownames(CAM)) # Removed samples
  # Pad filtered samples with 0
  CAM[t.filter, t.filter] <- 0
  CAM[is.na(CAM)] <- 0
  # Generate final result
  CAM <- as.matrix(CAM)
  diag(CAM) <- 1 # Set diagonal to 1
  CAM <- CAM[colnames(erMetrix), colnames(erMetrix)] # Sort samples
  return(list(optimalK = optimalK, CAM = CAM, i_clus = i_clus))
}

#####
# Function 2: Improve co-association matrix through two-level refinement at sample and cluster levels
#'
#' @description
#' Performs two-level weighting on the co-association matrix at point (sample) and cluster levels.
#' 1) First weighting at sample level
#' 2) Second weighting at cluster level
#' 
#' @param results.CAM: ConstructCAM results including co-association matrix, cluster information, and base clustering results
#' @param malignantObj: Seurat object for malignant cells
#' @param nfeatures: Number of highly variable genes for dimensionality reduction, default NULL (no reduction)
#' @param erMatrix: Gene expression matrix (genes as rows, samples as columns)
#' @param distance: Distance measure for sample weighting, default "euclidean". Options: "kendall", "manhattan", "spearman", "pearson"
#' 
#' @return Optimized co-association matrix
####
TwoLevelRefine <- function(results.CAM, malignantObj, nfeatures = NULL, erMatrix, distance = "euclidean"){
  #####
  # I Generate initial co-association matrix
  ####
  # Sum CAM results to obtain global co-association matrix
  CAM <- matrix(0, ncol(erMatrix), ncol(erMatrix))
  for(i in 1:length(results.CAM)){
    CAM <- CAM + results.CAM[[i]]$CAM
  }
  CAM <- CAM / length(results.CAM)
 
  #####
  # II First weighting at sample level
  ####
  # 1 Calculate weights using highly variable genes
  if(is.null(nfeatures)){
    weight_matrix <- as.matrix(erMatrix)
  } else {
    features <- FindVariableFeatures(malignantObj, selection.method = "vst", nfeatures = nfeatures)
    weight_matrix <- erMatrix[VariableFeatures(features),] %>% as.matrix()
  }

  weight <- switch(distance,
                   euclidean = dist(t(weight_matrix), method = "euclidean") %>% as.matrix(),
                   manhattan = dist(t(weight_matrix), method = "manhattan") %>% as.matrix(),
                   kendall = cor(weight_matrix, method = "kendall"),
                   spearman = cor(weight_matrix, method = "spearman"),
                   pearson = cor(weight_matrix, method = "pearson")) # Fixed typo: person -> pearson
  W <- 1 - (weight / max(weight)) # Normalize to [0,1] range
  # 2 Apply weights to co-association matrix
  first_matrix <- CAM * W

  #####
  # III Second weighting at cluster level
  ####
  # 1 Integrate base cluster information
  CLS <- lapply(results.CAM, function(x){
    cluster <- as.data.frame(x$i_clus)
    colnames(cluster) <- "cluster"
    cluster <- tibble::rownames_to_column(cluster, "cellID")
    # Calculate intra-cluster similarity weight
    t.cls <- NULL
    for(i in 1:max(x$i_clus)) {
      Num <- sum(cluster$cluster == i)
      t.cls[i] <- 1 / (Num * (Num - 1) / 2) # Inverse of pairwise combinations
    }
    return(list(cls = t.cls))
  })
  
  # 2 Find global min and max weights
  MAX <- max(unlist(CLS))
  MIN <- min(unlist(CLS))
  
  # 3 Standardize weights
  End_matrix <- matrix(0, ncol(erMatrix), ncol(erMatrix))
  for (i in 1:length(results.CAM)) {
    scal <- (map(CLS, 1)[[i]] - MIN) / (MAX - MIN) # Min-max normalization
    scal <- as.data.frame(scal)
    weight <- tibble::rownames_to_column(scal, "cluster")
    cluster <- as.data.frame(map(results.CAM, 3)[[i]])
    colnames(cluster) <- "cluster"
    cluster <- tibble::rownames_to_column(cluster, "cellID")
    
    # 4 Apply cluster-level weights to matrix
    L <- join(cluster, weight, by = "cluster")
    sort <- as.data.frame(colnames(erMatrix))
    colnames(sort) <- "cellID"
    SORT <- join(sort, L, by = "cellID")
    # 5 Replace NA with 0
    SORT <- mutate_all(SORT, ~replace(., is.na(.), 0))
    SORT <- SORT[, -c(1:2)] # Keep only weights
    
    # 6 Apply cluster weights to the matrix
    end_matrix <- sweep(first_matrix, MARGIN = 2, SORT, `*`)
    End_matrix <- End_matrix + end_matrix
  }
  
  # Final weighted matrix
  # End_matrix <- End_matrix / length(results.CAM) # Optional normalization
  End_matrix[lower.tri(End_matrix, diag = TRUE)] <- 0 # Zero lower triangle
  End_matrix <- End_matrix + t(End_matrix) # Make symmetric
  return(End_matrix)
}