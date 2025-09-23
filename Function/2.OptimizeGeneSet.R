###############################################
# TODO: Optimize gene sets
#
#' @description: Optimizes a gene set by removing genes can't distinguish metastasis
# input
#' @param geneSet: Gene set to be optimized.
#' @param erMatrix: Gene expression matrix. First column contains gene names, subsequent columns contain sample expressions.
#' @param pData: Phenotype label information corresponding to expression matrix samples. First column: sample ID; Second column: phenotype labels (two types, e.g., metastatic-nonmetastatic); Third column: patient ID for sample pairing.
#' @param geneNum: Number of genes to randomly select. Default: 30 genes.
#' @param sampleNum: Percentage of samples to randomly select. Default: 70% of total samples.
#' @param randomG: Number of gene randomizations. Default: 1,000,000 times.
#' @param randomS: Number of sample randomizations. Default: 36 times.
#' @param prob: Probability threshold for significant gene sets. Default: 0.9.
#
# output
#' @export optimalGS: Optimized gene set
###########################
OptimizeGeneSet <- function(geneSet, erMatrix, pData, geneNum = 30, sampleNum = 0.7,
                            randomG = 1000000, randomS = 36, prob=0.9){
  # Load required packages
  require(vcd)
  require(Matrix)
  require(dplyr)
  require(mclust)
  require(foreach)
  require(doParallel)
  
  # Create and register parallel cluster
  cl <- makeCluster(4)
  registerDoParallel(cl)
  
  ####
  # STEP 1: Preprocess expression matrix and parameters
  ####
  # 1.1 Extract expressed genes from gene set
  # Intersect gene set with expression matrix genes
  geneSet2 <- unique(intersect(geneSet, erMatrix[,1]))
  # Terminate if expressed genes are fewer than geneNum
  if(length(geneSet2) < geneNum) stop("Expressed gene set smaller than geneNum\n")
  # 1.2 Extract expression matrix for gene set
  erMatrix2 <- erMatrix[match(geneSet2, erMatrix[,1]),]
  
  ####
  # STEP 2: Build score matrix to evaluate gene-phenotype association
  ####
  # 2.1 Initialize variables
  # Score matrix
  scoreMatrix <- Matrix(0, nrow = randomS, ncol = randomG, sparse = TRUE)
  # Gene score vector
  geneScore <- rep(0, length(geneSet2))
  names(geneScore) <- geneSet2
  # Generate random gene sets
  randomGeneList <- foreach(i = 1:randomG) %dopar% {
    sample(geneSet2, geneNum)
  }
  
  # 2.2 Evaluate gene-phenotype discrimination capability
  # 2.2.1 Build randomized expression matrices
  samplMatrix <- lapply(1:randomS, function(x){
    # Randomly sample while maintaining phenotype proportions and pairing
    subData <- split(pData, pData[,2])
    t.sample <- sample(subData[[1]]$Num, nrow(subData[[1]])*sampleNum)
    t.pData <- filter(pData, Num %in% t.sample)
    t.erMatrix <- erMatrix2[,t.pData[,1]]
    rownames(t.erMatrix) <- erMatrix2[,1]
    return(t.erMatrix)
  })
  
  # 2.2.2 Evaluate classification performance of random gene sets
  ScoreMatrix <- foreach(i = 1:randomS, .packages = 'foreach', .combine = rbind) %:%
    foreach(x = randomGeneList, .combine = c) %dopar% {
      # Extract random gene set expression profile
      t.erMatrix2 <- samplMatrix[[i]][x,]
      
      # Fuzzy clustering
      res.cluster <- cluster::fanny(t(t.erMatrix2), k = 2, diss = F)$clustering
      # Count cluster-phenotype correspondences
      clusterTable <- as.data.frame(res.cluster)
      clusterTable[pData[,1],"label"] <- pData[,2]
      clusterTable <- na.omit(clusterTable)
      
      # Align cluster labels with phenotype labels (choose best match)
      t.table1 <- table(clusterTable[,2], clusterTable[,1])
      t.table2 <- table(clusterTable[,2], factor(clusterTable[,1], levels = c(2,1)))
      crossTable <- if(t.table1[1,1] > t.table2[1,1]) t.table1 else t.table2
      
      # Calculate Kappa agreement between clustering and true labels
      res.kappa <- vcd::Kappa(crossTable)
      # Compute significance p-value using Z-test
      value <- res.kappa[[1]][[1]]
      ASE <- res.kappa[[1]][[2]]
      Z <- value/ASE
      t.pvalue <- 2*pnorm(q=Z, lower.tail=F)
      # Record effectiveness of random gene set
      t.value <- ifelse(t.pvalue < 0.05, 1, 0)
      return(t.value)
    }
  
  # Shut down parallel cluster
  stopImplicitCluster()
  stopCluster(cl)
  
  # 2.2.3 Extract significant random gene sets and count genes
  t.count <- apply(ScoreMatrix,2,sum)
  for(x in which(t.count > randomS*prob)){
    t.gene <- randomGeneList[[x]]
    geneScore[t.gene] <- geneScore[t.gene] + 1
  } 
  
  return(sort(geneScore, decreasing = T))
}