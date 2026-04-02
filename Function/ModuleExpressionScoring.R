#########################################################
# TODO: Calculate module expression scores for each cell
#
# Main processes
#  1) Input modules and generate 1000 random gene sets with similar expression levels
#  2) For each cell, calculate the mean centered expression of the module and control gene sets
#  3) Compute the proportion of control gene sets whose mean expression exceeds that of the module (defined as p)
#  4) Define the score as -log10(p), then linearly rescale to [0,1]
# Ref: Cancer cell states recur across tumor types and form specific interactions with the tumor microenvironment, Nature Genetics, 2022
#
#' @param moduleList List of gene sets (modules), each element is a vector of genes, named by module
#' @param seuratObj Seurat object
#' @param mc.cores Number of cores for parallel computation (default = 1, no parallelization)
#' @param bin_num Number of bins used to group genes (default = 25)
#' @param assay Expression data used for scoring ("SCT" or "RNA")
#
#' @return A list of module scores, each element is a vector of scores for one module
#
#' @note After rescaling, scores are comparable across datasets and can be used to define activation states (score > 0.5 indicates activation)
#
############################################################################

ModuleExpressionScoring <- function(moduleList, seuratObj, mc.cores = 1, bin_num = 25, assay = "SCT"){
  library(Seurat)
  library(ggplot2)
  library(parallel)
  set.seed(123)
  
  #####
  # 0 Data preprocessing
  ####
  # 0.1 Background gene preprocessing
  # Sort genes based on mean expression
  if (assay == "SCT") {
    geneSorted <- sort(rowMeans(as.matrix(seuratObj@assays$SCT@data)))
  } else if (assay == "RNA") {
    geneSorted <- sort(rowMeans(as.matrix(seuratObj@assays$RNA@data)))
  } else {
    stop("Invalid assay type. Please choose 'SCT' or 'RNA'.")
  }
  
  # Divide all genes into bin_num bins
  data.cut <- cut_number(x = geneSorted + rnorm(n = length(geneSorted))/1e+15, 
                         n = bin_num, labels = FALSE, right = FALSE)
  names(data.cut) <- names(geneSorted)
  
  # Generate gene bins (list of genes grouped by bins)
  geneBins <- split(names(data.cut), data.cut)
  
  # 0.2 Centered expression matrix (for all genes)
  seuratObj.centered <- ScaleData(seuratObj, features = rownames(seuratObj))
  
  # 0.3 Reconstruct module information
  # Intersect with background genes
  moduleList <- lapply(moduleList, function(x) intersect(x, names(geneSorted)))
  
  # Check module names; assign numeric names if missing
  if(is.null(names(moduleList))) names(moduleList) <- 1:length(moduleList)
  
  # Store module information as a data frame
  ## geneset: genes in each module
  ## name: module name
  ## setIndex: module index
  ## bin: bin assignment for each gene
  module.df <- data.frame(
    geneset = unlist(moduleList), 
    name = rep(names(moduleList), sapply(moduleList, length)),
    setIndex = rep(1:length(moduleList), sapply(moduleList, length))
  )
  rownames(module.df) <- NULL
  module.df$bin <- data.cut[module.df$geneset]
  
  #####
  # 1 Generate control gene sets for each module
  ####
  controlGenes <- sapply(1:nrow(module.df), function(i){
    # Extract all genes in the corresponding module
    t.geneset <- moduleList[[module.df[i,"setIndex"]]]
    
    # Candidate genes from the same bin
    t1 <- setdiff(geneBins[[module.df[i,"bin"]]], t.geneset) # Avoid selecting module genes as controls
    
    # Sample 1000 random genes (with replacement)
    return(sample(t1, 1000, replace = TRUE))
  })
  
  # Transpose: each column represents one random sampling
  controlGenes <- t(controlGenes)
  
  #####
  # 2 Compute module scores
  ####
  moduleScoreList <- sapply(1:length(moduleList), function(i){
    cat("Calculating: ", names(moduleList)[i], "\n")
    module <- moduleList[[i]]
    t.index <- module.df$setIndex == i # Indices of module genes in module.df
    
    # 2.1 Compute mean expression of module and control gene sets
    if (assay == "SCT") {
      exprMean.module <- colMeans(seuratObj.centered@assays$SCT@scale.data[module.df$geneset[t.index], ])
      exprMean.control <- mclapply(1:1000, function(j){
        colMeans(seuratObj.centered@assays$SCT@scale.data[controlGenes[t.index, j], ])
      }, mc.cores = mc.cores)
    } else if (assay == "RNA") {
      exprMean.module <- colMeans(seuratObj.centered@assays$RNA@scale.data[module.df$geneset[t.index], ])
      exprMean.control <- mclapply(1:1000, function(j){
        colMeans(seuratObj.centered@assays$RNA@scale.data[controlGenes[t.index, j], ])
      }, mc.cores = mc.cores)
    }
    exprMean.control <- do.call(cbind, exprMean.control)
    
    # 2.2 Compute the proportion of control sets exceeding module expression
    t.matrix <- (exprMean.control - exprMean.module) >= 0
    p <- rowSums(t.matrix)/1000
    
    # 2.3 Define score as -log10(p) and normalize
    moduleScore <- -log10(p)
    moduleScore[is.infinite(moduleScore)] <- -log10(0.001) # Replace Inf with upper bound (3)
    moduleScore <- moduleScore / (-log10(0.001))
    
    # Return scores
    return(moduleScore)
  })
  
  colnames(moduleScoreList) = names(moduleList)
  rownames(moduleScoreList) = seuratObj$CellID
  cat("Done!", "\n")
  
  return(moduleScoreList)
}
