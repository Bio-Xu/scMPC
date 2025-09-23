##################################
# TODO: Evaluate features using sliding window structure with three classifiers
#
#' @param expr: Gene expression matrix. First column contains gene names, subsequent columns contain sample expressions.
#' @param pData: Phenotype label information corresponding to expression matrix samples. First column: sample ID; Second column: phenotype labels (two types, e.g., metastatic-nonmetastatic)
#' @param genelist: Rank-integrated gene set sorted by importance
#' @details
#' Uses three classifiers (logistic regression; random forest; support vector machine) and two classification metrics (precision, kappa)
#' SVM uses RBF kernel based on test results showing better performance than linear kernel
##############
Classifier <- function(expr, pData, genelist){
  
  library(dplyr)
  
  # Source classifier functions
  source("Function\\SVM_tidymodels.R")
  source("Function\\RF_tidymodels.R")
  source("Function\\LR_tidymodels.R")
  
  ## Prepare expression matrix based on phenotype data
  rownames(expr) <- expr[,1]
  expr <- expr[,-1] 
  # Match samples between expression and phenotype data
  sam <- intersect(pData[,1], colnames(expr))
  expr <- expr[,sam]
  expr <- as.data.frame(t(expr))
  pData <- dplyr::select(pData, -3) 
  # Combine phenotype information with expression matrix
  expr <- cbind(pData, expr)
  expr <- dplyr::select(expr, -1) 
  expr$group <- factor(expr$group)
  
  ## Initialize data frame to store classifier evaluation results
  SumResult <- data.frame()
  
  ## Sliding window over signature genes - add one gene at a time starting from top 2
  for(i in 2:length(genelist)){ 
    print(i)
    # Select top i genes from the ranked list
    sam1 <- intersect(genelist[1:i], colnames(expr))
    expr1 <- cbind(expr[,1], expr[,sam1])
    colnames(expr1)[1] <- "group"
    
    # Train classifiers (without hyperparameter tuning)
    # SVM with RBF kernel
    svm_r <- SVM_tidymodels(data = expr1, class = "group", features = ".", train.p = 3/4, 
                           strata = "group", cv = 10, method = "rbf", normalize = FALSE, 
                           tune.param = FALSE)
    
    # Random Forest
    rf <- RF_tidymodels(data = expr1, class = "group", features = ".", train.p = 3/4, 
                       strata = "group", cv = 10, method = "rand_forest", normalize = FALSE, 
                       tune.param = FALSE)
    
    # Logistic Regression
    lr <- LR_tidymodels(data = expr1, class = "group", features = ".", train.p = 3/4, 
                       strata = "group", cv = 10, method = "logistic_reg", normalize = FALSE, 
                       tune.param = FALSE)
    
    ## Performance metric calculation function
    performance <- function(Matrix){
      cfMatrix <- table(Matrix$.pred_class, Matrix$group)
      tp <- cfMatrix[1,1]
      fp <- cfMatrix[1,2]
      tn <- cfMatrix[2,2]
      fn <- cfMatrix[2,1]
      N <- tp + fn + fp + tn
      precision <- tp/(tp + fp)
      res.kappa <- vcd::Kappa(cfMatrix)
      kappa <- res.kappa[[1]][[1]]
      result <- list(precision = precision, kappa = kappa)
      return(result)
    }   
    
    ## Get evaluation results for different classifiers
    rf_result <- performance(rf)
    svm_result <- performance(svm_r)
    lr_result <- performance(lr)
    
    # Combine results for current gene set size
    Result <- data.frame(RF_precision = rf_result$precision, 
                         RF_kappa = rf_result$kappa,
                         SVM_precision = svm_result$precision, 
                         SVM_kappa = svm_result$kappa,
                         LR_precision = lr_result$precision, 
                         LR_kappa = lr_result$kappa)
    
    SumResult <- rbind(SumResult, Result)
  }
  return(SumResult)
}