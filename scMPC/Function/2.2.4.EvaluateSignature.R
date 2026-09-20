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
    performance <- function(Matrix) {
      Matrix$group <- factor(Matrix$group, levels = c("PT", "LNM"))
      Matrix$.pred_class <- factor(Matrix$.pred_class, levels = c("PT", "LNM"))

      cfMatrix <- table(Matrix$.pred_class, Matrix$group)
      
      if (!all(dim(cfMatrix) == c(2,2))) {
        stop("Confusion matrix is not 2x2. Check class distribution.")
      }
      
      tp <- cfMatrix["LNM","LNM"]
      fp <- cfMatrix["LNM","PT"]
      fn <- cfMatrix["PT","LNM"]
      tn <- cfMatrix["PT","PT"]
      N <- tp + fp + fn + tn
      
      accuracy <- (tp + tn) / N
      recall <- ifelse((tp + fn) == 0, NA, tp / (tp + fn))   # sensitivity
      precision <- ifelse((tp + fp) == 0, NA, tp / (tp + fp))
      F1 <- ifelse(is.na(precision) | is.na(recall) | (precision + recall) == 0,
                   NA,
                   2 * precision * recall / (precision + recall))
      specificity <- ifelse((tn + fp) == 0, NA, tn / (tn + fp))
      balanced_accuracy <- mean(c(recall, specificity), na.rm = TRUE)
      
      # -------- MCC --------
      mcc_denom <- sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
      MCC <- ifelse(mcc_denom == 0, NA,
                    (tp*tn - fp*fn) / mcc_denom)
      
      # -------- Kappa --------
      res.kappa <- vcd::Kappa(cfMatrix)
      kappa <- res.kappa[[1]][[1]]
      
      # -------- AUC --------
      if (!".pred_LNM" %in% colnames(Matrix)) {
        AUC <- NA
      } else {
        roc_obj <- pROC::roc(
          response = Matrix$group,
          predictor = Matrix$.pred_LNM,
          levels = c("PT", "LNM")
        )
        AUC <- as.numeric(pROC::auc(roc_obj))
      }
      
      result <- list(
        accuracy = accuracy,
        balanced_accuracy = balanced_accuracy,
        precision = precision,
        recall = recall,
        F1 = F1,
        MCC = MCC,
        kappa = kappa,
        AUC = AUC
      )
      return(result)
    }   
    
    ## Get evaluation results for different classifiers
    rf_result <- performance(rf)
    svm_result <- performance(svm_r)
    lr_result <- performance(lr)
    
    # Combine results for current gene set size
    Result <- data.frame(RF_precision = rf_result$precision, 
                         RF_kappa = rf_result$kappa,
                         RF_accuracy = rf_result$accuracy,
                         RF_recall = rf_result$recall,
                         RF_F1 = rf_result$F1,
                         RF_balanced_accuracy  = rf_result$balanced_accuracy ,
                         RF_MCC = rf_result$MCC,
                         RF_AUC = rf_result$AUC,
                         
                         SVM_precision = svm_result$precision, 
                         SVM_kappa = svm_result$kappa,
                         SVM_accuracy = svm_result$accuracy,
                         SVM_recall = svm_result$recall,
                         SVM_F1 = svm_result$F1,
                         SVM_balanced_accuracy  = svm_result$balanced_accuracy ,
                         SVM_MCC = svm_result$MCC,
                         SVM_AUC = svm_result$AUC,
                         
                         LR_precision = lr_result$precision, 
                         LR_kappa = lr_result$kappa),
                         LR_accuracy = lr_result$accuracy,
                         LR_recall = lr_result$recall,
                         LR_F1 = lr_result$F1,
                         LR_balanced_accuracy  = lr_result$balanced_accuracy ,
                         LR_MCC = lr_result$MCC,
                         LR_AUC = lr_result$AUC
    
    SumResult <- rbind(SumResult, Result)
  }
  return(SumResult)
}