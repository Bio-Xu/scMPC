###################################################
# TODO: Feature gene selection using caret package with custom RF function
#
#' @detail Custom RF function is required to use Precision as the evaluation metric
#' 
#' @param x: Data frame containing gene expression profiles (features as columns, samples as rows)
#' @param y: Class label variable name in the data
#' @param size: Vector of feature subset sizes to evaluate
#' @param metric: Metric for selecting best model, default Precision (custom)
#' @param method: Resampling method, default "cv" (cross-validation)
#' @param number: Number of folds for cross-validation, default 10
#' @param repeats: For repeated k-fold CV: number of complete fold sets to compute, default 1
#' @param functions: Predefined caret functions that can be directly called:
#'                   - Linear regression: lmFuncs
#'                   - Random forest: rfFuncs
#'                   - Naive Bayes: nbFuncs
#'                   - Bootstrap aggregating: treebagFuncs
#' @param multiClassSummary: Function that includes overall measures, but without class probabilities (classProbs),
#'                           logloss and AUC will not appear in results
#' @param importance: When all predictors are in the model, feature importance doesn't need to be recalculated each iteration.
#'                   Set importance = "first" to use only the first importance ranking
#' 
#' @return RFE results object containing feature selection metrics and optimal features
######################################################################################################

RF_RFE <- function(x, y, size, metric = "Precision", method = "cv", number = 10) {
  library(randomForest)
  library(mlbench)
  library(caret)
  
  ## Custom RF functions (modified from rfFuncs to use Precision metric)
  rfRFE <- list(
    ### Calculate performance metrics from true vs predicted values
    summary = multiClassSummary,
    
    ### Model fitting function
    fit = function(x, y, first, last, ...) {
      library(randomForest)
      randomForest(x, y, importance = first, ...)
    },
    
    ### Prediction function - ensures factor levels match input data
    pred = function(object, x) predict(object, x),
    
    ### Feature importance ranking function
    rank = function(object, x, y) {
      vimp <- varImp(object)
      if (is.factor(y)) {
        if (all(levels(y) %in% colnames(vimp))) {
          avImp <- apply(vimp[, levels(y), drop = TRUE], 1, mean)
          vimp$Overall <- avImp
        }
      }
      vimp <- vimp[order(vimp$Overall, decreasing = TRUE), , drop = FALSE]
      if (ncol(x) == 1) {
        vimp$var <- colnames(x)
      } else {
        vimp$var <- rownames(vimp)
      }
      vimp
    },
    
    ### Determine optimal feature set size based on resampling
    selectSize = pickSizeBest,
    
    ### Identify best features across all resampling iterations
    selectVar = pickVars
  )

  # Set up RFE control parameters
  control <- rfeControl(
    functions = rfRFE, 
    method = method, 
    number = number
  )
  
  # Run the RFE algorithm
  results <- rfe(
    x = x, 
    y = y, 
    sizes = size, 
    rfeControl = control,
    metric = metric
  )
  
  # Print summary of results
  print(results)
  
  return(results)
}