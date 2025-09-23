SVM_tidymodels <- function(data, class, features = ".", train.p = 3/4, strata = NULL, cv = 10,
                           method = c("linear", "rbf"), normalize = FALSE, tune.param = TRUE){
  library(tidyverse)
  library(tidymodels)
  
  # Convert character variables to factors
  if (identical(features, ".")){
    data.model <- data %>%
      mutate_if(is.character, as.factor)
  } else {
    data.model <- data %>%
      select(all_of(c(class, features, strata))) %>% 
      mutate_if(is.character, as.factor)
  }
  
  #####
  # 1 Split data into training and testing sets
  ####
  data.split <- initial_split(data.model, prop = train.p, strata = all_of(strata))
  # Training set
  data.train <- training(data.split)
  # Testing set
  data.test <- testing(data.split)
  
  #####
  # 2 Create modeling workflow
  ####
  # 2.1 Define formula
  t.formula <- as.formula(paste(class, "~", paste(features, collapse = "+")))
  
  # 2.2 Create recipe with optional normalization
  if (normalize){
    svm.rec <- data.train %>%
      recipe(formula = t.formula) %>% 
      step_center(all_predictors(), -all_outcomes()) %>% 
      step_scale(all_predictors(), -all_outcomes())
  } else {
    svm.rec <- data.train %>%
      recipe(formula = t.formula)
  }
  
  # 2.3 Define SVM model specification
  if (tune.param){
    svm.model <- switch(method,
                        linear = svm_linear(cost = tune()) %>% 
                          set_engine("kernlab") %>%
                          set_mode("classification"),
                        rbf = svm_rbf(cost = tune(), rbf_sigma = tune()) %>%
                          set_engine("kernlab") %>%
                          set_mode("classification"))
  } else {
    svm.model <- switch(method,
                        linear = svm_linear() %>% 
                          set_engine("kernlab") %>%
                          set_mode("classification"),
                        rbf = svm_rbf() %>%
                          set_engine("kernlab") %>%
                          set_mode("classification"))
  }
  
  # 2.4 Create workflow
  svm.wf <- workflow() %>% 
    add_model(svm.model) %>% 
    add_recipe(svm.rec)
  
  #####
  # 3 Set up cross-validation
  ####
  set.seed(234)
  svm.cv <- vfold_cv(data.train, v = cv, strata = all_of(strata))
  
  #####
  # 4 Train model
  ####
  if (tune.param){ # Tune hyperparameters and select best model
    # Create parameter grid
    svm.grid <- grid_regular(
      extract_parameter_set_dials(svm.model), levels = 5)
    
    # Evaluate model performance across parameters
    svm.res <- svm.wf %>% 
      tune_grid(
        resamples = svm.cv,
        grid = svm.grid
      )
    
    # Select best model
    svm.res <- select_best(svm.res) %>% 
      finalize_workflow(svm.wf, .)
  } else {
    # Fit model without tuning
    svm.res <- 
      svm.wf %>% 
      fit_resamples(svm.cv, control = control_resamples(save_pred = TRUE))
  }
  
  #####
  # 5 Return prediction results
  if (tune.param){
    # If using tuned parameters, predict on test set
    test.res <- svm.res %>% 
      fit(data = data.test)
    
    # Format prediction results
    pred.res <- predict(test.res, data.test) %>% 
      bind_cols(select(data.test, all_of(class))) %>% 
      bind_cols(predict(test.res, data.test, type = "prob"))
  } else {
    # Extract predictions from resampling
    pred.res <- svm.res %>% 
      collect_predictions()
  }
  
  return(pred.res)
}