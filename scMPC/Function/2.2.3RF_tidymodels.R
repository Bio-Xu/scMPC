RF_tidymodels <- function(data, class, features = ".", train.p = 3/4, strata = NULL, cv = 10,
                          method = "rand_forest", normalize = FALSE, tune.param = TRUE){
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
  set.seed(123)
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
    rf.rec <- data.train %>%
      recipe(formula = t.formula) %>% 
      step_center(all_predictors(), -all_outcomes()) %>% 
      step_scale(all_predictors(), -all_outcomes())
  } else {
    rf.rec <- data.train %>%
      recipe(formula = t.formula)
  }
  
  # 2.3 Define random forest model specification
  if (tune.param){
    rf.model <- switch(method,
                       rand_forest = rand_forest(mtry = tune(), trees = tune()) %>% 
                         set_engine("ranger") %>%
                         set_mode("classification"))
  } else {
    rf.model <- switch(method,
                       rand_forest = rand_forest() %>% 
                         set_engine("ranger") %>%
                         set_mode("classification"))
  }
  
  # 2.4 Create workflow
  rf.wf <- workflow() %>% 
    add_model(rf.model) %>% 
    add_recipe(rf.rec)
  
  #####
  # 3 Set up cross-validation
  ####
  set.seed(234)
  rf.cv <- vfold_cv(data.train, v = cv, strata = all_of(strata))
  
  #####
  # 4 Train model
  ####
  if (tune.param){ # Tune hyperparameters and select best model
    # Create parameter grid
    rf.grid <- grid_regular(
      mtry(range = c(5, 10)), trees(range = c(500, 1000)), levels = c(5, 6))
    
    # Evaluate model performance across parameters
    rf.res <- rf.wf %>% 
      tune_grid(
        resamples = rf.cv,
        grid = rf.grid
      )
    
    # Select best model
    rf.res <- select_best(rf.res) %>% 
      finalize_workflow(rf.wf, .)
  } else {
    # Fit model without tuning
    rf.res <- 
      rf.wf %>% 
      fit_resamples(rf.cv, control = control_resamples(save_pred = TRUE))
  }
  
  #####
  # 5 Return prediction results
  if (tune.param){
    # If using tuned parameters, predict on test set
    test.res <- rf.res %>% 
      fit(data = data.test)
    
    # Format prediction results
    pred.res <- predict(test.res, data.test) %>% 
      bind_cols(dplyr::select(data.test, all_of(class))) %>% 
      bind_cols(predict(test.res, data.test, type = "prob"))
  } else {
    # Extract predictions from resampling
    pred.res <- rf.res %>% 
      collect_predictions()
  }
  
  return(pred.res)
}