LR_tidymodels <- function(data, class, features = ".", train.p = 3/4, strata = NULL, cv = 10,
                           method = "logistic_reg", normalize = FALSE, tune.param = TRUE){
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
    lr.rec <- data.train %>%
      recipe(formula = t.formula) %>% 
      step_center(all_predictors(), -all_outcomes()) %>% 
      step_scale(all_predictors(), -all_outcomes())
  } else {
    lr.rec <- data.train %>%
      recipe(formula = t.formula)
  }
  
  # 2.3 Define model specification
  if (tune.param){
    lr.model <- switch(method,
                        logistic_reg = logistic_reg(penalty = tune(), mixture = tune()) %>% 
                          set_engine("glmnet") %>%
                          set_mode("classification"))
  } else {
    lr.model <- switch(method,
                        logistic_reg = logistic_reg(penalty = double(1), mixture = double(1)) %>% 
                          set_engine("glmnet") %>%
                          set_mode("classification"))
  }
  
  # 2.4 Create workflow
  lr.wf <- workflow() %>% 
    add_model(lr.model) %>% 
    add_recipe(lr.rec)
  
  #####
  # 3 Set up cross-validation
  ####
  set.seed(234)
  lr.cv <- vfold_cv(data.train, v = cv, strata = all_of(strata))
  
  #####
  # 4 Train model
  ####
  if (tune.param){ # Tune hyperparameters and select best model
    # Create parameter grid
    lr.grid <- grid_regular(
      penalty(range = c(-3, 0)), mixture(range = c(0.01, 0.99)), levels = c(20, 5))
    
    # Evaluate model performance across parameters
    lr.res <- lr.wf %>% 
      tune_grid(
        resamples = lr.cv,
        grid = lr.grid
      )
    
    # Select best model
    lr.res <- select_best(lr.res) %>% 
      finalize_workflow(lr.wf, .)
  } else {
    # Fit model without tuning
    lr.res <- 
      lr.wf %>% 
      fit_resamples(lr.cv, control = control_resamples(save_pred = TRUE))
  }
  
  #####
  # 5 Return prediction results
  if (tune.param){
    # If using tuned parameters, predict on test set
    test.res <- lr.res %>% 
      fit(data = data.test)
    
    # Format prediction results
    pred.res <- predict(test.res, data.test) %>% 
      bind_cols(dplyr::select(data.test, all_of(class))) %>% 
      bind_cols(predict(test.res, data.test, type = "prob"))
  } else {
    # Extract predictions from resampling
    pred.res <- lr.res %>% 
      collect_predictions()
  }
  
  return(pred.res)
}