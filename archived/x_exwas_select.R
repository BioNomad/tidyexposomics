exwas_select <- function(
    expOmicSet, 
    exposures, 
    outcome, 
    confounders = NULL, 
    family = "gaussian", 
    cv_folds = 10, 
    penalty_range = c(-6, 0), 
    mixture_range = c(0, 1), 
    penalty_levels = 100, 
    mixture_levels = 10
) {
  library(tidymodels)
  library(tidyverse)
  library(doParallel)
  
  # extract and preprocess coldata
  message("Extracting exposure data...")
  data <- colData(expOmicSet) |>
    as.data.frame()
  
  # combine exposures, confounders, and outcome
  predictors <- c(exposures, confounders)
  data <- data |>
    dplyr::select(all_of(c(predictors, outcome)))
  
  # remove rows with missing data
  message("Removing rows with missing data...")
  data <- data[complete.cases(data), ]
  
  # remove zero-variance predictors
  message("Removing zero-variance predictors...")
  nzv <- caret::nearZeroVar(data, saveMetrics = TRUE)
  if (any(nzv$nzv)) {
    removed_vars <- rownames(nzv[nzv$nzv, ])
    message(" - removed zero-variance predictors: ", paste(removed_vars, collapse = ", "))
    data <- data |>
      dplyr::select(-all_of(removed_vars))
    predictors <- setdiff(predictors, removed_vars)  # update predictors list
  }
  
  # scale numeric variables and convert character variables to numeric
  message("Scaling numeric variables and converting character variables to numeric...")
  data <- data |>
    mutate(across(where(is.numeric), ~ as.numeric(scale(.)))) |>
    mutate(across(where(is.character), ~ as.numeric(factor(.))))
  
  # prepare formula
  message("Preparing formula...")
  predictor_formula <- paste(predictors, collapse = " + ")
  formula <- as.formula(paste(outcome, "~", predictor_formula))
  
  # create recipe
  message("Creating recipe...")
  recipe <- recipe(formula, data = data) |>
    step_normalize(all_predictors())
  
  # set up cross-validation
  message("Setting up cross-validation (", cv_folds, " folds)...")
  validation_split <- vfold_cv(data, v = cv_folds, strata = outcome)
  
  # define elastic net model
  message("Defining elastic net model...")
  tune_spec <- linear_reg(penalty = tune(), mixture = tune()) |>
    set_engine("glmnet") |>
    set_mode(ifelse(family == "gaussian", "regression", "classification"))
  
  # define parameter grid
  message("Defining parameter grid...")
  param_grid <- grid_regular(
    penalty(range = penalty_range), 
    mixture(range = mixture_range), 
    levels = c(penalty = penalty_levels, mixture = mixture_levels)
  )
  
  # create workflow
  message("Creating workflow...")
  workflow <- workflow() |>
    add_recipe(recipe) |>
    add_model(tune_spec)
  
  # tune elastic net model
  message("Tuning elastic net model...")
  doParallel::registerDoParallel()
  tune_result <- workflow |>
    tune_grid(
      resamples = validation_split,
      grid = param_grid,
      metrics = metric_set(rmse)
    )
  
  # select best model
  message("Selecting best elastic net model...")
  tune_best <- tune_result |>
    select_best(metric = "rmse")
  
  message("Best parameters: penalty = ", tune_best$penalty, ", mixture = ", tune_best$mixture)
  
  # fit final model
  message("Fitting final elastic net model...")
  elastic_model <- linear_reg(penalty = tune_best$penalty, mixture = tune_best$mixture) |>
    set_engine("glmnet") |>
    set_mode(ifelse(family == "gaussian", "regression", "classification"))
  
  elastic_fit <- elastic_model |>
    fit(formula, data = data)
  
  # extract coefficients
  message("Extracting coefficients...")
  coefficients <- tidy(elastic_fit) |>
    filter(estimate != 0)  # retain non-zero coefficients
  
  # save results to metadata
  message("Saving results to expOmicSet metadata...")
  metadata(expOmicSet)$variable_selection <- list(
    best_params = tune_best,
    coefficients = coefficients,
    tune_result = tune_result,
    elastic_fit = elastic_fit
  )
  
  message("Variable selection exwas completed.")
  return(expOmicSet)
}
