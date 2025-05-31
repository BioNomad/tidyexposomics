#' Run Elastic Net with Tidymodels
#'
#' Fits an elastic net model using tidymodels, with automatic handling of classification or regression,
#' hyperparameter tuning via cross-validation, and variable importance (VIP) extraction.
#'
#' @param data A data frame containing both predictors and the outcome variable.
#' @param outcome Character string specifying the name of the outcome variable. Should be a factor (for classification) or numeric (for regression).
#' @param predictors Optional character vector of predictor variable names. If `NULL`, all columns except the outcome are used.
#' @param penalty_grid Integer; number of penalty values to evaluate in the regularization grid. Default is 10.
#' @param n_folds Integer; number of folds for cross-validation. Default is 10.
#' @param seed Integer; random seed for reproducibility. Default is 42.
#' @param mixture Numeric value between 0 and 1; elastic net mixing parameter. `0` is ridge, `1` is lasso. Default is 0.5.
#'
#' @return A list containing:
#' \item{model}{Fitted final model.}
#' \item{best_penalty}{Best penalty value selected via cross-validation.}
#' \item{predictions}{Predicted values on the training data.}
#' \item{confusion}{(Classification only) Confusion matrix.}
#' \item{roc_auc}{(Classification only) ROC AUC score.}
#' \item{rsq}{R-squared or McFadden R^2 depending on outcome type.}
#' \item{rmse}{(Regression only) Root mean squared error.}
#' \item{vip}{Variable importance scores based on absolute model coefficients.}
#'
#' @examples
#' \dontrun{
#' data <- your_dataframe
#' result <- run_enet(data, outcome = "outcome_variable")
#' }
#'
#' @import tidymodels performance
#' @importFrom dplyr select filter mutate arrange all_of bind_cols
#' @importFrom broom tidy
#' @importFrom stats glm as.formula
#' @export

run_enet <- function(
    data,
    outcome,
    predictors = NULL,
    penalty_grid = 10,
    n_folds = 10,
    seed = 42,
    mixture = 0.5 # 0 = ridge, 1 = lasso
) {
  require(tidymodels)
  require(performance)

  set.seed(seed)

  message("Detecting outcome type...")

  outcome_var <- data[[outcome]]

  if (is.ordered(outcome_var)) {
    outcome_var <- factor(outcome_var, ordered = FALSE)
    data[[outcome]] <- outcome_var
  }

  if (is.factor(outcome_var)) {
    model_mode <- "classification"
  } else if (is.numeric(outcome_var)) {
    model_mode <- "regression"
  } else {
    stop("Outcome variable must be either a factor or numeric.")
  }

  outcome_sym <- rlang::sym(outcome)
  if (is.null(predictors)) {
    predictors <- setdiff(names(data), outcome)
  }

  message("Building preprocessing recipe...")

  recipe_obj <- recipes::recipe(
    as.formula(paste(outcome, "~ .")),
    data = data) |>
    recipes::update_role(!!outcome_sym, new_role = "outcome") |>
    recipes::step_zv(recipes::all_predictors()) |>
    recipes::step_novel(recipes::all_nominal_predictors()) |>
    recipes::step_dummy(recipes::all_nominal_predictors()) |>
    recipes::step_zv(recipes::all_predictors()) |>
    recipes::step_normalize(recipes::all_numeric_predictors())

  message("Creating model specification...")

  model_spec <- if (model_mode == "classification") {
    parsnip::logistic_reg(
      penalty = tune(),
      mixture = tune()) |>
      parsnip::set_engine("glmnet") |>
      parsnip::set_mode("classification")

  } else {

    parsnip::linear_reg(
      penalty = tune(),
      mixture = tune()) |>
      parsnip::set_engine("glmnet") |>
      parsnip::set_mode("regression")
  }

  message("Creating workflow...")

  wf <- workflows::workflow() |>
    workflows::add_model(model_spec) |>
    workflows::add_recipe(recipe_obj)

  message("Running ", n_folds, "-fold cross-validation...")

  folds <- rsample::vfold_cv(
    data,
    v = n_folds,
    strata = if (model_mode == "classification") outcome else NULL)

  grid <- dials::grid_regular(
    dials::penalty(),
    dials::mixture(),
    levels = penalty_grid
  )

  message("Tuning hyperparameters...")

  tuned <- tune::tune_grid(
    wf,
    resamples = folds,
    grid = grid,
    metrics = if (model_mode == "classification") {
      yardstick::metric_set(
        yardstick::roc_auc,
        yardstick::accuracy)

    } else {

      yardstick::metric_set(
        yardstick::rmse,
        yardstick::rsq)
    }
  )

  best_metric <- ifelse(
    model_mode == "classification",
    "roc_auc",
    "rmse")

  best_penalty <- tune::select_best(
    tuned,
    metric = best_metric)

  message("Best penalty selected: ", signif(best_penalty$penalty, 3))

  final_wf <- tune::finalize_workflow(
    wf,
    best_penalty)

  message("Fitting final model to all data...")

  final_fit <- parsnip::fit(
    final_wf,
    data = data)

  # Extract VIP
  message("Extracting variable importance...")

  vip_tbl <- broom::tidy(final_fit$fit$fit) |>
    dplyr::filter(term != "(Intercept)") |>
    dplyr::mutate(vip = abs(estimate)) |>
    dplyr::arrange(desc(vip))

  if (model_mode == "classification") {
    preds <- predict(
      final_fit,
      data,
      type = "prob") |>
      dplyr::bind_cols(predict(final_fit, data)) |>
      dplyr::bind_cols(
        data |>
          dplyr::select(dplyr::all_of(outcome)))

    positive_class <- levels(data[[outcome]])[2]

    glm_model <- stats::glm(
      stats::as.formula(
        paste(outcome, "~", paste(predictors, collapse = "+"))),
      data = data,
      family = "binomial"
    )

    return(list(
      model = final_fit,
      best_penalty = best_penalty,
      predictions = preds,
      confusion = yardstick::conf_mat(
        preds,
        truth = !!outcome_sym,
        estimate = .pred_class),
      roc_auc = yardstick::roc_auc(
        preds,
        truth = !!outcome_sym,
        !!rlang::sym(paste0(".pred_", positive_class))),
      rsq = performance::r2_mcfadden(glm_model),
      vip = vip_tbl
    ))

  } else {
    preds <- predict(final_fit, data) |>
      dplyr::bind_cols(data |>
                         dplyr::select(dplyr::all_of(outcome)))

    return(list(
      model = final_fit,
      best_penalty = best_penalty,
      predictions = preds,
      rmse = yardstick::rmse(
        preds,
        truth = !!outcome_sym,
        estimate = .pred),
      rsq = yardstick::rsq(
        preds,
        truth = !!outcome_sym,
        estimate = .pred),
      vip = vip_tbl
    ))
  }
}
