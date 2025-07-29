test_that("run_sensitivity_analysis returns results with action = 'get'", {
  skip_if_not_installed("limma")

  # create dummy data
  dummy <- make_dummy_data(n_samples = 20)
  mae <- create_expomicset(
    codebook = dummy$codebook,
    exposure = dummy$exposure,
    omics = dummy$omics,
    row_data = dummy$row_data
  )

  # Run differential abundance
  mae <- run_differential_abundance(
    expomicset = mae,
    formula = ~ smoker + sex,
    abundance_col = "counts",
    method = "limma_voom",
    action = "add"
  )

  # Run the sensitivity analysis
  mae <- run_sensitivity_analysis(
    expomicset = mae,
    base_formula = ~ smoker + sex,
    methods = c("limma_voom"),
    scaling_methods = c("none"),
    covariates_to_remove = "sex",
    pval_col = "P.Value",
    logfc_col = "logFC",
    pval_threshold = 0.05,
    logFC_threshold = 0,
    bootstrap_n = 3,
    action = "add"
  )

  result <- mae@metadata$differential_analysis$sensitivity_analysis

  # check that the results are a list
  expect_type(result, "list")

  # check that each item is the expected data type
  expect_s3_class(result$sensitivity_df, "data.frame")
  expect_s3_class(result$feature_stability, "data.frame")
  expect_type(result$score_thresh, "double")

  # grab rowdata column names
  rowdata_cols <- mae |>
    pivot_feature() |>
    dplyr::select(-c(.exp_name,.feature)) |>
    colnames()

  # check that the expected rowdata column names are present
  expect_contains(colnames(result$sensitivity_df),rowdata_cols)

  # check that the condition column names are present
  expect_contains(colnames(result$sensitivity_df),
                  c("feature","contrast","method",
                    "scaling" ,"model","exp_name",
                    "bootstrap_id"))

  # ensure that the scores were calculated
  scores <- c('presence_rate','effect_consistency','stability_score','mean_log_p','logp_weighted_score','sd_logFC','iqr_logFC','cv_logFC','sign_flip_freq','sd_log_p')
  lapply(scores,function(score){
    expect_type(result$feature_stability[[score]],"double")
  })

})
