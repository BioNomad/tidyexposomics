# Plot Sensitivity Analysis Summary

Generates a ridge plot and bar chart summarizing feature stability
scores across assays.

## Usage

``` r
plot_sensitivity_summary(
  exposomicset,
  stability_score_thresh = NULL,
  stability_metric = "stability_score",
  title = "Distribution of Stability Scores"
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing sensitivity analysis
  results in `metadata(exposomicset)$sensitivity_analysis`.

- stability_score_thresh:

  A numeric threshold for stability scores. Default is `NULL`, which
  uses the threshold stored in
  `metadata(exposomicset)$sensitivity_analysis$score_thresh`.

- stability_metric:

  A character string specifying which stability metric to plot (e.g.,
  "stability_score", "logp_weighted_score"). Default is
  "stability_score".

- title:

  A character string specifying the title of the ridge plot. Default is
  "Distribution of Stability Scores".

## Value

A `patchwork` object combining a ridge plot and a bar chart.

## Details

This function:

- Extracts feature stability scores from
  `metadata(exposomicset)$sensitivity_analysis$feature_stability`.

- Displays a **ridge plot** of stability score distributions per assay.

- Displays a **bar chart** of the number of features per assay.

- Prints the number of features with stability scores above the
  threshold.

## Examples

``` r
# create example data
mae <- make_example_data(
    n_samples = 20,
    return_mae = TRUE
)
#> Ensuring all omics datasets are matrices with column names.
#> Creating SummarizedExperiment objects.
#> Creating MultiAssayExperiment object.
#> MultiAssayExperiment created successfully.


# Run differential abundance
mae <- run_differential_abundance(
    exposomicset = mae,
    formula = ~ smoker + sex,
    abundance_col = "counts",
    method = "limma_voom",
    action = "add"
)
#> Running differential abundance testing.
#> Processing assay: mRNA
#> Warning: All samples appear to belong to the same group.
#> tidybulk says: The design column names are "(Intercept), smokeryes, sexM"
#> tidybulk says: to access the DE object do `metadata(.)$tidybulk$limma_voom_object`
#> tidybulk says: to access the raw results (fitted GLM) do `metadata(.)$tidybulk$limma_voom_fit`
#> Processing assay: proteomics
#> Warning: All samples appear to belong to the same group.
#> tidybulk says: The design column names are "(Intercept), smokeryes, sexM"
#> tidybulk says: to access the DE object do `metadata(.)$tidybulk$limma_voom_object`
#> tidybulk says: to access the raw results (fitted GLM) do `metadata(.)$tidybulk$limma_voom_fit`
#> Differential abundance testing completed.

# Run the sensitivity analysis
mae <- run_sensitivity_analysis(
    exposomicset = mae,
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
#> Running bootstrap iteration 1 of 3
#> Warning: All samples appear to belong to the same group.
#> tidybulk says: The design column names are "(Intercept), smokeryes, sexM"
#> tidybulk says: to access the DE object do `metadata(.)$tidybulk$limma_voom_object`
#> tidybulk says: to access the raw results (fitted GLM) do `metadata(.)$tidybulk$limma_voom_fit`
#> Warning: All samples appear to belong to the same group.
#> tidybulk says: The design column names are "(Intercept), smokeryes, sexM"
#> tidybulk says: to access the DE object do `metadata(.)$tidybulk$limma_voom_object`
#> tidybulk says: to access the raw results (fitted GLM) do `metadata(.)$tidybulk$limma_voom_fit`
#> Running bootstrap iteration 2 of 3
#> Warning: All samples appear to belong to the same group.
#> tidybulk says: The design column names are "(Intercept), smokeryes, sexM"
#> tidybulk says: to access the DE object do `metadata(.)$tidybulk$limma_voom_object`
#> tidybulk says: to access the raw results (fitted GLM) do `metadata(.)$tidybulk$limma_voom_fit`
#> Warning: All samples appear to belong to the same group.
#> tidybulk says: The design column names are "(Intercept), smokeryes, sexM"
#> tidybulk says: to access the DE object do `metadata(.)$tidybulk$limma_voom_object`
#> tidybulk says: to access the raw results (fitted GLM) do `metadata(.)$tidybulk$limma_voom_fit`
#> Running bootstrap iteration 3 of 3
#> Warning: All samples appear to belong to the same group.
#> tidybulk says: The design column names are "(Intercept), smokeryes, sexM"
#> tidybulk says: to access the DE object do `metadata(.)$tidybulk$limma_voom_object`
#> tidybulk says: to access the raw results (fitted GLM) do `metadata(.)$tidybulk$limma_voom_fit`
#> Warning: All samples appear to belong to the same group.
#> tidybulk says: The design column names are "(Intercept), smokeryes, sexM"
#> tidybulk says: to access the DE object do `metadata(.)$tidybulk$limma_voom_object`
#> tidybulk says: to access the raw results (fitted GLM) do `metadata(.)$tidybulk$limma_voom_fit`
#> Computing feature stability across sensitivity conditions.
#> Feature stability analysis completed.
#> Number Of Features Above Threshold Of 0.59:
#> ----------------------------------------
#> mRNA: 0/128
#> proteomics: 3/80

# create the sensitivity summary plot
sens_sum_p <- mae |>
    plot_sensitivity_summary()
#> Number of Features with stability_score > 0.588452248517675:
#> mRNA: 0 / 128
#> proteomics: 0 / 80
```
