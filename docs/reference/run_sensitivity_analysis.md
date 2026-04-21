# Run Sensitivity Analysis for Differential Abundance

Performs sensitivity analysis by systematically varying statistical
methods, scaling strategies, and model formulas (with optional bootstrap
sampling) to assess the stability of differentially abundant features.

## Usage

``` r
run_sensitivity_analysis(
  exposomicset,
  base_formula,
  abundance_col = "counts",
  methods = c("limma_trend", "limma_voom", "DESeq2", "edgeR_quasi_likelihood"),
  scaling_methods = c("none", "TMM", "quantile"),
  contrasts = NULL,
  covariates_to_remove = NULL,
  pval_col = "adj.P.Val",
  logfc_col = "logFC",
  pval_threshold = 0.05,
  logFC_threshold = log2(1),
  score_thresh = NULL,
  score_quantile = 0.9,
  stability_metric = "stability_score",
  action = "add",
  bootstrap_n = 1
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` containing the assays to analyze.

- base_formula:

  The base model formula used for differential analysis.

- abundance_col:

  Character. Name of the column in the assays representing abundance.
  Default is `"counts"`.

- methods:

  Character vector of differential expression methods. Options include
  `"limma_trend"` ,`"limma_voom"`, `"DESeq2"`, and
  `"edgeR_quasi_likelihood"`.

- scaling_methods:

  Character vector of normalization methods to try. Options include
  `"none"`, `"TMM"`, and `"quantile"`.

- contrasts:

  Optional list of contrasts to apply for differential testing.

- covariates_to_remove:

  Optional character vector of covariates to remove from the base
  formula to generate model variants.

- pval_col:

  Name of the column containing p-values or adjusted p-values used to
  define significance.

- logfc_col:

  Name of the column containing log fold changes.

- pval_threshold:

  Numeric threshold for significance. Default is 0.05.

- logFC_threshold:

  Numeric threshold for absolute log fold change. Default is `log2(1)`
  (i.e., 0).

- score_thresh:

  Optional threshold for the selected stability metric. If not provided,
  calculated using `score_quantile`.

- score_quantile:

  Quantile used to define the threshold if `score_thresh` is not
  provided. Default is 0.9.

- stability_metric:

  Character. Name of the column in `feature_stability` to use as the
  scoring metric. Default is `"stability_score"`.

- action:

  Whether to `"add"` results to `metadata()` or `"get"` them as a list.
  Default is `"add"`.

- bootstrap_n:

  Integer. Number of bootstrap iterations. If 0, no resampling is
  performed. Default is 1.

## Value

If `action = "add"`, returns a `MultiAssayExperiment` with results
stored in
`metadata(exposomicset)$differential_analysis$sensitivity_analysis`. If
`action = "get"`, returns a list with three elements:

- `sensitivity_df`:

  Data frame of all differential results across model/method
  combinations.

- `feature_stability`:

  Data frame summarizing feature stability scores.

- `score_thresh`:

  The threshold used to define stable features.

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
#> Number Of Features Above Threshold Of 0.4:
#> ----------------------------------------
#> mRNA: 2/128
#> proteomics: 1/80
```
