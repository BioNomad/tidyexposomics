# Run Correlation Analysis

Computes correlations between exposures and feature types including
DEGs, omics, latent factors, top factor features, or principal
components (PCs). Optionally computes feature-feature correlations to
support network analysis.

## Usage

``` r
run_correlation(
  exposomicset,
  feature_type = c("degs", "omics", "factors", "factor_features", "exposures", "pcs"),
  exposure_cols = NULL,
  variable_map = NULL,
  n_pcs = NULL,
  feature_cors = FALSE,
  robust = FALSE,
  score_col = "stability_score",
  score_thresh = NULL,
  correlation_method = "spearman",
  correlation_cutoff = 0.3,
  cor_pval_column = "p.value",
  pval_cutoff = 0.05,
  deg_pval_col = "adj.P.Val",
  deg_logfc_col = "logFC",
  deg_pval_thresh = 0.05,
  deg_logfc_thresh = log2(1.5),
  batch_size = 1500,
  action = c("add", "get")
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object.

- feature_type:

  Type of features to correlate. One of `"degs"`, `"omics"`,
  `"factors"`, `"factor_features"`, `"exposures"`, or `"pcs"`.

- exposure_cols:

  Optional vector of exposure column names (from `colData`) to use.

- variable_map:

  Optional mapping of features to include by assay for `omics` mode.

- n_pcs:

  Number of PCs to use when `feature_type = "pcs"`.

- feature_cors:

  Logical; if `TRUE`, compute correlations between features rather than
  with exposures.

- robust:

  Logical; restrict DEGs to those passing sensitivity threshold.

- score_col:

  Column name in sensitivity analysis with feature stability score.

- score_thresh:

  Threshold for filtering robust features.

- correlation_method:

  One of `"pearson"`, `"spearman"`, or `"kendall"`.

- correlation_cutoff:

  Minimum absolute correlation to retain.

- cor_pval_column:

  Column in output to filter by p-value (default: `"p.value"`).

- pval_cutoff:

  Maximum p-value or FDR threshold to retain a correlation.

- deg_pval_col:

  Column with DEG adjusted p-values.

- deg_logfc_col:

  Column with DEG log fold-changes.

- deg_pval_thresh:

  P-value cutoff for DEGs.

- deg_logfc_thresh:

  Log fold-change cutoff for DEGs.

- batch_size:

  Number of features to process per batch (default: 1500).

- action:

  Whether to `"add"` results to metadata or `"get"` as a data frame.

## Value

If `action = "add"`, returns updated `MultiAssayExperiment` with results
added to metadata. If `action = "get"`, returns a tidy `data.frame` of
correlations.

## Examples

``` r
# create example data
mae <- make_example_data(
    n_samples = 10,
    return_mae = TRUE
)
#> Ensuring all omics datasets are matrices with column names.
#> Creating SummarizedExperiment objects.
#> Creating MultiAssayExperiment object.
#> MultiAssayExperiment created successfully.

# run correlation analysis
mae <- mae |>
    run_correlation(
        feature_type = "exposures",
        exposure_cols = c("exposure_pm25", "exposure_no2", "age", "bmi")
    )
```
