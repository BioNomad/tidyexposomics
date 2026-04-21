# Identify and Annotate Shared Top Features Across Integration Factors

Identifies top features shared across factors based on integration
method. For MOFA/MCIA, takes intersection across factors. For
DIABLO/RGCCA, takes features recurring in more than 2 block-specific
components.

## Usage

``` r
run_factor_overlap(
  exposomicset,
  robust = TRUE,
  stability_score = NULL,
  score_col = "stability_score",
  pval_thresh = 0.05,
  logfc_thresh = log2(1.5),
  pval_col = "padj",
  logfc_col = "logFC",
  action = "add"
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` with integration results and top factor
  features.

- robust:

  Logical; if `TRUE`, uses sensitivity score. Otherwise, uses DEG
  thresholds.

- stability_score:

  Optional numeric threshold (overrides default from metadata).

- score_col:

  Column name for sensitivity score. Default is `"stability_score"`.

- pval_thresh:

  DEG p-value threshold (if `robust = FALSE`). Default is `0.05`.

- logfc_thresh:

  DEG logFC threshold (if `robust = FALSE`). Default is `log2(1.5)`.

- pval_col:

  Column name for p-value. Default is `"padj"`.

- logfc_col:

  Column name for logFC. Default is `"logFC"`.

- action:

  `"add"` to return modified object, `"get"` to return data.frame.

## Value

Modified `MultiAssayExperiment` or `data.frame` of shared top features.

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

# perform multiomics integration
mae <- run_multiomics_integration(
    mae,
    method = "DIABLO",
    outcome = "smoker",
    n_factors = 3
)
#> Scaling each assay in MultiAssayExperiment.
#> Running multi-omics integration using DIABLO...
#> Applying DIABLO supervised integration.
#> Design matrix has changed to include Y; each block will be
#>             linked to Y.


# identify the features that contribute most to the factors
mae <- extract_top_factor_features(
    mae,
    factors = c("V1", "V2", "V3"),
    method = "percentile",
    percentile = 0.5,
    action = "add"
)
#> Extracting top contributing features for specified factors.
#> Using DIABLO loadings.
#> Applying percentile-based filtering (>50%).
#> Selected 0 features contributing to specified factors.

# perform differential abundance analysis
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

# determine the overlap in features
mae <- mae |>
    run_factor_overlap(
        robust = FALSE,
        pval_col = "adj.P.Val"
    )
```
