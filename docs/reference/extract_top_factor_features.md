# Extract Top Contributing Features for Factors

Identifies the most influential features for specified factors using
multiomics integration results. Features are selected based on either a
percentile cutoff or an absolute loading threshold.

## Usage

``` r
extract_top_factor_features(
  exposomicset,
  factors = NULL,
  pval_col = "p_adjust",
  pval_thresh = 0.05,
  method = "percentile",
  percentile = 0.9,
  threshold = 0.3,
  action = "add"
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing integration results.

- factors:

  A character vector specifying the factors of interest. If `NULL`,
  factors are automatically selected from the association results using
  the `pval_col` and `pval_thresh` criteria.

- pval_col:

  A string specifying the column name of the p-value or adjusted p-value
  used for factor selection if `factors` is `NULL`. Default is
  `"p_adjust"`.

- pval_thresh:

  A numeric value specifying the significance threshold for selecting
  factors from association results when `factors` is `NULL`. Default is
  `0.05`.

- method:

  A character string specifying the feature selection method
  (`"percentile"` or `"threshold"`). Default is `"percentile"`.

- percentile:

  A numeric value between 0 and 1 indicating the percentile threshold
  for feature selection when `method = "percentile"`. Default is `0.9`.

- threshold:

  A numeric value specifying the absolute loading cutoff for feature
  selection when `method = "threshold"`. Default is `0.3`.

- action:

  A character string indicating whether to return results (`"get"`) or
  add them to metadata (`"add"`). Default is `"add"`.

## Value

If `action = "add"`, returns the modified `exposomicset` with selected
features stored in metadata. If `action = "get"`, returns a data frame
containing:

- feature:

  The selected feature contributing to the factor.

- factor:

  The factor to which the feature contributes.

- loading:

  The factor loading value of the feature.

- exp_name:

  The experiment from which the feature originated.

## Details

The function extracts factor loadings from `metadata(exposomicset)`,
applies filtering based on the selected method, and identifies top
contributing features for each specified factor.

If `factors` is not provided, the function will automatically select
statistically significant factors from
`metadata(exposomicset)$association$assoc_factors$results_df` using the
specified `pval_col` and `pval_thresh` as criteria.

Features can be selected using:

- **Percentile-based filtering** (`method = "percentile"`): Selects
  features with absolute loadings above a specified percentile.

- **Threshold-based filtering** (`method = "threshold"`): Selects
  features with absolute loadings exceeding a fixed value.

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

top_feats <- extract_top_factor_features(
    mae,
    factors = c("V1", "V2", "V3"),
    method = "percentile",
    percentile = 0.9,
    action = "get"
)
#> Extracting top contributing features for specified factors.
#> Using DIABLO loadings.
#> Applying percentile-based filtering (>90%).
#> Selected 0 features contributing to specified factors.
```
