# Run Association Analysis

Perform GLM-based association testing between a specified outcome and
features from exposures, omics, or latent factors. Automatically adjusts
for covariates and supports both Gaussian and binomial models.

## Usage

``` r
run_association(
  exposomicset,
  outcome,
  source = c("omics", "exposures", "factors"),
  covariates = NULL,
  feature_set = NULL,
  log_trans = TRUE,
  top_n = NULL,
  family = "gaussian",
  correction_method = "fdr",
  action = "add",
  feature_col = NULL,
  mirna_assays = NULL
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing data and metadata.

- outcome:

  The outcome variable name (must be in `colData`).

- source:

  Source of features to test. One of `"omics"`, `"exposures"`,
  `"factors"`.

- covariates:

  Optional vector of covariate names to include in the model.

- feature_set:

  Optional character vector of exposure or GO terms to test.

- log_trans:

  Optional boolean value dictating whether or not to log2 transform
  omics features.

- top_n:

  Optional integer: if using omics source, select top `n` most variable
  features.

- family:

  GLM family; `"gaussian"` or `"binomial"`.

- correction_method:

  Method for p-value adjustment (default: `"fdr"`).

- action:

  If `"add"` (default), saves results to metadata; else returns results
  as list.

- feature_col:

  The column in `rowData` for matching gene symbols or IDs.

- mirna_assays:

  Optional character vector of assays to exclude when extracting GO
  terms.

## Value

If `action = "add"`, returns updated `MultiAssayExperiment`. Otherwise,
returns a list of:

- `results_df`: tidy summary of associations

- `covariates`: the covariates used

- `model_data`: model matrix used in the GLMs

## Examples

``` r
#' # create example data
mae <- make_example_data(
    n_samples = 10,
    return_mae = TRUE
)
#> Ensuring all omics datasets are matrices with column names.
#> Creating SummarizedExperiment objects.
#> Creating MultiAssayExperiment object.
#> MultiAssayExperiment created successfully.

# run association tests
mae <- mae |>
    run_association(
        source = "exposures",
        feature_set = c("exposure_pm25", "exposure_no2"),
        outcome = "smoker",
        covariates = c("age"),
        family = "binomial"
    )
#> Running GLMs.
```
