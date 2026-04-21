# Run Exposure-Omics Association

Test associations between each exposure and each omics feature using
limma's linear modeling framework.

## Usage

``` r
run_exposure_omics_association(
  exposomicset,
  exposures = NULL,
  exp_name = NULL,
  covariates = NULL,
  scaling_method = "none",
  correction_method = "fdr",
  top_pct = NULL,
  filter_by = c("variance", "mean"),
  action = "add"
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing exposome and omics data.

- exposures:

  Character vector of exposure variable names to test. If `NULL`, uses
  all variables from the codebook.

- exp_name:

  Name(s) of the omics assay(s) to test against. If `NULL`, uses all
  assays.

- covariates:

  Optional character vector of covariate names to include in the model.

- scaling_method:

  Character. Scaling method to apply before modeling. Options include
  `"none"` (default).

- correction_method:

  Method for p-value adjustment. Default is `"fdr"`.

- top_pct:

  Top X% of features to retain using either mean or variance which is
  specified by `filter_by`. If `NULL`, no features will be filtered.

- filter_by:

  Determination of how to filter omics features either by mean or
  variance.

- action:

  If `"add"` (default), saves results to metadata, if `"get"`, returns
  results as a data frame.

## Value

If `action = "add"`, returns updated `MultiAssayExperiment`. Otherwise,
returns a tibble with association results.

## Details

This function uses limma to test associations between multiple exposures
and omics features. For each exposure, a linear model is fit with the
exposure as the predictor and each omics feature as the outcome,
adjusting for covariates.

`omics_feature ~ exposure + covariate1 + covariate2 + ...`

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

# run exposure-omics association
mae <- mae |>
    run_exposure_omics_association(
        exposures = c("exposure_pm25", "exposure_no2"),
        covariates = c("age", "sex")
    )
#> Testing 2 exposures across 2 assays
#> Processing assay: mRNA
#> Processing assay: proteomics
```
