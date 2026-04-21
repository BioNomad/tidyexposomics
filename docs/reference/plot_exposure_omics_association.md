# Plot Exposure-Omics Associations

Plots the number of significant exposure-omics associations, grouped
either by exposure or the exposure category.

## Usage

``` r
plot_exposure_omics_association(
  exposomicset,
  plot_type = c("exposures", "category"),
  pval_col = "p_adjust",
  pval_thresh = 0.05
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing association results.

- plot_type:

  Character. One of `"exposures"` or `"category"`. Controls whether
  associations are summarized per exposure or per exposure category.
  Defaults to `"exposures"`.

- pval_col:

  Character. Name of the column used for p-value filtering. Defaults to
  `"p_adjust"`.

- pval_thresh:

  Numeric. Significance threshold applied to `pval_col`. Rows with
  values below this threshold are retained. Defaults to `0.05`.

## Value

A `ggplot` object.

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

# run exposure-omics association
mae <- mae |>
    run_exposure_omics_association(
        exposures = c("exposure_pm25", "exposure_no2"),
        covariates = c("age", "sex")
    )
#> Testing 2 exposures across 2 assays
#> Processing assay: mRNA
#> Processing assay: proteomics

plot_exposure_omics_association(
    exposomicset = mae,
    plot_type = "exposures"
)

```
