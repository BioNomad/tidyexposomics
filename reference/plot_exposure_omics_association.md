# Plot Exposure-Omics Associations

Plots the number of significant exposure-omics associations, grouped by
exposure, exposure category, or omics experiment.

## Usage

``` r
plot_exposure_omics_association(
  exposomicset,
  plot_type = c("exposures", "category", "exp_name"),
  pval_col = "p_adjust",
  pval_thresh = 0.05,
  estimate_thresh = 0
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing association results.

- plot_type:

  Character. One of `"exposures"`, `"category"`, or `"exp_name"`.
  Controls whether associations are summarized per exposure, per
  exposure category, or per omics experiment (assay). Defaults to
  `"exposures"`.

- pval_col:

  Character. Name of the column used for p-value filtering. Defaults to
  `"p_adjust"`.

- pval_thresh:

  Numeric. Significance threshold applied to `pval_col`. Rows with
  values below this threshold are retained. Defaults to `0.05`.

- estimate_thresh:

  Numeric. Minimum absolute effect size. Rows with `abs(estimate)`
  greater than or equal to this value are retained. Defaults to `0` (no
  filtering on estimate).

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
        formula = ~ .feature + age + sex,
        exposures = c("exposure_pm25", "exposure_no2")
    )
#> Testing 2 exposures across 2 assays using model: ~.feature + age + sex
#> Processing assay: mRNA
#> Processing assay: proteomics

plot_exposure_omics_association(
    exposomicset = mae,
    plot_type = "exp_name",
    estimate_thresh = 0.1
)

```
