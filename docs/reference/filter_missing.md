# Filter Features and Variables with High Missingness

Removes exposure variables and omics features with missing values above
a specified threshold. Generates missing data summaries and quality
control (QC) plots.

## Usage

``` r
filter_missing(exposomicset, na_thresh = 20, na_plot_thresh = 5)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing exposure and omics data.

- na_thresh:

  A numeric value specifying the percentage of missing data allowed
  before a variable or feature is removed. Default is `20`.

- na_plot_thresh:

  A numeric value specifying the minimum missing percentage for
  inclusion in QC plots. Default is `5`.

## Value

A `MultiAssayExperiment` object with filtered exposure variables and
omics features. QC results, including missingness summaries and plots,
are stored in `metadata(exposomicset)$na_qc`.

## Details

The function assesses missingness in both `colData(exposomicset)`
(exposure data) and `experiments(exposomicset)` (omics data).

- Exposure variables with more than `na_thresh`% missing values are
  removed.

- Omics features (rows in assay matrices) exceeding `na_thresh`% missing
  values are filtered.

- Missingness summaries and QC plots are generated using
  [`naniar::gg_miss_var()`](https://naniar.njtierney.com/reference/gg_miss_var.html)
  and stored in metadata.

## Examples

``` r
# Create example data
mae <- make_example_data(
    n_samples = 20,
    return_mae = TRUE
)
#> Ensuring all omics datasets are matrices with column names.
#> Creating SummarizedExperiment objects.
#> Creating MultiAssayExperiment object.
#> MultiAssayExperiment created successfully.

# Introduce some missingness
MultiAssayExperiment::colData(mae)$exposure_pm25[sample(1:20, 5)] <- NA

# Filter features and exposures with high missingness
mae_filtered <- filter_missing(
    exposomicset = mae,
    na_thresh = 20,
    na_plot_thresh = 5
)
#> Missing Data Filter threshold: 20%
#> Filtered metadata variables: exposure_pm25
#> Filtered rows with high missingness in mRNA: 0
#> Filtered rows with high missingness in proteomics: 0
```
