# Pivot a selected omics dataset from a MultiAssayExperiment into tidybulk format

Extracts a specified omics dataset from a `MultiAssayExperiment`,
optionally filters by feature (row) names, and returns a tidy tibble.
The output includes assay values along with sample metadata and feature
metadata.

## Usage

``` r
pivot_exp(exposomicset, exp_name, features = NULL)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing one or more omics assays.

- exp_name:

  A character string. The name of the omics dataset to extract (e.g.,
  "Proteomics").

- features:

  Optional character vector of row (feature) names to retain. If `NULL`,
  all features are included.

## Value

A tibble in tidy format with one row per feature/sample pair, including
all metadata and a new column `exp_name` indicating the assay source.
Assay values are provided in separate columns named after the assay
slot(s).

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

# pivot experiment
exp_data <- pivot_exp(
    exposomicset = mae,
    exp_name = "mRNA",
    features = c("feat_42")
)
```
