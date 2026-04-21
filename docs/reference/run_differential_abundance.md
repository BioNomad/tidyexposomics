# Run Differential Abundance Analysis

Performs differential abundance testing across all assays in a
`MultiAssayExperiment` object using a specified statistical method. The
function updates each assay with its corresponding `colData`, fits the
model using the provided formula, and combines the results into a
unified table.

## Usage

``` r
run_differential_abundance(
  exposomicset,
  formula,
  abundance_col = "counts",
  method = "limma_trend",
  contrasts = NULL,
  scaling_method = "none",
  action = "add"
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` containing assays to test.

- formula:

  A model formula for the differential analysis (e.g., ~ group + batch).

- abundance_col:

  Character. The name of the assay matrix to use for abundance values.
  Default is `"counts"`.

- method:

  Character. Differential analysis method to use. Currently supports
  `"limma_trend"` (default).

- contrasts:

  A named list of contrasts for pairwise comparisons. Default is `NULL`
  (uses default group comparisons).

- scaling_method:

  Character. Scaling method to apply before modeling. Options include
  `"none"` (default), `"zscore"`, etc.

- action:

  Character. Whether to `"add"` results to `exposomicset` metadata or
  `"get"` the results as a data frame. Default is `"add"`.

## Value

Either the updated `MultiAssayExperiment` (if `action = "add"`) or a
tibble with differential abundance results (if `action = "get"`).

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

# perform differential abundance analysis
mae <- run_differential_abundance(
    exposomicset = mae,
    formula = ~ smoker + sex,
    abundance_col = "counts",
    method = "limma_trend",
    action = "add"
)
#> Running differential abundance testing.
#> Processing assay: mRNA
#> Processing assay: proteomics
#> Differential abundance testing completed.
```
