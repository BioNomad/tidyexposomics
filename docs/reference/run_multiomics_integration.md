# Run Multi-Omics Integration

Performs multi-omics integration using one of several available methods:
MOFA, MCIA, RGCCA, or DIABLO. This function takes a
`MultiAssayExperiment` object with two or more assays and computes
shared latent factors across omics layers.

## Usage

``` r
run_multiomics_integration(
  exposomicset,
  method = "MCIA",
  n_factors = 10,
  scale = TRUE,
  outcome = NULL,
  max.iter = 500,
  near.zero.var = TRUE,
  action = "add"
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object with at least two assays.

- method:

  Character. Integration method to use. Options are `"MOFA"`, `"MCIA"`,
  `"RGCCA"`, or `"DIABLO"`.

- n_factors:

  Integer. Number of latent factors/components to compute. Default is
  10.

- scale:

  Logical. Whether to scale each assay before integration. Default is
  `TRUE`.

- outcome:

  Character. Required if `method = "DIABLO"`. Name of outcome variable
  in `colData` used for supervised integration.

- max.iter:

  numeric. Option to increase the number of iterations for
  [`mixOmics::block.splsda`](https://rdrr.io/pkg/mixOmics/man/block.splsda.html)
  the default is 500.

- near.zero.var:

  Logical. Option to remove variables with near zero variance for
  [`mixOmics::block.splsda`](https://rdrr.io/pkg/mixOmics/man/block.splsda.html),
  default is `TRUE` .

- action:

  Character. Whether to `"add"` results to the metadata or `"get"` them
  as a list. Default is `"add"`.

## Value

If `action = "add"`, returns a `MultiAssayExperiment` with integration
results stored in
`metadata(exposomicset)$multiomics_integration$integration_results`. If
`action = "get"`, returns a list with integration `method` and `result`.

## Details

- `"MOFA"` runs Multi-Omics Factor Analysis using the `MOFA2` package
  and returns a trained model.

- `"MCIA"` runs multi-co-inertia analysis using the `nipalsMCIA`
  package.

- `"RGCCA"` runs Regularized Generalized Canonical Correlation Analysis
  using the `RGCCA` package.

- `"DIABLO"` performs supervised integration using the `mixOmics`
  package and a specified outcome.

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
```
