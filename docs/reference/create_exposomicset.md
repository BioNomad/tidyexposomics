# Create an Exposomicset Object

Constructs a `MultiAssayExperiment` object from exposure data and
optionally omics datasets, ensuring proper formatting and alignment of
samples and features. For epidemiology-only workflows, omics data can be
omitted.

## Usage

``` r
create_exposomicset(codebook, exposure, omics = NULL, row_data = NULL)
```

## Arguments

- codebook:

  A data frame containing variable information metadata.

- exposure:

  A data frame containing exposure data, with rows as samples and
  columns as variables.

- omics:

  An optional list of matrices or a single matrix representing omics
  data. Each matrix should have samples as columns and features as rows.
  If `NULL`, creates an exposure-only exposomicset. Default is `NULL`.

- row_data:

  An optional list of `DataFrame` objects providing feature metadata for
  each omics dataset. If `NULL`, row metadata is generated
  automatically. Default is `NULL`.

## Value

A `MultiAssayExperiment` object containing the formatted exposure and
optionally omics datasets.

## Details

The function validates inputs and creates a `MultiAssayExperiment`
object. If omics data is provided, it converts matrices into
`SummarizedExperiment` objects with proper sample alignment. If omics is
`NULL`, the function creates an exposure-only object suitable for
epidemiological analyses using
[`run_association()`](https://BioNomad.github.io/tidyexposomics/reference/run_association.md)
with `source = "exposures"`.

## Examples

``` r
# Epi user workflow
# so no omics data
epi_data <- data.frame(
    pm25 = rnorm(10),
    outcome = rbinom(10, 1, 0.5),
    age = rnorm(10, 45, 10),
    row.names = paste0("subj_", 1:10)
)

codebook <- data.frame(
    variable = c("pm25", "outcome", "age"),
    category = c("exposure", "outcome", "covariate")
)

mae <- create_exposomicset(
    codebook = codebook,
    exposure = epi_data
)
#> No omics data provided. Creating exposure-only exposomicset.
#> MultiAssayExperiment created successfully.

# Multi-omics workflow
tmp <- make_example_data(n_samples = 10)

mae <- create_exposomicset(
    codebook = tmp$codebook,
    exposure = tmp$exposure,
    omics = tmp$omics,
    row_data = tmp$row_data
)
#> Ensuring all omics datasets are matrices with column names.
#> Creating SummarizedExperiment objects.
#> Creating MultiAssayExperiment object.
#> MultiAssayExperiment created successfully.
```
