# Cluster Exposures Based on Scaled Exposure Data

Performs hierarchical clustering of exposure variables using scaled
exposure data from `colData(exposomicset)` and records each exposure's
cluster in the codebook (`metadata(exposomicset)$codebook`).

## Usage

``` r
run_cluster_exposures(
  exposomicset,
  exposure_cols = NULL,
  dist_method = "correlation",
  user_k = NULL,
  cluster_method = "ward.D",
  clustering_approach = "diana",
  codebook_var_col = "variable",
  action = "add"
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing omics and exposure data.

- exposure_cols:

  A character vector of column names in `colData(exposomicset)` to
  cluster.

- dist_method:

  A character string specifying the distance metric. `"correlation"`
  (default) uses `1 - abs(cor)`, so exposures that are strongly
  associated (positively or negatively) cluster together. Any method
  accepted by [`stats::dist()`](https://rdrr.io/r/stats/dist.html) can
  also be used.

- user_k:

  An integer specifying the number of clusters. If `NULL`, an optimal
  `k` is determined.

- cluster_method:

  A character string specifying the hierarchical clustering method.
  Default is `"ward.D"`.

- clustering_approach:

  A character string specifying the method for determining `k`
  (`"diana"`, `"gap"`, `"elbow"`, `"dynamic"`, or `"density"`). Default
  is `"diana"`.

- codebook_var_col:

  A character string naming the codebook column that holds exposure
  names. Default is `"variable"`.

- action:

  A character string specifying `"add"` (store results in metadata) or
  `"get"` (return clustering results). Default is `"add"`.

## Value

If `action="add"`, returns the updated `exposomicset` with an
`exposure_cluster` column in `metadata(exposomicset)$codebook`. If
`action="get"`, returns a list with:

- exposure_cluster:

  A hierarchical clustering object (`hclust`).

- exposure_groups:

  A named vector of exposure cluster assignments.

## Examples

``` r
mae <- make_example_data(
    n_samples = 10,
    return_mae = TRUE
)
#> Ensuring all omics datasets are matrices with column names.
#> Creating SummarizedExperiment objects.
#> Creating MultiAssayExperiment object.
#> MultiAssayExperiment created successfully.

mae <- run_cluster_exposures(
    exposomicset = mae,
    exposure_cols = c("exposure_pm25", "exposure_no2", "age", "bmi"),
    clustering_approach = "diana"
)
#> Starting exposure clustering analysis...
#> Optimal number of clusters for exposures: 2
```
