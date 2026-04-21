# Cluster Samples Based on Exposure Data

Performs hierarchical clustering of samples using exposure data from
`colData(exposomicset)`.

## Usage

``` r
run_cluster_samples(
  exposomicset,
  exposure_cols = NULL,
  dist_method = NULL,
  user_k = NULL,
  cluster_method = "ward.D",
  clustering_approach = "diana",
  action = "add"
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing omics and exposure data.

- exposure_cols:

  A character vector of column names in `colData(exposomicset)` to use
  for clustering.

- dist_method:

  A character string specifying the distance metric (`"euclidean"`,
  `"gower"`, etc.). If `NULL`, it is automatically determined.

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

- action:

  A character string specifying `"add"` (store results in metadata) or
  `"get"` (return clustering results). Default is `"add"`.

## Value

If `action="add"`, returns the updated `exposomicset`. If
`action="get"`, returns a list with:

- sample_cluster:

  A hierarchical clustering object (`hclust`).

- sample_groups:

  A named vector of sample cluster assignments.

## Details

This function:

- Extracts **numeric exposure data** from `colData(exposomicset)`.

- Computes a **distance matrix** (`"gower"` for mixed data,
  `"euclidean"` for numeric).

- Determines the **optimal number of clusters (`k`)** using the
  specified method.

- Performs **hierarchical clustering** (`hclust`) and assigns samples to
  clusters.

- Generates a **heatmap** of scaled exposure values.

- Stores results in `metadata(exposomicset)$sample_clustering` when
  `action="add"`.

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

# determine sample clusters
mae <- run_cluster_samples(
    exposomicset = mae,
    exposure_cols = c("exposure_pm25", "exposure_no2", "age", "bmi"),
    clustering_approach = "diana"
)
#> Starting clustering analysis...
#> Optimal number of clusters for samples: 6
```
