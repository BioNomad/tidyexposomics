# Compute Composite Exposome Scores

Calculates a summary exposome score per sample using one of several
methods including mean, sum, median, PCA (first principal component),
IRT (Item Response Theory), quantile binning, or row-wise variance. The
resulting score is added to the `colData` of the `MultiAssayExperiment`
object.

## Usage

``` r
run_exposome_score(
  exposomicset,
  score_type,
  exposure_cols = NULL,
  scale = TRUE,
  score_column_name = NULL
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing exposure data in its
  `colData`.

- score_type:

  Character. The method used to compute the score. Options are:
  `"mean"`, `"sum"`, `"median"`, `"pca"`, `"irt"`, `"quantile"`,
  `"var"`.

- exposure_cols:

  Optional character vector. Specific exposure column names to include.
  If `NULL`, all numeric columns are used.

- scale:

  Logical. Whether to scale the exposures before computing the score.
  Default is `TRUE`.

- score_column_name:

  Optional name for the resulting score column. If `NULL`, an automatic
  name is used (e.g., `"exposome_score_pca"`).

## Value

A `MultiAssayExperiment` object with the exposome score added to
`colData()`.

## Details

- `"pca"` uses the first principal component from
  [`prcomp()`](https://rdrr.io/r/stats/prcomp.html).

- `"irt"` uses the `mirt` package to fit a graded response model to
  discretized exposures.

- `"quantile"` assigns decile bins (1-10) to each variable and sums them
  row-wise.

- `"var"` computes the row-wise variance across exposures.

## Examples

``` r
# create the example data
mae <- make_example_data(
    n_samples = 10,
    return_mae = TRUE
)
#> Ensuring all omics datasets are matrices with column names.
#> Creating SummarizedExperiment objects.
#> Creating MultiAssayExperiment object.
#> MultiAssayExperiment created successfully.

# create the air pollution score
mae <- run_exposome_score(
    mae,
    score_type = "pca",
    exposure_cols = c("exposure_pm25", "exposure_no2"),
    scale = TRUE,
    score_column_name = "air_pollution_score"
)
#> Extracting exposure data...
#> Calculating PCA exposure scores...
```
