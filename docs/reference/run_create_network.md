# Create Correlation Network from Feature Data

Constructs an undirected feature-feature or feature-exposure correlation
network from correlation results stored in a `MultiAssayExperiment`
object. The function supports multiple correlation formats depending on
`feature_type`, and stores or returns an `igraph` object with associated
node and edge metadata.

## Usage

``` r
run_create_network(
  exposomicset,
  feature_type = c("degs", "omics", "factors", "factor_features", "exposures",
    "degs_feature_cor", "omics_feature_cor", "factor_features_feature_cor"),
  action = c("add", "get")
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object containing correlation results in
  metadata.

- feature_type:

  Type of correlation result to convert to a network. One of: `"degs"`,
  `"omics"`, `"factors"`, `"factor_features"`, `"exposures"`,
  `"degs_feature_cor"`, `"omics_feature_cor"`, or
  `"factor_features_feature_cor"`.

- action:

  Whether to `"add"` the network to the object or `"get"` it as a list.

## Value

If `action = "add"`, returns the updated `MultiAssayExperiment` with a
new `network` entry in metadata. If `action = "get"`, returns a list
with `graph` (an `igraph` object) and `summary` (a tibble).

## Details

The function detects the appropriate edge and node structure based on
column names in the correlation results. Edge weights are based on
correlation coefficients and include FDR values.

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

# perform correlation analyses
# correlate with exposures
mae <- mae |>
    run_correlation(
        feature_type = "omics",
        variable_map = mae |>
            pivot_feature() |>
            dplyr::select(
                variable = .feature,
                exp_name = .exp_name
            ),
        exposure_cols = c("exposure_pm25", "exposure_no2", "age", "bmi")
    ) |>
    # correlate omics features with themselves
    run_correlation(
        feature_type = "omics",
        variable_map = mae |>
            pivot_feature() |>
            dplyr::select(
                variable = .feature,
                exp_name = .exp_name
            ),
        feature_cors = TRUE,
        exposure_cols = c("exposure_pm25", "exposure_no2", "age", "bmi")
    )
#> Log2-Transforming each assay in MultiAssayExperiment.
#> Log2-Transforming each assay in MultiAssayExperiment.

# create the networks
mae <- mae |>
    run_create_network(
        feature_type = "omics_feature_cor",
        action = "add"
    ) |>
    run_create_network(
        feature_type = "omics",
        action = "add"
    )
#> Creating network from correlation results.
#> Network added to metadata as: network_omics_feature_cor
#> Creating network from correlation results.
#> Network added to metadata as: network_omics
```
