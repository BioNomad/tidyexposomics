# Calculate Exposure Impact from Feature-Exposure Correlation Networks

Generalized centrality-based exposure impact analysis using DEG, omics,
or factor features.

## Usage

``` r
run_exposure_impact(
  exposomicset,
  feature_type = c("degs", "omics", "factor_features"),
  pval_col = "adj.P.Val",
  pval_thresh = 0.1,
  action = c("add", "get")
)
```

## Arguments

- exposomicset:

  A `MultiAssayExperiment` object with correlation and network metadata.

- feature_type:

  One of `"degs"`, `"omics"`, or `"factor_features"`.

- pval_col:

  Column in differential abundance results to filter DEGs. Default =
  `"adj.P.Val"`.

- pval_thresh:

  DEG p-value threshold. Ignored unless `feature_type == "degs"`.

- action:

  Either `"add"` (store in metadata) or `"get"` (return list).

## Value

Either an updated MultiAssayExperiment (if action = "add") or a list.

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

# perform impact analysis
mae <- mae |>
    run_exposure_impact(
        feature_type = "omics"
    )
```
