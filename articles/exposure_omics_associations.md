# Exposure-Omics Associations

[Previous: Omics
Preprocessing](https://BioNomad.github.io/tidyexposomics/articles/omics_preprocessing.md)

**Objectives**

- Test for differentially abundant omics features associated with an
  outcome using.
- Visualize differential abundance results across omics layers.
- Test associations between every exposure and every omics feature.

### Prerequisites

We covered [omics/exposure
preprocessing](https://BioNomad.github.io/tidyexposomics/articles/omics_preprocessing.md)
in a previous vignette. Here is the quick start code to catch you up to
speed!

Click to expand code

``` r
# Load Libraries
library(tidyverse)
library(tidyexposomics)

# Load example data
data("tidyexposomics_example")

# Create exposomic set object
expom <- create_exposomicset(
    codebook = tidyexposomics_example$annotated_cb,
    exposure = tidyexposomics_example$meta,
    omics = list(
        "Gene Expression" = tidyexposomics_example$exp_filt,
        "Methylation" = tidyexposomics_example$methyl_filt
    ),
    row_data = list(
        "Gene Expression" = tidyexposomics_example$exp_fdata,
        "Methylation" = tidyexposomics_example$methyl_fdata
    )
)

# Grab exposure variables
exp_vars <- tidyexposomics_example$annotated_cb |>
    filter(category %in% c(
        "aerosol",
        "main group molecular entity",
        "polyatomic entity"
    )) |>
    pull(variable) |>
    as.character()

# Impute missing values
expom <- run_impute_missing(
    exposomicset = expom,
    exposure_impute_method = "missforest",
    exposure_cols = exp_vars
)

# Check variable normality
expom <- run_normality_check(
    exposomicset = expom,
    action = "add"
)

# Transform variables
expom <- transform_exposure(
    exposomicset = expom,
    transform_method = "boxcox_best",
    exposure_cols = exp_vars
)
```

### Differential Abundance

We provide functionality to test for differentially abundant features
associated with an outcome across multiple omics layers. This is done
using the `run_differential_abundance` function, which fits a model
defined by the user (using the `formula` argument) and supports several
methods. Here we apply the `limma_trend` method, a widely used approach
for analyzing omics data. Users can also specify how features are scaled
(e.g. none, quantile, TMM) before fitting.

``` r
# Run differential abundance analysis
expom <- run_differential_abundance(
    exposomicset = expom,
    formula = ~hs_asthma,
    method = "limma_trend",
    scaling_method = "none",
    action = "add"
)
```

We can summarize the results of the differential abundance analysis with
a volcano plot, which highlights features with a high log fold change
and that are statistically significant. The `plot_volcano` function
generates this visualization, with options to set thresholds for
p-values and log fold changes, and to label a subset of top-ranked
features. In this example, we use the `feature_clean` column to display
interpretable feature names.

**Note:** we set the `pval_col` to `P.Value` for the purposes of this
example, but we recommend keeping the default of `adj.P.Val` to use the
adjusted p-values.

``` r
# Plot the volcano plot
plot_volcano(
    exposomicset = expom,
    top_n_label = 3,
    feature_col = "feature_clean",
    logFC_thresh = log2(1),
    pval_thresh = 0.05,
    pval_col = "P.Value",
    logFC_col = "logFC",
    nrow = 1
)
```

![Volcano plot of differentially abundant features across omics
layers.](exposure_omics_associations_files/figure-html/volcano-plot-1.png)

Volcano plot of differentially abundant features across omics layers.

## Exposure-Omics Association

------------------------------------------------------------------------

### Exposure-Omics Association

Above we saw that there are not too many omics features associated with
asthma. Which may be due to the subsampling in this example or because
exposures are driving different biology. Let’s examine what omics
features exposures are associated with.

``` r
# Run association testing between every exposure and omics feature
expom <- run_exposure_omics_association(
    exposomicset = expom,
    exposures = exp_vars,
    action = "add"
)
```

Now let’s see how many omics features each exposure is associated with
using the `plot_exposure_omics_association` function where we can either
plot by the individual exposure or exposure category:

``` r
# Plot the number of exposure-omics associations
plot_exposure_omics_association(
  exposomicset = expom,
  plot_type = "category",
  pval_col = "p.value",
  pval_thresh = 0.05)
```

![Barplot of the number of omics features associated with
exposures.](exposure_omics_associations_files/figure-html/plot-exp-omics-assoc-1.png)

Barplot of the number of omics features associated with exposures.

Here we can see that there are omics features associated with exposures,
while there are fewer that are associated with asthma directly.

[Next: Multi-omics
integration](https://BioNomad.github.io/tidyexposomics/articles/multi_omics_integration.md)
