# Omics Preprocessing

[Previous: Exposome
Scores](https://BioNomad.github.io/tidyexposomics/articles/exposome_scores.md)

**Objectives**

- Learn how to add in omics data to the `exposomicset` object.
- Filter omics data to remove features with low expression/abundance or
  variance.

## Adding in the omics data

We covered [exposure
preprocessing](https://BioNomad.github.io/tidyexposomics/articles/exposure_preprocessing.md)
in a previous vignette. However, this time, we are going to add in the
omics data too!

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
```

Here you will see we have two lists: (1) the omics matrices, where the
rows are omics features (e.g. genes, proteins, metabolites) and the
columns are samples, (2) row data with information about these omics
features. Now we will process our exposure data just like we did in the
[exposure preprocessing
tutorial](https://BioNomad.github.io/tidyexposomics/articles/exposure_preprocessing.md)

``` r
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

## Filtering Omics Features

We can filter omics features based on variance or expression levels. The
`filter_omics` function is used to filter omics features. The method
argument is used to set the method for filtering. Here we can use
either:

- **Variance**: Filters features based on variance. We recommend this
  for omics based on continuous measurements, such as log-transformed
  counts, M-values, protein intensities, or metabolite concentrations.

- **Expression**: Filters features based on expression levels. We
  recommend this for omics where many values may be near-zero or zero,
  such as RNA-seq data.

The `assays` argument is used to specify which omics data to filter. The
`assay_name` argument is used to specify which assay to filter. The
`min_var`, `min_value`, and `min_prop` arguments are used to set the
minimum variance, minimum expression value, and minimum proportion of
samples exceeding the minimum value, respectively.

``` r
# filter omics layers by variance and expression
# Methylation filtering
expom <- filter_omics(
    exposomicset = expom,
    method = "variance",
    assays = "Methylation",
    assay_name = 1,
    min_var = 0.05
)

# Gene expression filtering
expom <- filter_omics(
    exposomicset = expom,
    method = "expression",
    assays = "Gene Expression",
    assay_name = 1,
    min_value = 1,
    min_prop = 0.3
)
```

[Next: Exposure-Omics
Associations](https://BioNomad.github.io/tidyexposomics/articles/exposure_omics_associations.md)
