<!-- badges: start -->
  [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
<!-- badges: end -->
  
# tidyexposomics <a href="#"><img src="./vignettes/logo.png" align="right" height="200" /></a>

<br>
<br>
<br>
<br>


- **Website:** https://bionomad.github.io/tidyexposomics/index.html

## Overview

The `tidyexposomics` package is designed to facilitate the integration of exposure and omics data to identify exposure-omics associations. Functions follow the tidyverse framework, where commands are designed to be simplified and intuitive. The tidyexposomics package provides functionality to perform quality control, sample and exposure association analysis, differential abundance analysis, multi-omics integration, and functional enrichment analysis.

<p align="center" width="100%">
    <img width="80%" src="vignettes/overview.png">
</p>


## Command Structure

To make the package more user-friendly, we have named our functions to be more intuitive. For example, we use the following naming conventions:

<p align="center" width="100%">
    <img width="80%" src="vignettes/command_str.png">
</p>

Results can be added to the `MultiAssayExperiment` object or returned directly with `action = 'get'`. We suggest adding results, given that pipeline steps are tracked and can be output to the R console, plotted as a workflow diagram, or exported to an excel worksheet.


## Quick Start

The following code is an example of a basic tidyexposomics workflow. It includes loading example data, performing basic quality control, running exposure-wide association studies (ExWAS), differential abundance analysis, correlating differentially expressed genes (DEGs) with exposures, and functional enrichment analysis. However, there is so much more to the `tidyexposomics` package! So check out the [Get Started](articles/tidyexposomics.html) page for a more detailed walkthrough of the package's functionality.


## Installation

The `tidyexposomics` package depends on R (>= 4.4.0) and can be installed using the following code:

```R
# Install and Load Packages
BiocManager("tidyexposomics")
library(tidyexposomics)
library(tidyverse)
```

## Load Example Data

We provide example data based off the [ISGlobal Exposome data challenge 2021](https://www.sciencedirect.com/science/article/pii/S016041202200349X?via%3Dihub). Here, we will examine how exposures and omics features relate to asthma status.

```R
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
```

## Quality Control

We provide several quality control functions including those that handle filtering missing data, imputation, variable normality checks, and variable transformation.

```R
# Filter & impute exposures
expom <- expom |>
  filter_missing(na_thresh = 5) |>
  run_impute_missing(exposure_impute_method = "missforest")

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

## ExWAS 

Here we model the association between exposures and asthma status and adjust our model for child age, biological sex, and cohort.

```R
# Perform ExWAS Analysis
# Perform ExWAS Analysis
expom <- run_association(
    exposomicset = expom,
    source = "exposures",
    outcome = "hs_asthma",
    feature_set = exp_vars,
    action = "add",
    family = "binomial"
)

# Visualize associations
plot_association(
    exposomicset = expom,
    source = "exposures",
    terms = exp_vars,
    filter_thresh = 0.05,
    filter_col = "p.value",
    r2_col = "r2"
)
```

![](vignettes/exwas_assoc.png)

## Differential Abundance

Differentially abundance analysis is supported in `tidyexposomics`. Here we use `limma_trend` to identify features associated with asthma status.

```R
# Run differential abundance analysis
expom <- run_differential_abundance(
    exposomicset = expom,
    formula = ~hs_asthma,
    method = "limma_trend",
    scaling_method = "none",
    action = "add"
)
    
# Plot Differential Abundance Results
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

![](vignettes/volcano.png)

## Multi-Omics Integration

Multi-omics integration is supported to derive insights across omics layers. Here we use the `DIABLO` method and set the outcome variable of interest to asthma status.

```R
# Perform multi-omics integration
expom <- run_multiomics_integration(
    exposomicset = expom,
    method = "DIABLO",
    n_factors = 5,
    outcome = "hs_asthma",
    action = "add"
)
                             
# Identify factors that are associated with the outcome
expom <- run_association(
    exposomicset = expom,
    source = "factors",
    outcome = "hs_asthma",
    feature_set = exp_vars,
    action = "add",
    family = "binomial"
)
    
# Extract top features that contribute to a factor
expom <- extract_top_factor_features(
    exposomicset = expom,
    method = "percentile",
    pval_col = "p_adjust",
    pval_thresh = 0.05,
    percentile = 0.7,
    action = "add"
)
```


## Exposure-Omics Association

Now that we have our multi-omics features associated with asthma status, we can correlate these with our exposures. This helps identify how exposure classes may affect asthma biology.

```R
# Grab top common factor features and ensure
# feature is renamed to variable for the variable_map
top_factor_features <- expom |>
    extract_results(result = "multiomics_integration") |>
    pluck("top_factor_features") |>
    dplyr::select(variable = feature, exp_name)

# Correlate top factor features with exposures
# Perform correlation analysis between factor features and exposures
expom <- run_correlation(
    exposomicset = expom,
    feature_type = "omics",
    variable_map = top_factor_features,
    exposure_cols = exp_vars,
    action = "add",
    correlation_cutoff = 0.3,
    pval_cutoff = 0.05,
    cor_pval_column = "p.value"
)

# Perform correlation analysis between factor features
expom <- run_correlation(
    exposomicset = expom,
    feature_type = "omics",
    variable_map = top_factor_features,
    exposure_cols = exp_vars,
    feature_cors = TRUE,
    action = "add",
    correlation_cutoff = 0.3,
    pval_cutoff = 0.05,
    cor_pval_column = "p.value"
)
```


## Enrichment Analysis

After identifying features associated with asthma and exposures, we can perform functional enrichment analysis to understand what biological processes are affected.

```R
# Run enrichment analysis on factor features correlated with exposures
expom <- run_enrichment(
    exposomicset = expom,
    feature_type = c("omics_cor"),
    feature_col = "feature_clean",
    db = c("GO"),
    species = "goa_human",
    fenr_col = "gene_symbol",
    padj_method = "none",
    pval_thresh = 0.1,
    min_set = 1,
    max_set = 800,
    clustering_approach = "diana",
    action = "add"
)

# Plot enrichment term network plot
plot_enrichment(
    exposomicset = expom,
    feature_type = "omics_cor",
    plot_type = "network",
    label_top_n = 1
)
```

![](vignettes/enr_network.png)
