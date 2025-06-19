# tidyexposomics <a href="#"><img src="./vignettes/logo.png" align="right" height="200" /></a>

<br>
<br>
<br>
<br>

## Overview

The `tidyexposomics` package is designed to facilitate the integration of exposure and omics data to identify exposure-omics associations. We structure our commands to fit into the tidyverse framework, where commands are designed to be simplified and intuitive. Here we provide functionality to perform quality control, sample and exposure association analysis, differential abundance analysis, multi-omics integration, and functional enrichment analysis.

<p align="center" width="100%">
    <img width="80%" src="vignettes/overview.png">
</p>

## Installation

You can install the development version of `tidyexposomics` from GitHub with:

```r
# Install directly through GitHub
devtools::install_github("BioNomad/tidyexposomics")
```

## Command Structure

To make the package more user-friendly, we have named our functions to be more intuitive. For example, we use the following naming conventions:

<p align="center" width="100%">
    <img width="80%" src="vignettes/command_str.png">
</p>

We provide functionality to either `add` results to the existing object storing the omics/exposure data or to return the direct results using the `get` option using the `action` argument. We suggest adding results, as we also include a `run_pipeline_summary()` function to generate a diagram of the workflow. This is useful for keeping track of the pipeline steps. 

## Quick Start

The following code is a great way to get started with the package. It includes loading example data, performing basic quality control, running exposure-wide association studies (ExWAS), differential abundance analysis, correlating differentially expressed genes (DEGs) with exposures, and functional enrichment analysis.

<div class="callout">
<h3>More to tidyexposomics!</h3>
<p>
However, there is so much more to the `tidyexposomics` package! So check out the [Get Started](articles/tidyexposomics.html) page for a more detailed walkthrough of the package's functionality.
</p>
</div>


```R
# Install and Load Packages
# -------------------------

devtools::install_github("BioNomad/tidyexposomics")
library(tidyexposomics)
library(tidyverse)


# Load Example Data

omics_list <- list(
  "Gene Expression" = exp,
  "Metabolomics" = met,
  "Proteomics" = prot
)

fdata <- list(
  "Gene Expression" = exp_fdata,
  "Metabolomics" = met_fdata,
  "Proteomics" = prot_fdata
)

# Create MultiAssayExperiment container
expom <- create_expomicset(
  var_info = des,
  exposure = meta,
  omics = omics_list,
  row_data = fdata
)

# Identify exposure variables
exp_vars <- des |>
  filter(grepl("exposure|chemical", category, ignore.case = TRUE)) |>
  pull(variable) |>
  as.character()

# Basic Quality Control
# ---------------------

expom <- expom |>
  filter_missing(na_thresh = 1) |>
  run_impute_missing(
    exposure_impute_method = "median",
    omics_impute_method = "median"
  ) |>
  run_pca(action = "add") |>
  filter_sample_outliers(outliers = c("s411","s378", "s588", "s764", "s857", "s918", "s936")) |>
  run_normality_check(action = "add") |>
  transform_exposure(transform_method = "boxcox_best", cols_of_interest = exp_vars)

# Exposure-Wide Association (ExWAS)

expom <- expom |>
  associate_exposure_outcome(
    outcome = "hs_asthma",
    exposures = exp_vars[!exp_vars %in% c("hs_asthma", "hs_child_age", "e3_sex")],
    covariates = c("hs_child_age", "e3_sex"),
    family = "binomial",
    action = "add"
  )

# Plot significant associations
expom |>
  plot_associate_exposure_outcome(
    subtitle = "Covariates: Age, Sex",
    result = 1,
    filter_thresh = 0.05
  )


# Differential Abundance Analysis
# -------------------------------

expom <- expom |>
  run_differential_abundance(
    formula = ~ hs_asthma + hs_child_age + e3_sex,
    method = "limma_voom",
    minimum_counts = 1,
    minimum_proportion = 0.1,
    scaling_method = "none",
    action = "add"
  )

# Plot results as volcano plot
expom |>
  plot_volcano(
    top_n_label = 3,
    logFC_thresh = log2(1),
    pval_thresh = 0.05,
    nrow = 1
  )


# Correlate DEGs with Exposures
# -----------------------------

expom <- expom |>
  correlate_exposures_degs(
    exposure_cols = exp_vars,
    robust = TRUE,
    score_thresh = 0.25,
    deg_logfc_thresh = log2(1),
    deg_pval_thresh = 0.05,
    correlation_cutoff = 0.01,
    pval_cutoff = 0.1,
    action = "add"
  )


# Functional Enrichment Analysis
# ------------------------------

expom <- expom |>
  run_enrichment(
    geneset = "deg_exp_cor",
    feature_col = "gene",
    clustering_approach = "dynamic",
    pval_threshold = 0.05,
    pvalueCutoff = 0.1,
    pAdjustMethod = "none",
    qvalueCutoff = 1,
    logfc_threshold = log2(1),
    action = "add"
  )

# Plot enriched terms (dotplot)
expom |>
  plot_dotplot_enrichment(
    geneset = "deg_exp_cor",
    top_n = 5,
    n_per_group = 5,
    add_top_genes = TRUE,
    top_n_genes = 5
  )

```
