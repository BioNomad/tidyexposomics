## ----setup,include=FALSE------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load-libraries,message=FALSE,warning=FALSE-------------------------------
# Load Libraries
library(tidyverse)
library(tidyexposomics)

## ----load-data----------------------------------------------------------------
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

## ----exp-vars-----------------------------------------------------------------
# Grab exposure variables
exp_vars <- tidyexposomics_example$annotated_cb |>
    filter(category %in% c(
        "aerosol",
        "main group molecular entity",
        "polyatomic entity"
    )) |>
    pull(variable) |>
    as.character()

## ----impute-missing-----------------------------------------------------------

# Impute missing values
expom <- run_impute_missing(
    exposomicset = expom,
    exposure_impute_method = "missforest",
    exposure_cols = exp_vars
)

## ----trasform-vars------------------------------------------------------------

# Transform variables
expom <- transform_exposure(
    exposomicset = expom,
    transform_method = "boxcox_best",
    exposure_cols = exp_vars
)

## ----calc-exposome, results='hide',eval=requireNamespace("mirt", quietly = TRUE)----
# 
# # determine which aerosol variables to use
# aerosols <- c("h_pm25_ratio_preg_None", "h_pm10_ratio_preg_None")
# 
# # Create exposome scores
# expom <- expom |>
#     run_exposome_score(
#         exposure_cols = aerosols,
#         score_type = "median",
#         score_column_name = "exposome_median_score"
#     ) |>
#     run_exposome_score(
#         exposure_cols = aerosols,
#         score_type = "pca",
#         score_column_name = "exposome_pca_score"
#     ) |>
#     run_exposome_score(
#         exposure_cols = aerosols,
#         score_type = "irt",
#         score_column_name = "exposome_irt_score"
#     ) |>
#     run_exposome_score(
#         exposure_cols = aerosols,
#         score_type = "quantile",
#         score_column_name = "exposome_quantile_score"
#     ) |>
#     run_exposome_score(
#         exposure_cols = aerosols,
#         score_type = "var",
#         score_column_name = "exposome_var_score"
#     )

