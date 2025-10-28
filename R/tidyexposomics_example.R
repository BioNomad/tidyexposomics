#' Example exposome multi-omics dataset
#'
#' A downsampled version of the ISGlobal Exposome Data Challenge 2021
#' dataset (Maitre *et al.*, *Environment International*, 2022; DOI: 10.1016/j.envint.2022.107422).
#' This example dataset is included with **tidyexposomics** for vignette
#' evaluation. This subset represents samples from children with
#' low socioeconomic status in the first cohort of the original dataset.
#'
#' @details
#' The data have been filtered and transformed to illustrate
#' the `tidyexposomics` workflow. Only a small subset of variables, and
#' the top 500 most variable features per omic layer are retained.
#'
#' **Contents**
#' \describe{
#'   \item{meta}{A data frame of selected exposure, demographic, and
#'   outcome variables for a subset of participants.}
#'   \item{annotated_cb}{A data frame providing ontology-linked annotation
#'   for the exposures in `meta`.}
#'   \item{exp_filt}{A numeric matrix of gene expression values
#'   (500 features, 48 samples) representing the top-variance transcripts.}
#'   \item{exp_fdata}{Feature-level metadata for `exp_filt`, including
#'   cleaned gene symbols.}
#'   \item{methyl_filt}{A numeric matrix of DNA methylation M-values
#'                      (500 CpG sites, 48 samples).}
#'   \item{methyl_fdata}{Feature-level metadata for `methyl_filt`.}
#' }
#'
#' **Source**
#'
#' Derived from the ISGlobal Exposome Data Challenge 2021
#' (Maitre *et al.*, *Environment International*, 2022;
#' DOI: 10.1016/j.envint.2022.107422),licensed under CC-BY 4.0.
#' The original data are available on Figshare (Project 98813)
#' and GitHub (`isglobal-exposomeHub/ExposomeDataChallenge2021`).
#' This example dataset was processed and downsampled by the
#' tidyexposomics authors and is not a replacement for the full dataset.
#'
#' @docType data
#' @usage data("tidyexposomics_example")
#'
#' @examples
#' data("tidyexposomics_example")
#'
"tidyexposomics_example"
