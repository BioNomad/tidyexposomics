% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_association.R
\name{run_association}
\alias{run_association}
\title{Generalized Association Analysis Across Data Sources}
\usage{
run_association(
  expomicset,
  outcome,
  source = c("omics", "exposures", "factors", "go_pcs"),
  covariates = NULL,
  feature_set = NULL,
  top_n = NULL,
  family = "gaussian",
  correction_method = "fdr",
  action = "add",
  min_genes = 10,
  feature_col = NULL,
  mirna_assays = NULL,
  print = NULL
)
}
\arguments{
\item{expomicset}{A MultiAssayExperiment object.}

\item{outcome}{A character string naming the outcome variable.}

\item{source}{One of 'omics', 'exposures', 'factors', or 'go_pcs'.}

\item{covariates}{Optional character vector of covariates.}

\item{feature_set}{Optional vector of feature names (e.g., exposures, GO terms).}

\item{top_n}{For omics: select top N variable features per assay.}

\item{family}{GLM family ('gaussian' or 'binomial').}

\item{correction_method}{Method for p-value adjustment. Default is 'fdr'.}

\item{action}{Return results ('get') or attach to metadata ('add').}

\item{min_genes}{Minimum genes per GO group (for GO PC source).}

\item{feature_col}{Optional column to map GO terms to features.}

\item{mirna_assays}{Assays to exclude in GO PC mode.}
}
\value{
If action = 'add', returns updated MultiAssayExperiment. Else returns result list.
}
\description{
Performs GLM-based association testing between an outcome and features
from omics, exposures, latent factors, or GO-based PCs in a MultiAssayExperiment.
}
