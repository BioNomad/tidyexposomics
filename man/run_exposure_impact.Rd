% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_exposure_impact.R
\name{run_exposure_impact}
\alias{run_exposure_impact}
\title{Calculate Exposure Impact from Feature-Exposure Correlation Networks}
\usage{
run_exposure_impact(
  expomicset,
  feature_type = c("degs", "omics", "factors"),
  pval_col = "adj.P.Val",
  pval_thresh = 0.1,
  action = c("add", "get")
)
}
\arguments{
\item{expomicset}{A \code{MultiAssayExperiment} object with correlation and network metadata.}

\item{feature_type}{One of \code{"degs"}, \code{"omics"}, or \code{"factors"}.}

\item{pval_col}{Column in differential abundance results to filter DEGs. Default = \code{"adj.P.Val"}.}

\item{pval_thresh}{DEG p-value threshold. Ignored unless \code{feature_type == "degs"}.}

\item{action}{Either \code{"add"} (store in metadata) or \code{"get"} (return list).}
}
\value{
Either an updated MultiAssayExperiment (if action = "add") or a list.
}
\description{
Generalized centrality-based exposure impact analysis using DEG, omics, or factor features.
}
