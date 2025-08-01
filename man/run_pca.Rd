% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_pca.R
\name{run_pca}
\alias{run_pca}
\title{Perform Principal Component Analysis (PCA)}
\usage{
run_pca(
  expomicset,
  log_trans_exp = FALSE,
  log_trans_omics = TRUE,
  action = "add"
)
}
\arguments{
\item{expomicset}{A \code{MultiAssayExperiment} object containing omics
and exposure data.}

\item{log_trans_exp}{A boolean value specifying whether to log2
transform the exposure data}

\item{log_trans_omics}{a boolean value specifying whether to log2
transform the omics data}

\item{action}{A character string specifying whether to store
(\code{"add"}) or return (\code{"get"}) the results. Default is \code{"add"}.}
}
\value{
A \code{MultiAssayExperiment} object with PCA results added to
metadata (if \code{action = "add"}) or a list with:
\item{pca_df}{A tibble of the transformed input data.}
\item{pca_feature}{A \code{prcomp} object containing PCA results for features.}
\item{pca_sample}{A \code{prcomp} object containing PCA results for samples.}
\item{outliers}{A character vector of detected sample outliers.}
}
\description{
Runs PCA on the feature and sample spaces of a \code{MultiAssayExperiment} object,
identifying outliers based on Mahalanobis distance.
}
\details{
This function:
\itemize{
\item Identifies \strong{common samples} across all assays and exposure data.
\item Performs \strong{PCA on features} (transformed and standardized).
\item Performs \strong{PCA on samples} and computes Mahalanobis distance to
detect outliers.
\item \strong{Output Handling}:
\itemize{
\item \code{"add"}: Stores results in \code{metadata(expomicset)$pca} and
updates \code{colData} with PCs.
\item \code{"get"}: Returns a list containing the PCA results.
}
}
}
\examples{
# create example data
mae <- make_example_data(
    n_samples = 10,
    return_mae = TRUE
)

# run pca
mae <- mae |>
    run_pca()

}
