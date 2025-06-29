% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_exposome_score.R
\name{run_exposome_score}
\alias{run_exposome_score}
\title{Calculate Exposome Score for a MultiAssayExperiment}
\usage{
run_exposome_score(
  expomicset,
  score_type,
  exposure_cols = NULL,
  scale = TRUE,
  score_column_name = NULL
)
}
\arguments{
\item{expomicset}{A \code{MultiAssayExperiment} object containing exposure data in its \code{colData}.}

\item{score_type}{Character string indicating the scoring method to use. One of \code{"mean"}, \code{"median"}, \code{"pca"}, \code{"irt"}, \code{"quantile"}, or \code{"var"}.}

\item{exposure_cols}{Optional character vector specifying which exposure columns to use. If \code{NULL}, all numeric columns (excluding those starting with "PC") are used.}

\item{scale}{Logical; if \code{TRUE}, standardizes (z-scores) the selected exposures before scoring. Default is \code{TRUE}.}

\item{score_column_name}{Optional string to rename the resulting exposome score column in the output \code{colData}. If \code{NULL}, a default name based on \code{score_type} is used.}
}
\value{
A \code{MultiAssayExperiment} object with a new column added to its \code{colData}, containing the computed exposome score.
}
\description{
This function calculates exposome scores using various methods (mean, median, PCA, IRT, quantile, variance)
based on numeric exposure variables found in the \code{colData} of a \code{MultiAssayExperiment} object.
}
\details{
\itemize{
\item \code{"mean"}: Computes the row-wise mean across exposure values.
\item \code{"sum"}: Computes the row-wise sum across exposure values.
\item \code{"median"}: Computes the row-wise median across exposure values.
\item \code{"pca"}: Performs PCA on exposures and uses the first principal component.
\item \code{"irt"}: Fits a unidimensional graded IRT model using decile-binned exposures.
\item \code{"quantile"}: Converts each exposure to deciles (1–10) and sums them row-wise.
\item \code{"var"}: Computes the row-wise variance across exposure values.
}

Note: IRT-based scoring requires the \code{mirt} package.
}
\examples{
# Compute a PCA-based exposome score using specific exposure columns
# mae <- run_exposome_score(mae, exposure_cols = c("no2", "pm25"), score_type = "pca", score_column_name = "pca_score")

}
