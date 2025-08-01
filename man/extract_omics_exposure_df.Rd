% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_omics_exposure_df.R
\name{extract_omics_exposure_df}
\alias{extract_omics_exposure_df}
\title{Extract Merged Omics and Exposure Data Frame}
\usage{
extract_omics_exposure_df(expomicset, variable_map = NULL, log2_trans = TRUE)
}
\arguments{
\item{expomicset}{A \code{MultiAssayExperiment} object containing
omics and exposure data.}

\item{variable_map}{A data frame with columns \code{"variable"} and \code{"exp_name"},
indicating which variables belong to each omics or exposure domain.}

\item{log2_trans}{Logical; whether to log2-transform omics data.
Default is \code{TRUE}.}
}
\value{
A data frame where rows correspond to samples,
and columns contain exposure variables and log2-transformed omics features.
Columns from different omics types are disambiguated using prefixes.
}
\description{
This function extracts and merges exposure variables from \code{colData}
with selected features from omics datasets
in a \code{MultiAssayExperiment} object. Optionally applies log2 transformation
to omics data and restricts features based on a variable map.
}
\details{
If \code{variable_map} is provided, it is used to select variables from
both exposures and omics. If not provided, all numeric \code{colData} variables
are used as exposures (excluding variables matching \verb{^PC}),
and all omics features are included.
}
\examples{
# create example data
mae <- make_example_data(
    n_samples = 10,
    return_mae = TRUE
)
# export the omics exposure df
merged_df <- extract_omics_exposure_df(
    mae,
    log2_trans = TRUE
)

}
