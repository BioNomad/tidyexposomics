% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_correlation_tile.R
\name{plot_correlation_tile}
\alias{plot_correlation_tile}
\title{Plot Correlation Tilemap}
\usage{
plot_correlation_tile(
  expomicset,
  feature_type = c("pcs", "degs", "omics", "factors", "factor_features", "exposures"),
  pval_cutoff = 0.05,
  na_color = "grey100",
  fill_limits = c(-1, 1),
  midpoint = 0
)
}
\arguments{
\item{expomicset}{A \code{MultiAssayExperiment} object containing correlation
results in metadata.}

\item{feature_type}{Type of correlation results to plot.
One of \code{"pcs"}, \code{"degs"}, \code{"omics"},
\code{"factors"}, \code{"factor_features"}, or \code{"exposures"}.
Must match the key used in
\code{metadata(expomicset)$correlation[[feature_type]]}.}

\item{pval_cutoff}{Numeric p-value cutoff below which
correlations are displayed.
Nonsignificant tiles are rendered in the \code{na_color}.
Default is \code{0.05}.}

\item{na_color}{Color used to represent nonsignificant or
missing correlations. Default is \code{"grey100"}.}

\item{fill_limits}{Numeric vector of length 2 defining the scale
limits for correlation values. Default is \code{c(-1, 1)}.}

\item{midpoint}{Numeric value for centering the fill gradient.
Default is \code{0}.}
}
\value{
A \code{ggplot2} tile plot (or a combined tile + barplot if
\code{feature_type = "pcs"}).
}
\description{
Visualizes a correlation matrix as a heatmap tile plot using
correlation results stored in the metadata of a
\code{MultiAssayExperiment} object. When \code{feature_type = "pcs"},
the function forces PCs to appear on the x-axis and exposures on the y-axis,
and it adds a barplot showing how many PCs are significantly
associated with each exposure. It also suppresses nonsignificant tiles
based on a specified p-value threshold.
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

# correlate with exposures
mae <- mae |>
    run_correlation(
        feature_type = "pcs",
        exposure_cols = c("exposure_pm25", "exposure_no2", "age", "bmi")
    )

# make the correlation tile plot
cor_tile_p <- mae |>
    plot_correlation_tile(
        feature_type = "pcs"
    )

}
