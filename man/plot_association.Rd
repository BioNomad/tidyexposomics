% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_association.R
\name{plot_association}
\alias{plot_association}
\title{Plot Association Results (Unified Forest Plot)}
\usage{
plot_association(
  expomicset,
  source = c("omics", "exposures", "factors", "go_pcs"),
  terms = NULL,
  filter_col = "p.value",
  filter_thresh = 0.05,
  direction_filter = "all",
  add_r2_tile = TRUE,
  r2_col = "adj_r2",
  facet = FALSE,
  nrow = 1,
  subtitle = NULL
)
}
\arguments{
\item{expomicset}{A \code{MultiAssayExperiment} object containing association
results in metadata.}

\item{source}{Character string indicating the association source.
One of \code{"omics"}, \code{"exposures"}, \code{"factors"}, or \code{"go_pcs"}.}

\item{terms}{Optional character vector of term names to subset the plot to.
Default is \code{NULL} (include all).}

\item{filter_col}{Column used to assess statistical significance
(default: \code{"p.value"}).}

\item{filter_thresh}{Numeric threshold applied to \code{filter_col}
(default: \code{0.05}).}

\item{direction_filter}{Direction of associations to retain.
One of \code{"all"} (default), \code{"up"}, or \code{"down"}.}

\item{add_r2_tile}{Logical; if \code{TRUE}, includes a tile plot for \code{r2_col}
(default: \code{TRUE}).}

\item{r2_col}{Column used for coloring the tile plot (default: \code{"adj_r2"}).}

\item{facet}{Logical; if \code{TRUE} and \code{source == "go_pcs"},
apply nested faceting by experiment and GO cluster (default: \code{FALSE}).}

\item{nrow}{Integer; number of rows for facet layout if enabled (default: \code{1}).}

\item{subtitle}{Optional subtitle for the plot. If \code{NULL},
automatically generated from covariates used in the model.}
}
\value{
A \code{ggplot2} object: either a single forest plot or a composite plot
with an R^2 tile strip.
}
\description{
Generates a forest plot for association results from any source stored
in the metadata of a \code{MultiAssayExperiment} object.
Supports faceting and visual augmentation with R^2 tiles when available.
}
\details{
This function visualizes effect size estimates and confidence intervals from
association analyses. It allows filtering by
direction (\code{"up"} for positive, \code{"down"} for negative) and by significance.
For \code{source = "go_pcs"}, it supports special
formatting by splitting term labels into nested facets.

The R^2 tile (if enabled) adds a side heatmap indicating model fit for
each association. This can be useful for model diagnostics.
}
\examples{

# create example data
mae <- make_example_data(
    n_samples = 10,
    return_mae = TRUE
)

# run association tests
mae <- mae |>
    run_association(
        source = "exposures",
        feature_set = c("exposure_pm25", "exposure_no2"),
        outcome = "smoker",
        covariates = c("age"),
        family = "binomial"
    )

assoc_plot <- mae |>
    plot_association(
        source = "exposures"
    )

}
