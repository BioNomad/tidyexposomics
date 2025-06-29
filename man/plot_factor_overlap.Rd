% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_factor_overlap.R
\name{plot_factor_overlap}
\alias{plot_factor_overlap}
\title{Plot Overlap of Top Factor Features Across Experiments}
\usage{
plot_factor_overlap(
  expomicset,
  venn_text_size = 3,
  venn_stroke_size = 0.1,
  venn_set_name_size = 3.5,
  venn_show_percent = FALSE,
  venn_colors = NULL,
  shared_bar_colors = NULL,
  da_bar_colors = NULL,
  da_bar_facet_cols = NULL,
  da_bar_facet_alpha = 0.5,
  left_heights = c(2, 1),
  col_widths = c(1, 1)
)
}
\arguments{
\item{expomicset}{A \code{MultiAssayExperiment} object with \code{top_factor_features} and \code{common_top_factor_features} in metadata.}

\item{venn_text_size}{Text size for Venn diagram labels. Default is \code{3}.}

\item{venn_stroke_size}{Outline stroke width in Venn diagram. Default is \code{0.1}.}

\item{venn_set_name_size}{Text size for set names in the Venn diagram. Default is \code{3.5}.}

\item{venn_show_percent}{Logical, whether to show percentages in the Venn diagram. Default is \code{FALSE}.}

\item{venn_colors}{Optional named vector of colors for each factor in the Venn diagram.}

\item{shared_bar_colors}{Optional named vector of colors for the shared/unique barplot fill.}

\item{da_bar_colors}{Optional named vector of colors for DEG (differentially abundant) barplots.}

\item{da_bar_facet_cols}{Optional vector of colors for facet strip backgrounds in DEG-related barplots.}

\item{da_bar_facet_alpha}{Numeric from 0–1 for facet background transparency. Default is \code{0.5}.}

\item{left_heights}{A numeric vector of relative heights for the Venn + shared DEG barplot panel. Default is \code{c(2, 1)}.}

\item{col_widths}{A numeric vector of relative widths for left (Venn/shared) vs right (DEG barplot) columns. Default is \code{c(1, 1)}.}
}
\value{
A \code{patchwork} object showing the combined plots: barplots and Venn diagram.
}
\description{
Visualizes the \strong{overlap of top features} across multi-omics integration factors using barplots and a Venn diagram.
}
\details{
This function:
\itemize{
\item Displays overlap of top features across integration factors using a Venn diagram (\code{ggvenn}).
\item Visualizes barplots of:
\itemize{
\item Shared vs unique top features across experiments.
\item Differentially abundant top features (per factor).
\item DEGs among only the shared features (overall).
}
\item Adds optional customization for layout spacing, colors, and transparency.
}
}
\examples{
\dontrun{
plot_factor_overlap(expomicset, venn_show_percent = TRUE)
}

}
