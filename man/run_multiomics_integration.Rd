% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_multiomics_integration.R
\name{run_multiomics_integration}
\alias{run_multiomics_integration}
\title{Perform Multi-Omics Integration}
\usage{
run_multiomics_integration(
  expomicset,
  method = "MCIA",
  n_factors = 10,
  scale = TRUE,
  action = "add"
)
}
\arguments{
\item{expomicset}{A \code{MultiAssayExperiment} object containing at least two omics assays.}

\item{method}{A character string specifying the integration method. Options: \code{"MOFA"}, \code{"MCIA"}. Default is \code{"MCIA"}.}

\item{n_factors}{An integer specifying the number of latent factors/components to compute. Default is \code{10}.}

\item{scale}{A logical indicating whether to standardize (\code{Z-score}) the assays before integration. Default is \code{TRUE}.}

\item{action}{A character string specifying whether to store (\code{"add"}) or return (\code{"get"}) the results. Default is \code{"add"}.}
}
\value{
A \code{MultiAssayExperiment} object with integration results added to metadata (if \code{action = "add"}) or a list with:
\item{method}{The integration method used (\code{"MOFA"} or \code{"MCIA"}).}
\item{result}{The trained integration model.}
}
\description{
Integrates multiple assays in a \code{MultiAssayExperiment} object using either MOFA+ or MCIA.
}
\details{
This function:
\itemize{
\item \strong{MOFA+ (\code{method = "MOFA"})}:
\itemize{
\item Creates a MOFA model from \code{expomicset}.
\item Trains the model using default MOFA parameters.
\item Saves and loads the trained model.
}
\item \strong{MCIA (\code{method = "MCIA"})}:
\itemize{
\item Applies \strong{NIPALS MCIA} via \code{nipalsMCIA::nipals_multiblock()}.
}
\item \strong{Scaling}: If \code{scale = TRUE}, all assays are standardized using \code{.scale_multiassay()}.
\item \strong{Output Handling}:
\itemize{
\item \code{"add"}: Stores results in \code{metadata(expomicset)$integration_results}.
\item \code{"get"}: Returns a list containing \code{method} and \code{result} (trained model).
}
}
}
\examples{
\dontrun{
expom <- run_multiomics_integration(
  expomicset = expom,
  method = "MCIA",
  n_factors = 10,
  scale = TRUE
)
}

}
