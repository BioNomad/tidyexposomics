% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transform_exposure.R
\name{transform_exposure}
\alias{transform_exposure}
\title{Transform Exposure Variables in a MultiAssayExperiment}
\usage{
transform_exposure(expomicset, exposure_cols = NULL, transform_method = "best")
}
\arguments{
\item{expomicset}{A \code{MultiAssayExperiment} object.}

\item{exposure_cols}{A character vector of column names from \code{colData(expomicset)} to transform.
If \code{NULL}, the function will attempt to use metadata from \code{check_normality()}.}

\item{transform_method}{One of \code{"none"}, \code{"log2"}, \code{"x_1_3"}, \code{"sqrt"},
\code{"boxcox_best"}, or \code{"best"} (default). \code{"best"} evaluates all available
methods and chooses the one that most improves normality.}
}
\value{
A \code{MultiAssayExperiment} object with updated \code{colData} and transformation
results saved in \code{metadata(expomicset)$transformation}.
}
\description{
Applies various transformation methods to numeric exposure variables in the
\code{colData} of a \code{MultiAssayExperiment} object to improve normality.
}
\details{
Supported transformations include:
\itemize{
\item \code{"none"}: No transformation
\item \code{"log2"}: Log base 2 with minimum value shift
\item \code{"x_1_3"}: Cube root with minimum value shift
\item \code{"sqrt"}: Square root with minimum value shift
\item \code{"boxcox_best"}: Box-Cox power transformation based on optimal lambda
\item \code{"best"}: Automatically selects the transformation method that best
normalizes the data (based on Shapiro-Wilk test)
}

The \code{"boxcox_best"} transformation is based on the Box-Cox method from the
\code{MASS} package. Variables are shifted to ensure strictly positive values before
transformation. Lambda intervals are interpreted as follows:
\itemize{
\item \eqn{\lambda < -1.5}: \eqn{1 / x^2}
\item \eqn{-1.5 \leq \lambda < -0.75}: \eqn{1 / x}
\item \eqn{-0.75 \leq \lambda < -0.25}: \eqn{1 / sqrt(x)}
\item \eqn{-0.25 \leq \lambda < 0.25}: \eqn{log(x)}
\item \eqn{0.25 \leq \lambda < 0.75}: \eqn{sqrt(x)}
\item \eqn{0.75 \leq \lambda < 1.5}: identity
\item \eqn{\lambda \geq 1.5}: \eqn{x^2}
}
}
\examples{
\dontrun{
transformed_mae <- transform_exposure(my_mae,
                                      exposure_cols = c("pm25", "no2"),
                                      transform_method = "best")
}
}
