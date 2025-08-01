% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_loading.R
\name{download_dataset}
\alias{download_dataset}
\title{Download and load a zipped .RData file from the tidyexposomics GitHub release}
\usage{
download_dataset(
  name = c("omics_list", "fdata", "meta", "annotated_cb", "chebi", "ecto", "hpo"),
  dest_dir = tempdir(),
  verbose = TRUE
)
}
\arguments{
\item{name}{The name of the dataset (without .zip or .RData extension).
Valid options include: "omics_list", "fdata", "meta", "annotated_cb",
"chebi", "ecto", "hpo".}

\item{dest_dir}{Destination directory to save
and extract the file (default: tempdir()).}

\item{verbose}{Print messages? (default: TRUE)}
}
\value{
Invisibly returns the loaded object.
}
\description{
Download and load a zipped .RData file from the tidyexposomics GitHub release
}
\examples{
# load the example meta data
download_dataset(
    name = "meta"
)

}
