% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_loading.R
\name{load_annotation_data}
\alias{load_annotation_data}
\title{Load Ontology Data}
\usage{
load_annotation_data(
  to_load = c("all", "chebi", "ecto", "hpo"),
  dest_dir = tempdir()
)
}
\arguments{
\item{to_load}{Character vector indicating which ontology to load.}

\item{dest_dir}{Destination directory to store temporary data files.}
}
\value{
Invisibly returns \code{TRUE} after the ontology data is loaded into the global environment.
}
\description{
Downloads and loads \code{chebi}, \code{ecto}, and \code{hpo} into the global environment.
}
\examples{
# load ontology
load_annotation_data(
    to_load = "ecto"
)

}
