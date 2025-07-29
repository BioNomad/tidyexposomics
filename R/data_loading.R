#' Download and load a zipped .RData file from the tidyexposomics GitHub release
#'
#' @param name The name of the dataset (without .zip or .RData extension).
#'   Valid options include: "omics_list", "fdata", "meta", "annotated_cb", "chebi", "ecto", "hpo".
#' @param dest_dir Destination directory to save and extract the file (default: tempdir()).
#' @param verbose Print messages? (default: TRUE)
#'
#' @return Invisibly returns the loaded object.
#' @export
# download_dataset <- function(name = c("omics_list",
#                                       "fdata",
#                                       "meta",
#                                       "annotated_cb",
#                                       "chebi",
#                                       "ecto",
#                                       "hpo"),
#                              dest_dir = tempdir(),
#                              verbose = TRUE) {
#   name <- match.arg(name)
#   zip_name <- paste0(name, ".zip")
#   rdata_name <- paste0(name, ".RData")
#
#   url <- paste0("https://github.com/bionomad/tidyexposomics/releases/download/data-v1/", zip_name)
#   zip_path <- file.path(dest_dir, zip_name)
#
#   if (verbose) message("Downloading ", zip_name, "...")
#   utils::download.file(url, destfile = zip_path, mode = "wb", quiet = !verbose)
#
#   if (verbose) message("Unzipping to ", dest_dir)
#   utils::unzip(zip_path, exdir = dest_dir)
#
#   if (verbose) message("Loading ", rdata_name)
#   load(file.path(dest_dir, rdata_name), envir = globalenv())
#
#   invisible(get(name, envir = globalenv()))
# }

download_dataset <- function(name = c("omics_list", "fdata", "meta", "annotated_cb", "chebi", "ecto", "hpo"),
                             dest_dir = tempdir(),
                             verbose = TRUE) {
  name <- match.arg(name)
  zip_name <- paste0(name, ".zip")
  rdata_name <- paste0(name, ".RData")

  url <- paste0("https://github.com/bionomad/tidyexposomics/releases/download/data-v1/", zip_name)
  zip_path <- file.path(dest_dir, zip_name)

  if (verbose) message("Downloading ", zip_name, " using curl...")
  curl::curl_download(url, destfile = zip_path, mode = "wb", quiet = !verbose)

  if (verbose) message("Unzipping to ", dest_dir)
  utils::unzip(zip_path, exdir = dest_dir)

  if (verbose) message("Loading ", rdata_name)
  load(file.path(dest_dir, "data", rdata_name), envir = globalenv())

  invisible(get(name, envir = globalenv()))
}


#' Load example omics data for vignette/demo use
#'
#' Downloads and loads `omics_list`, `annotated_cb`, `fdata`, and `meta` into the global environment.
#'
#' @param dest_dir Destination directory to store temporary data files.
#' @export
load_example_data <- function(dest_dir = tempdir()) {
  download_dataset("omics_list", dest_dir)
  download_dataset("annotated_cb", dest_dir)
  download_dataset("fdata", dest_dir)
  download_dataset("meta", dest_dir)
}


#' Load Ontology Data
#'
#' Downloads and loads `chebi`, `ecto`, and `hpo` into the global environment.
#'
#' @param dest_dir Destination directory to store temporary data files.
#' @export
load_annotation_data <- function(dest_dir = tempdir()) {
  download_dataset("chebi", dest_dir)
  download_dataset("ecto", dest_dir)
  download_dataset("hpo", dest_dir)
}
