#' Download tidyexposomics zenodo data
#'
#' Download and load a zipped .RData file from the tidyexposomics
#' Zenodo release: https://doi.org/10.5281/zenodo.17049349. The data bundle
#' contains ontology data from ECTO, CHEBI, and HPO, and processed data from
#' the ISGlobal Exposome Data Challenge 2021
#' (Maitre et al., Environment International, 2022;
#' DOI: 10.1016/j.envint.2022.107422). Please cite the original study
#' and data creators when using this bundle.
#'
#'
#' @param name The name of the dataset (without .zip or .RData extension).
#'   Valid options include: "omics_list", "fdata", "meta", "annotated_cb",
#'    "chebi", "ecto", "hpo".
#' @param dest_dir Destination directory to save
#' and extract the file (default: tempdir()).
#' @param verbose Print messages? (default: TRUE)
#'
#' @return Invisibly returns the loaded object.
#' @examples
#' # load the example meta data
#' download_dataset(
#'     name = "meta"
#' )
#'
#' @export
download_dataset <- function(
    name = c(
        "omics_list",
        "fdata",
        "meta",
        "annotated_cb",
        "chebi",
        "ecto",
        "hpo"
    ),
    dest_dir = tempdir(),
    verbose = TRUE) {
    .check_suggested("curl")
    name <- match.arg(name)
    zip_name <- paste0(name, ".zip")
    rdata_name <- paste0(name, ".RData")

    url <- paste0(
        "https://zenodo.org/records/17049350/files/",
        zip_name
    )
    zip_path <- file.path(dest_dir, zip_name)

    if (verbose) message("Downloading ", zip_name, " using curl...")
    curl::curl_download(url,
        destfile = zip_path,
        mode = "wb",
        quiet = !verbose
    )

    if (verbose) message("Unzipping to ", dest_dir)
    utils::unzip(zip_path, exdir = dest_dir)

    if (verbose) message("Loading ", rdata_name)
    load(file.path(dest_dir, "data", rdata_name), envir = globalenv())

    invisible(get(name, envir = globalenv()))
}


#' Load example omics data for vignette/demo use
#'
#' Downloads and loads `omics_list`, `annotated_cb`, `fdata`, and `meta`
#' into the global environment.
#'
#' @param dest_dir Destination directory to store temporary data files.
#' @param to_load Character vector indicating which data set to load.
#'
#' @return Invisibly returns `TRUE` after the example data is loaded into the global environment.
#'
#' @examples
#' # load example data
#' load_example_data(
#'     to_load = "annotated_cb"
#' )
#'
#' @export
load_example_data <- function(
    to_load = c("all", "omics_list", "annotated_cb", "fdata", "meta"),
    dest_dir = tempdir()) {
    to_load <- match.arg(to_load)

    if (to_load == "all") {
        download_dataset("omics_list", dest_dir)
        download_dataset("annotated_cb", dest_dir)
        download_dataset("fdata", dest_dir)
        download_dataset("meta", dest_dir)
    } else {
        download_dataset(to_load, dest_dir)
    }

    # returning a value
    invisible(TRUE)
}


#' Load Ontology Data
#'
#' Downloads and loads `chebi`, `ecto`, and `hpo` into the global environment.
#'
#' @param dest_dir Destination directory to store temporary data files.
#' @param to_load Character vector indicating which ontology to load.
#'
#' @return Invisibly returns `TRUE` after the ontology data is loaded into the global environment.
#'
#' @examples
#' # load ontology
#' load_annotation_data(
#'     to_load = "ecto"
#' )
#'
#' @export
load_annotation_data <- function(
    to_load = c("all", "chebi", "ecto", "hpo"),
    dest_dir = tempdir()) {
    to_load <- match.arg(to_load)

    if (to_load == "all") {
        download_dataset("chebi", dest_dir)
        download_dataset("ecto", dest_dir)
        download_dataset("hpo", dest_dir)
    } else {
        download_dataset(to_load, dest_dir)
    }

    # returning a value
    invisible(TRUE)
}
