#' Download and cache a tidyexposomics dataset
#'
#'
#' @param name Dataset name: one of "omics_list", "fdata", "meta",
#'   "annotated_cb", "expom_1", "chebi", "ecto", "hpo".
#' @param verbose Logical; print messages.
#' @return A list or object loaded from the cached .RData.
#' @importFrom BiocFileCache BiocFileCache bfcnew bfcrpath
#' @export
download_dataset <- function(
  name = c(
      "omics_list", "fdata", "meta", "annotated_cb",
      "expom_1", "chebi", "ecto", "hpo"
  ),
  verbose = TRUE
) {
    name <- match.arg(name)
    zip_name <- paste0(name, ".zip")
    url <- paste0("https://zenodo.org/records/17137560/files/", zip_name)

    # using BiocFileCache now
    bfc <- BiocFileCache::BiocFileCache()
    cached_zip <- BiocFileCache::bfcrpath(bfc, url)

    dest_dir <- tempfile(pattern = paste0("tidyexposomics_", name, "_"))
    utils::unzip(cached_zip, exdir = dest_dir)

    rdata_path <- file.path(dest_dir, "data", paste0(name, ".RData"))
    env <- new.env()
    load(rdata_path, envir = env)
    get(name, envir = env)
}


#' Load Ontology Data
#'
#' Downloads and loads `chebi`, `ecto`, and `hpo` into the global environment.
#'
#' @param to_load Character vector indicating which ontology to load.
#' @param verbose Logical; print messages.
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
  verbose = TRUE
) {
    to_load <- match.arg(to_load)
    names <- if (to_load == "all") c("chebi", "ecto", "hpo") else to_load

    res <- lapply(names, \(nm) download_dataset(nm, verbose = verbose))
    names(res) <- names
    res
}
