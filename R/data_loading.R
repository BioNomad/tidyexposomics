#' Download and cache a tidyexposomics dataset
#'
#' @param name Dataset name: one of "omics_list", "fdata", "meta",
#'   "annotated_cb", "expom_1", "chebi", "ecto", "hpo".
#' @param verbose Logical; print messages.
#' @param validate Logical; validate MD5 checksum.
#' @return A list or object loaded from the cached .RData.
#' @importFrom BiocFileCache BiocFileCache bfcrpath
#' @importFrom tools md5sum
#' @export
download_dataset <- function(
  name = c(
      "omics_list", "fdata", "meta", "annotated_cb",
      "expom_1", "chebi", "ecto", "hpo"
  ),
  verbose = TRUE,
  validate = TRUE
) {
    name <- match.arg(name)

    checksums <- c(
        omics_list = "0e0b4bc175da647363478954039b2b9d",
        fdata = "561e9eb4136e4d6116914adb9fc5a30f",
        meta = "db317363ba547eaeb6543acca850ac92",
        annotated_cb = "5ca6a27a29ea50004b1e8dcff50a5487",
        expom_1 = "5e84c508469a4db77a22f03f27ab9ccd",
        chebi = "8a5703b5ac54e63006c002c834a3f703",
        ecto = "0a267757ea47d5d99f5f4929f26186c6",
        hpo = "9948e6b99815da81daa2d9a2e586a9ae"
    )

    zip_name <- paste0(name, ".zip")
    url <- paste0("https://zenodo.org/records/17137560/files/", zip_name)

    bfc <- BiocFileCache::BiocFileCache()
    cached_zip <- BiocFileCache::bfcrpath(bfc, url)

    if (validate) {
        observed <- tools::md5sum(cached_zip)
        expected <- checksums[[name]]
        if (observed != expected) {
            stop(
                "Checksum mismatch for ", name, ".\n",
                "Expected: ", expected, "\n",
                "Observed: ", observed,
                call. = FALSE
            )
        }
        if (verbose) message("Checksum validated for ", name)
    }

    dest_dir <- tempfile(pattern = paste0("tidyexposomics_", name, "_"))
    utils::unzip(cached_zip, exdir = dest_dir)
    rdata_path <- file.path(dest_dir, "data", paste0(name, ".RData"))

    env <- new.env()
    load(rdata_path, envir = env)
    get(name, envir = env)
}

#' Load Ontology Data
#'
#' Downloads and loads `chebi`, `ecto`, and `hpo`.
#'
#' @param to_load Character vector indicating which ontology to load.
#' @param verbose Logical; print messages.
#' @param validate Logical; validate MD5 checksum.
#'
#' @return A named list of ontology objects.
#'
#' @examples
#' \dontrun{
#' # load single ontology
#' onts <- load_annotation_data(to_load = "ecto")
#'
#' # load all ontologies
#' onts <- load_annotation_data(to_load = "all")
#' }
#'
#' @export
load_annotation_data <- function(
  to_load = c("all", "chebi", "ecto", "hpo"),
  verbose = TRUE,
  validate = TRUE
) {
    to_load <- match.arg(to_load)
    names <- if (to_load == "all") c("chebi", "ecto", "hpo") else to_load
    res <- lapply(names, \(nm) download_dataset(nm, verbose = verbose, validate = validate))
    names(res) <- names
    res
}
