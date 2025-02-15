library(tidyverse)
library(tidybulk)
library(MultiAssayExperiment)

# Define a function to convert MAE to a tibble
tidyMultiAssayExperiment <- function(mae) {
  assays_as_tibbles <- lapply(names(experiments(mae)), function(experiment_name) {
    experiment <- experiments(mae)[[experiment_name]]
    
    lapply(assayNames(experiment), function(assay_name) {
      assay_data <- assay(experiment, assay_name) |>
        as.data.frame() |>
        rownames_to_column(".feature") |>
        pivot_longer(-.feature, names_to = ".sample", values_to = assay_name)
      
      return(assay_data)
    }) |>
      purrr::reduce(left_join, by = c(".feature", ".sample")) |>
      mutate(experiment_name = experiment_name) |>
      relocate(experiment_name, .after = last_col())
  }) |>
    bind_rows()
  
  col_metadata <- colData(mae) |>
    as.data.frame() |>
    rownames_to_column(".sample")
  
  tidy_mae <- assays_as_tibbles |>
    left_join(col_metadata, by = ".sample")
  
  gc()
  
  return(tidy_mae)
}


# Custom print method for tidyMultiAssayExperiment
print.tidyMultiAssayExperiment <- function(x, ..., n = NULL, width = NULL, n_extra = NULL) {
  require(tibble)
  require(cli)
  require(dplyr)
  
  # Convert to tidy tibble
  my_tibble <- tidyMultiAssayExperiment(x)
  
  # Gather summary details
  number_of_features <- my_tibble |> distinct(.feature) |> nrow()
  number_of_samples <- my_tibble |> distinct(.sample) |> nrow()
  number_of_experiments <- my_tibble |> distinct(experiment_name) |> nrow()
  assay_names <- my_tibble |> distinct(experiment_name) |> pull(experiment_name) |> paste(collapse = ", ")
  
  # Print header
  cat(cli::col_br_black("# A MultiAssayExperiment-tibble abstraction: ",sprintf(
    "%s × %s",
    nrow(my_tibble), ncol(my_tibble)
  ),"\n"))
  cat(cli::col_br_black(sprintf(
    "# Features=%s | Samples=%s | Experiments=%s\nAssays: %s\n",
    number_of_features,
    number_of_samples,
    number_of_experiments,
    assay_names
  )))
  
  # Print the tibble
  print(as_tibble(my_tibble), n = n, width = width, n_extra = n_extra)
  invisible(x)
}

# Overwrite the `show` method for MultiAssayExperiment
setMethod("show", "MultiAssayExperiment", function(object) {
  if (inherits(object, "MultiAssayExperiment")) {
    # Use the custom print method directly to avoid recursion
    print.tidyMultiAssayExperiment(object)
    gc()
  } else {
    # Fallback for other objects
    methods::show(object)
  }
})


