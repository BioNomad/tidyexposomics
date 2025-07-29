#' @title Export tidyexposomics Results to Excel
#'
#' @description
#' Exports selected results stored in a `MultiAssayExperiment` object created by the `tidyexposomics` pipeline to an Excel workbook. Users can select which result types to include, and optionally add placeholder sheets for missing data.
#'
#' @param expom A `MultiAssayExperiment` object with results stored in `@metadata`, typically created by the `tidyexposomics` pipeline.
#' @param file Character. Path to the output Excel file (e.g., `"tidyexposomics_results.xlsx"`).
#' @param result_types Character vector specifying which result categories to export. Options include:
#'   - `"correlation"`: Correlation results (e.g., exposure–PC, exposure–feature).
#'   - `"association"`: Association results (e.g., ExWAS, omics–outcome models).
#'   - `"differential_analysis"`: Differential abundance results, including sensitivity analysis if available.
#'   - `"multiomics_integration"`: Common top features contributing to latent factors from multi-omics integration.
#'   - `"network"`: Exposure impact metrics from network analyses.
#'   - `"enrichment"`: Enrichment results by omic and exposure category.
#'   - `"exposure_summary"`: Summary statistics for exposure variables.
#'   - `"pipeline"`: Overview of steps completed in the pipeline.
#'
#'   Use `"all"` to export all of the above categories.
#'
#' @param include_empty_tabs Logical. If `TRUE`, adds placeholder sheets for any missing result types. Default is `FALSE`.
#'
#' @return An Excel file is written to the specified path. A message is printed with the file location.
#'
#' @examples
#' \dontrun{
#' extract_results_excel(expom_1, file = "tidy_results.xlsx", result_types = "all")
#' extract_results_excel(expom_1, file = "enrichment_only.xlsx", result_types = "enrichment")
#' }
#'
#' @importFrom openxlsx createWorkbook addWorksheet writeData addStyle setColWidths saveWorkbook createStyle
#' @importFrom purrr iwalk
#' @importFrom tibble enframe
#' @export
extract_results_excel <- function(expom,
                                  file = "tidyexposomics_results.xlsx",
                                  result_types = c("correlation",
                                                   "association",
                                                   "differential_analysis",
                                                   "multiomics_integration",
                                                   "network",
                                                   "enrichment",
                                                   "exposure_summary",
                                                   "pipeline"),
                                  include_empty_tabs = FALSE) {
  stopifnot("MultiAssayExperiment" %in% class(expom))
  library(openxlsx)

  # Handle "all"
  all_types <- c("correlation",
                 "association",
                 "differential_analysis",
                 "multiomics_integration",
                 "network",
                 "enrichment",
                 "exposure_summary",
                 "pipeline")
  if ("all" %in% result_types) {
    result_types <- all_types
  }

  wb <- createWorkbook()

  safe_add_sheet <- function(name, data) {
    if (!is.null(data) && is.data.frame(data)) {
      addWorksheet(wb, name)
      writeData(wb, name, data, withFilter = TRUE)

      # Bold headers + border around header row
      header_style <- createStyle(
        textDecoration = "bold",
        border = "TopBottomLeftRight",
        borderStyle = "thin"
      )
      addStyle(wb, sheet = name, style = header_style, rows = 1, cols = 1:ncol(data), gridExpand = TRUE)

      # Auto column widths
      setColWidths(wb, sheet = name, cols = 1:ncol(data), widths = "auto")

    } else if (include_empty_tabs) {
      addWorksheet(wb, name)
      writeData(wb, name, data.frame(Note = "No data available."))
    }
  }

  # --- Adding Correlation Results ---------
  if ("correlation" %in% result_types) {
    message("Writing Correlation Results.")
    enr <- expom@metadata$correlation
    if (!is.null(enr) && is.list(enr)) {
      purrr::iwalk(enr, function(obj, groupname) {
        df <- obj
        safe_add_sheet(paste0("Correlation_", groupname), df)
      })
    }
  }

  # --- Adding Association Results -----

  if ("association" %in% result_types) {
    message("Writing Association Results.")
    enr <- expom@metadata$association
    if (!is.null(enr) && is.list(enr)) {
      purrr::iwalk(enr, function(obj, groupname) {
        df <- obj
        safe_add_sheet(paste0("association_", groupname), df)
      })
    }
  }

  # --- Adding Differential Abundance Results --------
  if ("differential_analysis" %in% result_types) {
    message("Writing Differential Abundance Results.")
    da <- expom@metadata$differential_analysis$differential_abundance
    if (!is.null(da)) {
      safe_add_sheet("Differential Abundance", da)
    }

    if ("sensitivity_analysis" %in% names(expom@metadata$differential_analysis)){
      message("Writing Sensitivity Analysis Results.")
      sens <- expom@metadata$differential_analysis$sensitivity_analysis$feature_stability
      if (!is.null(sens)) {
        safe_add_sheet("Sensitivity Analysis", sens)
      }
    }
  }

  # --- Adding Multiomics Integration Results ----------
  if ("multiomics_integration" %in% result_types) {
    message("Writing Multiomics Integration Results.")
    factors <- expom@metadata$multiomics_integration$common_top_factors
    if (!is.null(factors)) {
      safe_add_sheet("Common Factor Features", factors)
    }
  }

  # --- Adding Network Results --------------
  if ("network" %in% result_types) {
    message("Writing Network Impact Results.")
    exposure_impact <- expom@metadata$network$exposure_impact
    if (!is.null(exposure_impact) && is.list(exposure_impact)) {
      purrr::iwalk(exposure_impact, function(obj, groupname) {
        df <- obj
        safe_add_sheet(paste0("Network_Impact_", groupname), df)
      })
    }
  }

  # --- Adding Enrichment Results -------------
  if ("enrichment" %in% result_types) {
    message("Writing Enrichment Results.")
    enr <- expom@metadata$enrichment
    if (!is.null(enr) && is.list(enr)) {
      purrr::iwalk(enr, function(obj, groupname) {
        df <- obj
        safe_add_sheet(paste0("Enrichment_", groupname), df)
      })
    }
  }

  # --- Adding Pipeline summary -------
  if ("pipeline" %in% result_types) {
    message("Writing Pipeline Summary.")
    pipeline_info <- expom@metadata$summary$steps
    if (is.list(pipeline_info)) {
      steps_df <- tibble::enframe(pipeline_info, name = "Step", value = "Details")
      safe_add_sheet("Pipeline_Steps", steps_df)
    }
  }

  # --- Adding Exposure summaries -------
  if ("exposure_summary" %in% result_types) {
    message("Writing Exposure Summary Results.")
    safe_add_sheet("Exposure_Summary",
                   tryCatch({
                     run_summarize_exposures(expom, action = "get")
                   }, error = function(e) NULL))
  }

  # Save
  saveWorkbook(wb, file, overwrite = TRUE)
  message("Results written to: ", normalizePath(file))
}
