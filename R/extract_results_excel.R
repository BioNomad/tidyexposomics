#' @title Export tidyexposomics Results to Excel
#'
#' @description
#' Exports selected results stored in a `MultiAssayExperiment` object created
#' by the `tidyexposomics` pipeline to an Excel workbook.
#' Users can select which result types to include,
#' and optionally add placeholder sheets for missing data.
#'
#' @param exposomicset A `MultiAssayExperiment` object with results stored
#' in `@metadata`, typically created by the `tidyexposomics` pipeline.
#' @param file Character. Path to the output Excel file.
#' @param result_types Character vector specifying which result categories
#'  to export. Options include:
#'   - `"correlation"`: Correlation results.
#'   - `"association"`: Association results.
#'   - `"differential_analysis"`: Differential abundance results,
#'    including sensitivity analysis if available.
#'   - `"multiomics_integration"`: Common top features contributing to
#'    latent factors from multi-omics integration.
#'   - `"network"`: Exposure impact metrics from network analyses.
#'   - `"enrichment"`: Enrichment results by omic and exposure category.
#'   - `"exposure_summary"`: Summary statistics for exposure variables.
#'   - `"pipeline"`: Overview of steps completed in the pipeline.
#'
#'   Use `"all"` to export all of the above categories.
#'
#' @param include_empty_tabs Logical. If `TRUE`, adds placeholder sheets for
#' any missing result types. Default is `FALSE`.
#'
#' @return An Excel file is written to the specified path.
#' A message is printed with the file location.
#'
#' @examples
#' # Create example data
#' mae <- make_example_data(
#'     n_samples = 20,
#'     return_mae = TRUE
#' )
#'
#' # run correlation analysis
#' mae <- mae |>
#'     run_correlation(
#'         feature_type = "exposures",
#'         exposure_cols = c("exposure_pm25", "exposure_no2", "age", "bmi")
#'     )
#'
#' # file path of the output file
#' tmp <- tempfile(fileext = ".xlsx")
#'
#' # extract the correlation results
#' extract_results_excel(
#'     exposomicset = mae,
#'     result_types = "correlation",
#'     file = tmp
#' )
#'
#' @importFrom purrr iwalk
#' @importFrom tibble enframe
#' @export
extract_results_excel <- function(
    exposomicset,
    file = "tidyexposomics_results.xlsx",
    result_types = c(
        "correlation",
        "association",
        "differential_analysis",
        "multiomics_integration",
        "network",
        "enrichment",
        "exposure_summary",
        "pipeline"
    ),
    include_empty_tabs = FALSE) {
    .check_suggested(pkg = "openxlsx")
    stopifnot("MultiAssayExperiment" %in% class(exposomicset))

    # Handle "all"
    all_types <- c(
        "correlation",
        "association",
        "differential_analysis",
        "multiomics_integration",
        "network",
        "enrichment",
        "exposure_summary",
        "pipeline"
    )
    if ("all" %in% result_types) {
        result_types <- all_types
    }

    wb <- openxlsx::createWorkbook()

    safe_add_sheet <- function(name, data) {
        if (!is.null(data) && is.data.frame(data)) {
            openxlsx::addWorksheet(wb, name)
            openxlsx::writeData(wb, name, data, withFilter = TRUE)

            # Bold headers + border around header row
            header_style <- openxlsx::createStyle(
                textDecoration = "bold",
                border = "TopBottomLeftRight",
                borderStyle = "thin"
            )

            openxlsx::addStyle(
                wb,
                sheet = name, style = header_style,
                rows = 1, cols = seq_len(ncol(data)), gridExpand = TRUE
            )

            openxlsx::setColWidths(
                wb,
                sheet = name,
                cols = seq_len(ncol(data)), widths = "auto"
            )
        } else if (include_empty_tabs) {
            openxlsx::addWorksheet(wb, name)
            openxlsx::writeData(wb, name, data.frame(Note = "No data available."))
        }
    }

    # --- Adding Correlation Results ---------
    if ("correlation" %in% result_types) {
        message("Writing Correlation Results.")
        enr <- MultiAssayExperiment::metadata(exposomicset)$correlation
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
        enr <- MultiAssayExperiment::metadata(exposomicset)$association
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
        da <- MultiAssayExperiment::metadata(exposomicset)$differential_analysis$differential_abundance
        if (!is.null(da)) {
            safe_add_sheet("Differential Abundance", da)
        }

        if ("sensitivity_analysis" %in% names(
            MultiAssayExperiment::metadata(exposomicset)$differential_analysis
        )) {
            message("Writing Sensitivity Analysis Results.")
            sens <- MultiAssayExperiment::metadata(exposomicset)$differential_analysis$sensitivity_analysis$feature_stability
            if (!is.null(sens)) {
                safe_add_sheet("Sensitivity Analysis", sens)
            }
        }
    }

    # --- Adding Multiomics Integration Results ----------
    if ("multiomics_integration" %in% result_types) {
        message("Writing Multiomics Integration Results.")
        factors <- MultiAssayExperiment::metadata(exposomicset)$multiomics_integration$common_top_factors
        if (!is.null(factors)) {
            safe_add_sheet("Common Factor Features", factors)
        }
    }

    # --- Adding Network Results --------------
    if ("network" %in% result_types) {
        message("Writing Network Impact Results.")
        exposure_impact <- MultiAssayExperiment::metadata(exposomicset) |>
            purrr::pluck(
                "network",
                "exposure_impact"
            )
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
        enr <- MultiAssayExperiment::metadata(exposomicset)$enrichment
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
        pipeline_info <- MultiAssayExperiment::metadata(exposomicset)$summary$steps
        if (is.list(pipeline_info)) {
            steps_df <- tibble::enframe(pipeline_info,
                name = "Step",
                value = "Details"
            )
            safe_add_sheet("Pipeline_Steps", steps_df)
        }
    }

    # --- Adding Exposure summaries -------
    if ("exposure_summary" %in% result_types) {
        message("Writing Exposure Summary Results.")
        safe_add_sheet(
            "Exposure_Summary",
            tryCatch(
                {
                    run_summarize_exposures(exposomicset, action = "get")
                },
                error = function(e) NULL
            )
        )
    }

    # Save
    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    message("Results written to: ", normalizePath(file))
}
