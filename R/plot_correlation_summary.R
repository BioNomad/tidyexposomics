#' Plot Correlation Summary from Exposure-Feature Correlations
#'
#' Generates a bar plot summary of exposure-feature correlations using customizable modes.
#'
#' @param expomicset A `MultiAssayExperiment` object with correlation results in metadata.
#' @param feature_type One of `"degs"`, `"factors"`, `"omics"`, or `"exposures"`.
#' @param mode One of:
#'   - `"top_exposures"`: Top exposures by assay (like `plot_bar_correlate_exposure_fill`)
#'   - `"top_features"`: Top features by exposure category (like `plot_bar_correlate_feature_fill`)
#'   - `"exposure_category"`: Total associations by exposure category
#'   - `"assay"`: Total associations by omics assay
#'   - `"summary"`: Patchwork layout combining all
#' @param top_n Number of top exposures or features to display (for top modes). Default is `15`.
#'
#' @return A `ggplot2` object or a `patchwork` object (if `mode = "summary"`).
#' @export
plot_correlation_summary <- function(
    expomicset,
    feature_type = c("degs", "factors"),
    mode = c("top_exposures", "top_features", "exposure_category", "assay", "summary"),
    top_n = 15
) {
  requireNamespace("ggplot2")
  requireNamespace("patchwork")
  requireNamespace("janitor")

  feature_type <- match.arg(feature_type)
  mode <- match.arg(mode)

  # Select correlation data
  cor_data <- MultiAssayExperiment::metadata(expomicset)$correlation[[feature_type]]
  if (is.null(cor_data)) stop("No correlation results found for type: ", feature_type)

  # # Grab experiment names
  # exp_names <- names(MultiAssayExperiment::experiments(expomicset))
  #
  # # Change names of the columns to match column type
  # cor_data <- cor_data |>
  #   dplyr::rename(
  #     exposure = var1,
  #     feature = var2
  #   ) |>
  #   # Add in exp_name for the assay
  #   dplyr::mutate(
  #     exp_name = stringr::str_extract(
  #       feature,
  #       paste0("^(", paste0(stringr::str_replace_all(exp_names, " ", "[ _]"), collapse = "|"), ")")
  #     ),
  #     # Remove matched exp_name and any underscore or space after it
  #     feature = dplyr::case_when(
  #       !is.na(exp_name) ~ stringr::str_remove(feature, paste0("^", exp_name, "[ _]?")),
  #       TRUE ~ feature
  #     )
  #   )
  #
  #
  # # Merge with exposure metadata
  # cor_data <- cor_data |>
  #   dplyr::left_join(
  #     MultiAssayExperiment::metadata(expomicset)$var_info,
  #     by = c("exposure"="variable")
  #   )

  # Helper plots
  plot_top_exposures <- function() {
    cor_data |>
      dplyr::group_by(exposure, exp_name) |>
      dplyr::reframe(n_exp_name = dplyr::n()) |>
      dplyr::group_by(exposure) |>
      dplyr::mutate(total = sum(n_exp_name)) |>
      dplyr::ungroup() |>
      dplyr::filter(exposure %in% head(
        cor_data |> janitor::tabyl(exposure) |> dplyr::arrange(desc(n)) |> dplyr::pull(exposure),
        top_n)) |>
      ggplot2::ggplot(ggplot2::aes(
        x = n_exp_name,
        y = forcats::fct_reorder(exposure, total),
        fill = exp_name)) +
      ggplot2::geom_bar(stat = "identity", alpha = 0.7) +
      scale_fill_tidy_exp() +
      ggpubr::theme_pubr() +
      ggplot2::labs(
        title = "Top Exposures by Assay",
        x = "No. of Associations", y = "", fill = "Assay")+
      theme(plot.title = element_text(face = "bold.italic"))
  }

  plot_top_features <- function() {
    cor_data |>
      dplyr::group_by(feature, category) |>
      dplyr::reframe(n_category = dplyr::n()) |>
      dplyr::group_by(feature) |>
      dplyr::mutate(total = sum(n_category)) |>
      dplyr::ungroup() |>
      dplyr::filter(feature %in% head(
        cor_data |> janitor::tabyl(feature) |> dplyr::arrange(desc(n)) |> dplyr::pull(feature),
        top_n)) |>
      ggplot2::ggplot(ggplot2::aes(
        x = n_category,
        y = forcats::fct_reorder(feature, total),
        fill = category)) +
      ggplot2::geom_bar(stat = "identity", alpha = 0.7) +
      scale_fill_tidy_exp() +
      ggpubr::theme_pubr() +
      ggplot2::labs(
        title = "Top Features by Exposure Category",
        x = "No. of Associations", y = "", fill = "Category")+
      theme(plot.title = element_text(face = "bold.italic"))
  }

  plot_exposure_categories <- function() {
    cor_data |>
      janitor::tabyl(category) |>
      ggplot2::ggplot(ggplot2::aes(x = n, y = reorder(category, n), fill = category)) +
      ggplot2::geom_bar(stat = "identity", alpha = 0.7) +
      ggplot2::geom_segment(ggplot2::aes(
        xend = n, x = n,
        yend = as.numeric(reorder(category, n)) + 0.45,
        y = as.numeric(reorder(category, n)) - 0.45,
        color = category), size = 1) +
      scale_fill_tidy_exp() +
      scale_color_tidy_exp()+
      ggpubr::theme_pubr() +
      ggplot2::labs(title = "Exposure Associations by Category", x = "Count", y = "")+
      theme(plot.title = element_text(face = "bold.italic"))
  }

  plot_assays <- function() {
    cor_data |>
      janitor::tabyl(exp_name) |>
      ggplot2::ggplot(ggplot2::aes(x = n, y = reorder(exp_name, n), fill = exp_name)) +
      ggplot2::geom_bar(stat = "identity", alpha = 0.7) +
      ggplot2::geom_segment(ggplot2::aes(
        xend = n, x = n,
        yend = as.numeric(reorder(exp_name, n)) + 0.45,
        y = as.numeric(reorder(exp_name, n)) - 0.45,
        color = exp_name), size = 1) +
      scale_fill_tidy_exp() +
      scale_color_tidy_exp() +
      ggpubr::theme_pubr() +
      ggplot2::labs(title = "Omics Associations by Assay", x = "Count", y = "")+
      theme(plot.title = element_text(face = "bold.italic"))
  }

  # Return appropriate plot
  switch(mode,
         top_exposures = plot_top_exposures(),
         top_features = plot_top_features(),
         exposure_category = plot_exposure_categories(),
         assay = plot_assays(),
         summary = {
           (plot_top_features() / plot_exposure_categories() + patchwork::plot_layout(heights = c(3, 1))) |
             (plot_top_exposures() / plot_assays() + patchwork::plot_layout(heights = c(3, 1)))
         })
}
