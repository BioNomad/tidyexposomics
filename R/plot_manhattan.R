#' Plot a Manhattan-style ExWAS summary across omics categories
#'
#' This function generates a multi-faceted Manhattan plot from the results of
#' `run_association()`, visualizing the significance of associations
#' across omics features and exposures, grouped by category. Significant
#' features can be highlighted and labeled, and strip backgrounds can be
#' colored per facet.
#'
#' @param exposomicset A `MultiAssayExperiment` object that has
#' already been processed by `run_association()` for the `"omics"` source,
#' the `"exposures"` source, or both.
#' @param omics_model,exposures_model Names of the stored models to plot for
#' the omics and exposures sources (see [association_models()]). Can be left
#' as `NULL` when only one model is stored for a source.
#' @param effect Which effect to plot: `"main"` (default) is the feature's main
#' effect. For models with interactions, supply an interacting coefficient
#' name (e.g. `"sexmale"`).
#' @param pval_thresh Numeric threshold for significance (default = 0.05).
#' @param feature_col A character string indicating the column name to use
#' for feature labeling and highlighting (e.g., `"term"` or `"feature"`).
#' Default is `"term"`.
#' @param alpha Transparency applied to facet strip colors (default = 0.5).
#' @param min_per_cat Minimum number of features per category to be shown
#' (default = 1).
#' @param vars_to_label Optional character vector of variable names to label
#' explicitly, matched against the `feature_col` column (or `term`).
#' @param sig_color Color used for significant points (default = `"magenta2"`).
#' @param non_sig_cols Character vector of alternating colors for
#' non-significant points (default = `c("grey25", "grey75")`).
#' @param pval_thresh_line_col Color of the horizontal significance
#' threshold line (default = `"grey25"`).
#' @param panel_sizes Numeric vector passed to `ggh4x::force_panelsizes()`
#' to control panel widths (default = `c(1,1,1,1,1)`).
#' @param linetype Line type for the horizontal threshold (default = `"dashed"`).
#' @param facet_cols Optional vector of colors to use for
#' facet strip backgrounds.
#' @param label_size Numeric size of the feature label text (default = 3.5).
#' @param facet_angle Angle (in degrees) for strip text rotation (default = 90).
#' @param facet_text_face Font face for facet strip labels
#' (default = `"bold.italic"`).
#'
#' @return A `ggplot` object showing the Manhattan-style faceted plot.
#'
#' @details
#' - This function expects `run_association()` to have been run first, for
#'   `source = "omics"`, `source = "exposures"`, or both. Whichever is present
#'   is plotted.
#' - Facets represent categories (omics experiments, or the `category` column
#'   of the exposure codebook), and points represent features.
#' - Points below the significance threshold are colored using `non_sig_cols`,
#'   while significant ones are colored with `sig_color` and optionally labeled.
#' - Uses `ggrepel` to avoid overlapping labels and `ggh4x` for
#' enhanced faceting.
#' - The `feature_col` argument allows customization of which column is used
#'   to label or identify features, enabling compatibility with different
#'   result formats.
#'
#' @examples
#' # create example data
#' mae <- make_example_data(
#'     n_samples = 10,
#'     return_mae = TRUE
#' )
#'
#' # run association tests
#' mae <- mae |>
#'     run_association(
#'         source = "omics",
#'         formula = smoker ~ age,
#'         top_n = 20,
#'         family = "binomial"
#'     )
#'
#' # create the manhattan plot
#' manhattan_p <- mae |>
#'     plot_manhattan(feature_col = "term")
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_vline scale_color_identity labs theme element_text element_rect
#' @importFrom dplyr filter mutate group_by ungroup n bind_rows pull
#' @importFrom purrr pluck
#' @importFrom ggrepel geom_label_repel
#' @importFrom scales alpha
#' @importFrom ggpubr theme_pubr
#'
#' @export
plot_manhattan <- function(
  exposomicset,
  omics_model = NULL,
  exposures_model = NULL,
  effect = "main",
  pval_thresh = 0.05,
  feature_col = "term",
  alpha = 0.5,
  min_per_cat = 1,
  vars_to_label = NULL,
  sig_color = "magenta2",
  non_sig_cols = c("grey25", "grey75"),
  pval_thresh_line_col = "grey25",
  panel_sizes = c(1, 1, 1, 1, 1),
  linetype = "dashed",
  facet_cols = NULL,
  label_size = 3.5,
  facet_angle = 90,
  facet_text_face = "bold.italic"
) {
    .check_suggested(pkg = "forcats")
    .check_suggested(pkg = "ggh4x")

    if (!("association" %in% names(MultiAssayExperiment::metadata(exposomicset)))) {
        stop("Please run `run_association()` first.")
    }

    # Pull whichever sources have results (the omics or exposures source can
    # be missing), keep the requested effect, and keep only the columns used
    # here so differing annotation columns can't clash when binding.
    keep_cols <- unique(c("term", feature_col, "category", "p.value"))

    collect_source <- function(source, model) {
        res <- .get_assoc_results(
            exposomicset,
            source,
            model,
            arg = paste0(source, "_model")
        )

        if (is.null(res)) {
            return(NULL)
        }

        res <- res |>
            dplyr::filter(effect == !!effect)

        # exposures take their label from `term`; omics may use a rowData column
        if (source == "exposures" || !feature_col %in% names(res)) {
            res[[feature_col]] <- res$term
        }

        # exposure categories come from the codebook
        if (!"category" %in% names(res)) {
            res$category <- if (source == "exposures") "Exposures" else NA_character_
        }

        res |>
            dplyr::select(dplyr::any_of(keep_cols)) |>
            dplyr::mutate(category = as.character(category))
    }

    results <- dplyr::bind_rows(
        collect_source("omics", omics_model),
        collect_source("exposures", exposures_model)
    )

    if (nrow(results) == 0) {
        stop(
            "No association results found for sources 'omics' or 'exposures'. ",
            "Run `run_association()` first."
        )
    }

    manhattan_data <- results |>
        dplyr::filter(!is.na(category)) |>
        dplyr::mutate(category = gsub("_", " ", category)) |>
        dplyr::mutate(
            var = forcats::fct_reorder(as.character(term), category),
            thresh_met = ifelse(p.value < pval_thresh, "yes", "no")
        )

    if (nrow(manhattan_data) == 0) {
        stop(
            "No results have a `category`. Omics results get one from the ",
            "experiment name, and exposures from the `category` column of the ",
            "codebook in metadata."
        )
    }

    if (!is.null(vars_to_label)) {
        manhattan_data <- manhattan_data |>
            dplyr::mutate(label = ifelse(
                as.character(var) %in% vars_to_label |
                    as.character(!!dplyr::sym(feature_col)) %in% vars_to_label,
                as.character(!!dplyr::sym(feature_col)), NA
            )) |>
            dplyr::mutate(
                label = dplyr::case_when(
                    !is.na(label) ~ paste0("italic(", label, ")"),
                    .default = NA_character_
                )
            )
    } else {
        manhattan_data <- manhattan_data |>
            dplyr::mutate(
                label = ifelse(thresh_met == "yes",
                    as.character(!!dplyr::sym(feature_col)),
                    NA
                )
            )
    }

    # Filter out categories with fewer than min_per_cat features
    manhattan_data <- manhattan_data |>
        dplyr::group_by(category) |>
        dplyr::filter(dplyr::n() >= min_per_cat) |>
        dplyr::ungroup()

    # Get unique ordered categories
    ordered_cats <- manhattan_data$category |>
        unique() |>
        sort()

    # Alternating point colors for non-significant values
    category_grey_map <- stats::setNames(
        rep(non_sig_cols, length.out = length(ordered_cats)),
        ordered_cats
    )

    # Add color column to data
    manhattan_data <- manhattan_data |>
        dplyr::mutate(
            cols = ifelse(
                thresh_met == "yes",
                sig_color,
                category_grey_map[category]
            )
        )

    # Facet colors
    if (is.null(facet_cols)) {
        facet_cols <- tidy_exp_pal[seq_len(dplyr::n_distinct(manhattan_data$category))]
    }

    manhattan_data |>
        ggplot(aes(
            y = var,
            x = -log10(p.value),
            color = cols
        )) +
        geom_point() +
        geom_vline(
            xintercept = -log10(pval_thresh),
            linetype = linetype,
            color = pval_thresh_line_col
        ) +
        ggrepel::geom_label_repel(
            data = manhattan_data,
            aes(
                label = label,
                fontface = "bold",
                fill = "white"
            ),
            parse = TRUE,
            size = label_size,
            fill = "white",
            max.overlaps = 5,
            min.segment.length = .5,
            box.padding = 0.6,
            point.padding = 0.2,
            color = "black",
            na.rm = TRUE
        ) +
        ggh4x::facet_grid2(
            category ~ .,
            scales = "free_y",
            space = "free_y",
            strip = ggh4x::strip_themed(
                background_y = ggh4x::elem_list_rect(
                    fill = scales::alpha(
                        facet_cols,
                        alpha
                    )
                )
            )
        ) +
        ggh4x::force_panelsizes(rows = panel_sizes) +
        scale_color_identity() +
        ggpubr::theme_pubr(legend = "none") +
        theme(
            axis.ticks.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            strip.placement = "outside",
            strip.background = element_rect(linewidth = 0, colour = NA),
            strip.text.y = element_text(
                angle = facet_angle,
                face = facet_text_face
            ),
            panel.spacing = ggplot2::unit(0, "lines")
        ) +
        labs(
            y = "",
            x = expression(-Log[10] * "P")
        )
}
