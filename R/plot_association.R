#' Plot Association Results
#'
#' Generates a forest plot for association results from any source stored
#' in the metadata of a `MultiAssayExperiment` object.
#' Supports faceting and visual augmentation with R^2 tiles when available.
#'
#' @param exposomicset A `MultiAssayExperiment` object containing association
#' results in metadata.
#' @param source Character string indicating the association source.
#' One of `"omics"`, `"exposures"`, `"factors"`, or `"go_pcs"`.
#' @param model Name of the stored model to plot (see [association_models()]).
#' Can be left as `NULL` when only one model is stored for the source.
#' @param terms Optional character vector of term names to subset the plot to.
#'  Default is `NULL` (include all).
#' @param effect Character; which effect to plot. `"main"` (default) is the
#' feature's main effect. For models with interactions, supply the name of the
#' interacting coefficient (e.g. `"sexMale"`) to plot the interaction instead.
#' @param filter_col Column used to assess statistical significance
#' (default: `"p.value"`).
#' @param filter_thresh Numeric threshold applied to `filter_col`
#' (default: `0.05`).
#' @param direction_filter Direction of associations to retain.
#' One of `"all"` (default), `"up"`, or `"down"`.
#' @param add_r2_tile Logical; if `TRUE`, includes a tile plot for `r2_col`
#'  (default: `TRUE`).
#' @param r2_col Column used for coloring the tile plot (default: `"adj_r2"`).
#' @param up_col Character; color used for significant positive associations
#'  (default: `"#8E0152"`).
#' @param down_col Character; color used for significant negative associations
#'  (default: `"#006666"`).
#' @param ns_col Character; color used for non-significant associations
#'  (default: `"grey55"`).
#' @param r2_low_col Character; color for the low end of the R^2 tile gradient
#'  (default: `"#e8e1eb"`). Only used when `add_r2_tile = TRUE`.
#' @param r2_high_col Character; color for the high end of the R^2 tile
#'  gradient (default: `"#3f1b4f"`). Only used when `add_r2_tile = TRUE`.
#' @param facet Logical; if `TRUE` and `source == "go_pcs"`,
#'  apply nested faceting by experiment and GO cluster (default: `FALSE`).
#' @param nrow Integer; number of rows for facet layout if enabled (default: `1`).
#' @param subtitle Optional subtitle for the plot. If `NULL`,
#' automatically generated from the model formula and sample subset.
#'
#' @return A `ggplot2` object: either a single forest plot or a composite plot
#'  with an R^2 tile strip.
#'
#' @details
#' This function visualizes effect size estimates and confidence intervals from
#' association analyses. It allows filtering by
#' direction (`"up"` for positive, `"down"` for negative) and by significance.
#'
#' The R^2 tile (if enabled) adds a side heatmap indicating model fit for
#' each association. This can be useful for model diagnostics.
#'
#' To see several models together, use `compare_association()` and `plot_association_comparison()`.
#'
#' @examples
#' mae <- make_example_data(
#'     n_samples = 10,
#'     return_mae = TRUE
#' )
#'
#' mae <- mae |>
#'     run_association(
#'         formula = smoker ~ age,
#'         source = "exposures",
#'         feature_set = c("exposure_pm25", "exposure_no2"),
#'         family = "binomial"
#'     )
#'
#' assoc_plot <- mae |>
#'     plot_association(
#'         source = "exposures"
#'     )
#'
#' @export
plot_association <- function(
  exposomicset,
  source = c("omics", "exposures", "factors", "go_pcs"),
  model = NULL,
  terms = NULL,
  effect = "main",
  filter_col = "p.value",
  filter_thresh = 0.05,
  direction_filter = "all",
  add_r2_tile = TRUE,
  r2_col = "adj_r2",
  up_col = "#8E0152",
  down_col = "#006666",
  ns_col = "grey55",
  r2_low_col = "#e8e1eb",
  r2_high_col = "#3f1b4f",
  facet = FALSE,
  nrow = 1,
  subtitle = NULL
) {
    requireNamespace("ggplot2")
    requireNamespace("dplyr")
    requireNamespace("tidyr")
    .check_suggested(pkg = "patchwork")

    source <- match.arg(source)

    assoc_result <- .get_assoc_model(exposomicset, source, model)
    covariates <- assoc_result$covariates

    df <- assoc_result$results_df |>
        dplyr::filter(effect == !!effect) |>
        dplyr::filter(is.null(terms) | term %in% terms)

    # Mark significance
    df <- df |>
        dplyr::mutate(
            is_significant = !!rlang::sym(filter_col) < filter_thresh,
            direction = dplyr::case_when(
                (estimate > 0 & (is_significant == TRUE)) ~ "up",
                (estimate < 0 & (is_significant == TRUE)) ~ "down",
                .default = "ns"
            )
        )

    # Apply direction filter
    if (direction_filter == "up") {
        df <- dplyr::filter(df, direction == "up")
    } else if (direction_filter == "down") {
        df <- dplyr::filter(df, direction == "down")
    }

    # Extract y-axis variable and handle GO PC formatting
    if (source == "go_pcs") {
        df <- tidyr::separate(
            df,
            col = "term",
            into = c("PC", "exp_name", "Cluster", "go_group"),
            sep = "/",
            remove = FALSE
        )
        df$label <- df$go_group
        facet_formula <- exp_name + Cluster ~ .
    } else {
        df$label <- df$term
        facet_formula <- NULL
    }

    # Only proceed if add_r2_tile is TRUE and column exists
    if (add_r2_tile && r2_col %in% names(df)) {
        leg_label <- ifelse(
            r2_col == "r2",
            expression("R"^2),
            expression("Adj. R"^2)
        )

        tile_plot <- df |>
            ggplot2::ggplot(ggplot2::aes(
                x = 1,
                y = stats::reorder(label, estimate),
                fill = !!rlang::sym(r2_col)
            )) +
            ggplot2::geom_tile() +
            ggplot2::scale_fill_gradient(
                low = r2_low_col,
                high = r2_high_col
            ) +
            ggplot2::labs(
                x = NULL,
                y = NULL,
                fill = leg_label
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                axis.text.x = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank(),
                panel.grid = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank()
            )
    }

    # Auto-subtitle
    if (is.null(subtitle)) {
        sample_key <- assoc_result$sample_key

        subtitle <- if (!is.null(assoc_result$formula)) {
            paste(
                c(
                    paste("Model:", .formula_key(assoc_result$formula)),
                    if (!is.null(sample_key) && sample_key != "all") {
                        paste("Samples:", sample_key)
                    }
                ),
                collapse = " | "
            )
        } else if (!is.null(covariates) && length(covariates) > 0) {
            paste("Covariates:", paste(covariates, collapse = ", "))
        }
    }

    # Plot
    p <- ggplot2::ggplot(
        df,
        ggplot2::aes(
            x = estimate,
            y = stats::reorder(label, estimate),
            color = direction
        )
    ) +
        ggplot2::geom_vline(
            xintercept = 0,
            linetype = "dashed"
        ) +
        ggplot2::geom_errorbar(
            ggplot2::aes(
                xmin = estimate - std.error,
                xmax = estimate + std.error
            ),
            orientation = "y",
            color = "grey55",
            width = 0.2
        ) +
        ggplot2::geom_point(
            shape = 18,
            size = 5,
            alpha = 1
        ) +
        ggplot2::scale_color_manual(
            values = c(
                up = up_col,
                down = down_col,
                ns = ns_col
            )
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            legend.position = "none",
            plot.subtitle = ggplot2::element_text(face = "italic"),
            plot.title = ggplot2::element_text(face = "bold.italic")
        ) +
        ggplot2::labs(
            x = "Effect size",
            y = "",
            title = paste0(tools::toTitleCase(source), "-Outcome Associations"),
            subtitle = subtitle
        )

    if (facet && source == "go_pcs") {
        if (!requireNamespace("ggh4x", quietly = TRUE)) {
            warning("Faceting requires the `ggh4x` package.")
        } else {
            p <- p +
                ggh4x::facet_nested_wrap(
                    facet_formula,
                    nrow = nrow,
                    scales = "free_x"
                ) +
                ggplot2::theme(
                    strip.text = ggplot2::element_text(face = "bold.italic")
                )
        }
    }

    if (add_r2_tile && exists("tile_plot")) {
        return(patchwork::wrap_plots(
            p + ggplot2::theme(axis.text.y = ggplot2::element_text()),
            tile_plot,
            widths = c(3, 0.5)
        ))
    }

    p
}


#' Plot an Association Model Comparison
#'
#' Visualise the output of `compare_association()` either as a dodged forest
#' plot (one colour per model) or as a scatter of each model's estimates
#' against the reference model's.
#'
#' @param comparison Output of `compare_association()`.
#' @param type `"forest"` (default) or `"scatter"`. The scatter plots every
#' non-reference model against the reference, with the identity line; points on
#' the line are unchanged between models.
#' @param colors character vector of colors, used to color the points.
#' @param title Plot title.
#'
#' @return A `ggplot2` object.
#'
#' @examples
#' mae <- make_example_data(
#'     n_samples = 10,
#'     return_mae = TRUE
#' ) |>
#'     run_association(
#'         formula = list(
#'             crude = smoker ~ 1,
#'             adjusted = smoker ~ age
#'         ),
#'         source = "exposures",
#'         feature_set = c("exposure_pm25", "exposure_no2"),
#'         family = "binomial"
#'     )
#'
#' mae |>
#'     compare_association(source = "exposures") |>
#'     plot_association_comparison()
#' @export
plot_association_comparison <- function(
  comparison,
  type = c("forest", "scatter"),
  colors = "default",
  title = "Association Comparison"
) {
    requireNamespace("ggplot2")
    type <- match.arg(type)

    ref <- unique(as.character(comparison$reference))

    if (type == "forest") {
        if (colors == "default") {
            pal <- c("#1E0C1F", "#00c9c1", "#cf68d9", "#6d7ea9", "#e8d5b5", "#944470", "#67abff", "#b5a7b6", "#002f4f", "#f2b459", "#4f4350", "#51aabb", "#994a4f", "#00a3cd", "#240d00", "#48817d")
        } else {
            pal <- colors
        }

        return(
            comparison |>
                dplyr::filter(!is.na(estimate)) |>
                ggplot2::ggplot(ggplot2::aes(
                    x = estimate,
                    y = stats::reorder(term, ref_estimate, FUN = mean, na.rm = TRUE),
                    color = model,
                    group = model
                )) +
                ggplot2::geom_errorbar(
                    ggplot2::aes(
                        xmin = estimate - std.error,
                        xmax = estimate + std.error
                    ),
                    orientation = "y",
                    color = "grey55",
                    width = 0.2,
                    position = ggplot2::position_dodge(width = 0.7)
                ) +
                ggplot2::geom_point(
                    shape = 18,
                    size = 5,
                    alpha = 1,
                    position = ggplot2::position_dodge(width = 0.7)
                ) +
                ggplot2::geom_vline(
                    xintercept = 0,
                    linetype = "dashed"
                ) +
                ggplot2::theme_bw() +
                ggplot2::theme(
                    plot.title = ggplot2::element_text(face = "bold.italic")
                ) +
                ggplot2::scale_color_manual(values = pal) +
                ggplot2::labs(
                    x = "Effect size",
                    y = "",
                    color = "Model",
                    title = title
                )
        )
    }

    comparison |>
        dplyr::filter(
            as.character(model) != reference,
            !is.na(estimate),
            !is.na(ref_estimate)
        ) |>
        ggplot2::ggplot(ggplot2::aes(
            x = ref_estimate,
            y = estimate
        )) +
        ggplot2::geom_abline(
            slope = 1,
            intercept = 0,
            linetype = "dashed"
        ) +
        ggplot2::geom_hline(
            yintercept = 0,
            color = "grey70"
        ) +
        ggplot2::geom_vline(
            xintercept = 0,
            color = "grey70"
        ) +
        ggplot2::geom_point(
            shape = 18,
            size = 4
        ) +
        ggplot2::facet_wrap(~model) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold.italic"),
            strip.text = ggplot2::element_text(face = "bold.italic")
        ) +
        ggplot2::labs(
            x = paste("Effect size:", ref[[1]]),
            y = "Effect size: comparison model",
            title = title
        )
}
