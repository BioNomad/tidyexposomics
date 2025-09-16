#' Plot Correlation Tilemap
#'
#' Visualizes a correlation matrix as a heatmap tile plot using
#' correlation results stored in the metadata of a
#' `MultiAssayExperiment` object. When `feature_type = "pcs"`,
#' the function forces PCs to appear on the x-axis and exposures on the y-axis,
#' and it adds a barplot showing how many PCs are significantly
#' associated with each exposure. It also suppresses nonsignificant tiles
#'  based on a specified p-value threshold.
#'
#' @param exposomicset A `MultiAssayExperiment` object containing correlation
#' results in metadata.
#' @param feature_type Type of correlation results to plot.
#'  One of `"pcs"`, `"degs"`, `"omics"`,
#' `"factors"`, `"factor_features"`, or `"exposures"`.
#'  Must match the key used in
#' `metadata(exposomicset)$correlation[[feature_type]]`.
#' @param pval_cutoff Numeric p-value cutoff below which
#' correlations are displayed.
#' Nonsignificant tiles are rendered in the `na_color`.
#' Default is `0.05`.
#' @param na_color Color used to represent nonsignificant or
#'  missing correlations. Default is `"grey100"`.
#' @param fill_limits Numeric vector of length 2 defining the scale
#' limits for correlation values. Default is `c(-1, 1)`.
#' @param midpoint Numeric value for centering the fill gradient.
#' Default is `0`.
#'
#' @return A `ggplot2` tile plot (or a combined tile + barplot if
#'  `feature_type = "pcs"`).
#'
#' @examples
#'
#' # create example data
#' mae <- make_example_data(
#'     n_samples = 10,
#'     return_mae = TRUE
#' )
#'
#' # run pca
#' mae <- mae |>
#'     run_pca()
#'
#' # correlate with exposures
#' mae <- mae |>
#'     run_correlation(
#'         feature_type = "pcs",
#'         exposure_cols = c("exposure_pm25", "exposure_no2", "age", "bmi")
#'     )
#'
#' # make the correlation tile plot
#' cor_tile_p <- mae |>
#'     plot_correlation_tile(
#'         feature_type = "pcs"
#'     )
#'
#' @export
plot_correlation_tile <- function(
    exposomicset,
    feature_type = c(
        "pcs",
        "degs",
        "omics",
        "factors",
        "factor_features",
        "exposures"
    ),
    pval_cutoff = 0.05,
    na_color = "grey100",
    fill_limits = c(-1, 1),
    midpoint = 0) {
    .check_suggested(pkg = "patchwork")
    feature_type <- match.arg(feature_type)

    correlation_list <- MultiAssayExperiment::metadata(exposomicset)$correlation
    if (is.null(correlation_list) || !feature_type %in% names(correlation_list)) {
        stop(
            "Correlation results for feature_type '",
            feature_type,
            "' not found in metadata."
        )
    }

    correlation_df <- correlation_list[[feature_type]]

    if (!all(c("var1", "var2", "correlation", "p.value") %in%
        colnames(correlation_df))) {
        cols_needed <- "var1, var2, correlation, p.value."
        stop("Correlation data must contain columns: ", cols_needed)
    }

    extract_pc_number <- function(x) {
        as.numeric(stringr::str_extract(x, "\\d+"))
    }

    if (feature_type == "pcs") {
        # Identify PCs in each var
        is_pc1 <- grepl("^PC\\d+$", correlation_df$var1)
        is_pc2 <- grepl("^PC\\d+$", correlation_df$var2)

        # Keep only PC vs non-PC pairs
        keep <- xor(is_pc1, is_pc2)
        df <- correlation_df[keep, , drop = FALSE]

        # Assign PC to x-axis, exposure to y-axis
        df <- df |>
            dplyr::mutate(
                var_pc = dplyr::case_when(is_pc1 ~ var1,
                    is_pc2 ~ var2,
                    .default = NA_character_
                ),
                var_exp = dplyr::case_when(is_pc1 ~ var2,
                    is_pc2 ~ var1,
                    .default = NA_character_
                ),
                correlation = dplyr::case_when(
                    p.value < pval_cutoff ~ correlation,
                    .default = NA_real_
                )
            ) |>
            dplyr::filter(!is.na(var_pc), !is.na(var_exp)) |>
            dplyr::mutate(
                var_pc = factor(var_pc, levels = unique(var_pc)[
                    order(extract_pc_number(unique(var_pc)))
                ])
            )

        # Count significant PC associations per exposure
        bar_df <- df |>
            dplyr::filter(!is.na(correlation)) |>
            dplyr::count(var_exp, name = "n")

        # Order exposures by number of PC associations
        bar_df <- bar_df |>
            dplyr::arrange(n) |>
            dplyr::mutate(var_exp = factor(var_exp, levels = var_exp))

        # Reorder var_exp in df accordingly
        df <- df |>
            dplyr::filter(var_exp %in% bar_df$var_exp) |>
            dplyr::mutate(var_exp = factor(var_exp,
                levels = levels(bar_df$var_exp)
            ))

        # Create heatmap
        heatmap <- ggplot2::ggplot(df, ggplot2::aes(
            x = var_pc,
            y = var_exp,
            fill = correlation
        )) +
            ggplot2::geom_tile() +
            ggplot2::scale_fill_gradient2(
                low = "blue4", high = "red4", mid = "white",
                midpoint = midpoint, limits = fill_limits,
                na.value = na_color, name = "Correlation"
            ) +
            ggplot2::theme_bw() +
            ggplot2::labs(
                x = "",
                y = ""
            ) +
            ggpubr::rotate_x_text(90)

        # Create barplot
        bar_plot <- ggplot2::ggplot(
            bar_df,
            ggplot2::aes(
                x = n,
                y = var_exp
            )
        ) +
            ggplot2::geom_col(fill = "grey40") +
            ggplot2::theme_classic() +
            ggplot2::labs(
                x = "No. of PCs",
                y = NULL
            ) +
            ggplot2::theme(
                axis.text.y = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_blank()
            )

        # Combine plots
        combined_plot <- (heatmap + bar_plot) +
            patchwork::plot_layout(
                widths = c(4, 1),
                guides = "collect"
            ) &
            ggplot2::theme(legend.position = "right")

        return(combined_plot)
    } else {
        # Handle symmetric correlation matrices
        df <- correlation_df |>
            dplyr::filter(var1 != var2) |>
            dplyr::mutate(
                var_min = pmin(var1, var2),
                var_max = pmax(var1, var2)
            ) |>
            dplyr::select(-var1, -var2) |>
            dplyr::rename(var1 = var_min, var2 = var_max) |>
            dplyr::mutate(
                correlation = dplyr::case_when(
                    p.value < pval_cutoff ~ correlation,
                    .default = NA_real_
                )
            )

        p <- ggplot2::ggplot(df, ggplot2::aes(
            x = var2,
            y = var1,
            fill = correlation
        )) +
            ggplot2::geom_tile() +
            ggplot2::scale_fill_gradient2(
                low = "blue", high = "red", mid = "white",
                midpoint = midpoint, limits = fill_limits,
                na.value = na_color, name = "Correlation"
            ) +
            ggplot2::theme_bw() +
            ggplot2::labs(x = NULL, y = NULL) +
            ggpubr::rotate_x_text(angle = 90)

        return(p)
    }
}
