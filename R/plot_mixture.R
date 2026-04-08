#' Plot Mixture Analysis Results
#'
#' Create visualizations for mixture analysis results from qgcomp, WQS, or BKMR.
#'
#' @param exposomicset A `MultiAssayExperiment` object with mixture analysis results,
#'   or a mixture analysis result object directly.
#' @param method Which method's results to plot: `"qgcomp"`, `"wqs"`, or `"bkmr"`.
#' @param plot_type Type of plot to generate. Options depend on method:
#'
#'   For qgcomp: `"weights"` (default), `"effect"`
#'
#'   For WQS: `"weights"` (default), `"effect"`
#'
#'   For BKMR: `"pips"` (default), `"univariate"`, `"overall"`, `"interaction"`
#' @param threshold For weight plots, show reference line at 1/n_exposures.
#'   Default is `TRUE`.
#' @param top_n For weight plots, only show top n exposures by weight.
#'   Default is `NULL` (show all).
#' @param ... Additional arguments passed to underlying plot functions.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' mae <- make_example_data(n_samples = 100, return_mae = TRUE)
#'
#' mae <- mae |>
#'     run_mixture_analysis(
#'         outcome = "outcome_bmi",
#'         exposures = c("exposure_pm25", "exposure_no2"),
#'         method = "qgcomp",
#'         n_boot = 0
#'     )
#'
#' plot_mixture(mae, method = "qgcomp", plot_type = "weights")
#'
#' @export
plot_mixture <- function(
  exposomicset,
  method = c("qgcomp", "wqs", "bkmr"),
  plot_type = NULL,
  threshold = TRUE,
  top_n = NULL,
  ...
) {
    method <- match.arg(method)

    # Extract results from MAE or use directly
    if (inherits(exposomicset, "MultiAssayExperiment")) {
        results <- MultiAssayExperiment::metadata(exposomicset)$mixture_analysis[[method]]
        if (is.null(results)) {
            stop(sprintf("No %s results found. Run run_mixture_analysis() first.", method))
        }
    } else if (is.list(exposomicset)) {
        results <- exposomicset
        # Check if method matches
        if (!is.null(results$method) && results$method != method) {
            stop(sprintf(
                "Results are from %s but method = '%s' was specified.",
                toupper(results$method), method
            ))
        }
    } else {
        stop("Input must be a MultiAssayExperiment or mixture analysis result list.")
    }

    # Set default plot types
    if (is.null(plot_type)) {
        plot_type <- switch(method,
            qgcomp = "weights",
            wqs = "weights",
            bkmr = "pips"
        )
    }

    # Set appropriate plotting function
    switch(method,
        qgcomp = .plot_qgcomp(results, plot_type, threshold, top_n, ...),
        wqs = .plot_wqs(results, plot_type, threshold, top_n, ...),
        bkmr = .plot_bkmr(results, plot_type, ...)
    )
}

# ---- qgcomp plots ----

#' Plot qgcomp results
#' @keywords internal
#' @noRd
.plot_qgcomp <- function(
  results,
  plot_type,
  threshold,
  top_n,
  ...
) {
    switch(plot_type,
        weights = .plot_qgcomp_weights(results, threshold, top_n),
        effect = .plot_qgcomp_effect(results),
        stop("plot_type must be 'weights' or 'effect'")
    )
}

#' Plot qgcomp weights
#' @keywords internal
#' @noRd
.plot_qgcomp_weights <- function(
  results,
  threshold,
  top_n
) {
    weights <- results$weights
    n_exp <- nrow(weights)

    # Reshape for plotting, where positive weights go right and negative go left
    plot_data <- weights |>
        dplyr::mutate(
            weight = dplyr::case_when(
                direction == "positive" ~ weight_pos,
                direction == "negative" ~ -weight_neg,
                TRUE ~ 0
            )
        ) |>
        dplyr::arrange(weight) |>
        dplyr::mutate(exposure = factor(exposure, levels = exposure))

    # Filter to top_n if specified
    if (!is.null(top_n) && top_n < n_exp) {
        plot_data <- plot_data |>
            dplyr::slice(c(
                1:ceiling(top_n / 2),
                (dplyr::n() - floor(top_n / 2) + 1):dplyr::n()
            ))
    }

    # Extract mixture effect for subtitle
    mix_est <- results$mixture_effect$estimate[1]
    mix_pval <- results$mixture_effect$p_value[1]

    p <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(
            x = weight,
            y = exposure,
            fill = direction
        )
    ) +
        ggplot2::geom_col() +
        ggplot2::scale_fill_manual(
            values = c(
                "negative" = "blue4",
                "positive" = "red4",
                "none" = "grey50"
            ),
            name = "Direction"
        ) +
        ggplot2::geom_vline(
            xintercept = 0,
            linetype = "solid",
            color = "black"
        ) +
        ggplot2::labs(
            title = "Quantile G-Computation Weights",
            subtitle = sprintf(
                "Mixture effect (psi): %.3f (p = %.3f)",
                mix_est,
                mix_pval
            ),
            x = "Weight",
            y = "Exposure"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid.major.y = ggplot2::element_blank(),
            legend.position = "bottom"
        )

    # Add threshold lines
    if (threshold && n_exp > 1) {
        thresh_val <- 1 / n_exp
        p <- p +
            ggplot2::geom_vline(
                xintercept = c(-thresh_val, thresh_val),
                linetype = "dashed", color = "grey40", alpha = 0.7
            )
    }

    p
}

#' Plot qgcomp mixture effect
#' @keywords internal
#' @noRd
.plot_qgcomp_effect <- function(
  results
) {
    effect <- results$mixture_effect
    partial <- results$partial_effects

    effect_data <- dplyr::bind_rows(
        tibble::tibble(
            term = "Overall",
            estimate = effect$estimate[1],
            lower = effect$lower_ci[1],
            upper = effect$upper_ci[1],
            type = "overall"
        ),
        tibble::tibble(
            term = "Positive",
            estimate = partial$partial_effect[1],
            lower = NA_real_,
            upper = NA_real_,
            type = "partial"
        ),
        tibble::tibble(
            term = "Negative",
            estimate = partial$partial_effect[2],
            lower = NA_real_,
            upper = NA_real_,
            type = "partial"
        )
    ) |>
        dplyr::mutate(term = factor(
            term,
            levels = c("Negative", "Overall", "Positive")
        ))

    ggplot2::ggplot(
        effect_data,
        ggplot2::aes(
            x = term,
            y = estimate,
            fill = type
        )
    ) +
        ggplot2::geom_col(width = 0.6) +
        ggplot2::geom_errorbar(
            ggplot2::aes(
                ymin = lower,
                ymax = upper
            ),
            width = 0.2,
            na.rm = TRUE
        ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::scale_fill_manual(
            values = c(
                "overall" = "#3A3A52",
                "partial" = "#947085"
            ),
            guide = "none"
        ) +
        ggplot2::labs(
            title = "Mixture Effects",
            subtitle = sprintf(
                "Overall psi = %.3f [%.3f, %.3f]",
                effect$estimate[1],
                effect$lower_ci[1],
                effect$upper_ci[1]
            ),
            x = NULL,
            y = "Effect Estimate"
        ) +
        ggplot2::theme_bw() +
        ggplot2::coord_flip()
}

# ---- WQS plots ----

#' Plot WQS results
#' @keywords internal
#' @noRd
.plot_wqs <- function(
  results,
  plot_type,
  threshold,
  top_n,
  ...
) {
    switch(plot_type,
        weights = .plot_wqs_weights(
            results,
            threshold, top_n
        ),
        effect = .plot_wqs_effect(results),
        stop("plot_type must be 'weights' or 'effect'")
    )
}

#' Plot WQS weights
#' @keywords internal
#' @noRd
.plot_wqs_weights <- function(results, threshold, top_n) {
    weights <- results$weights
    n_exp <- nrow(weights)

    # Check if we have both directions
    has_both <- "weight_pos" %in% names(weights)

    if (has_both) {
        # Reshape for faceted plot
        plot_data <- weights |>
            tidyr::pivot_longer(
                cols = c(weight_pos, weight_neg),
                names_to = "direction",
                values_to = "weight"
            ) |>
            dplyr::mutate(
                direction = dplyr::case_when(
                    direction == "weight_pos" ~ "Positive",
                    direction == "weight_neg" ~ "Negative"
                )
            ) |>
            dplyr::group_by(direction) |>
            dplyr::arrange(dplyr::desc(weight)) |>
            dplyr::ungroup()

        # Filter to top_n per direction
        if (!is.null(top_n) && top_n < n_exp) {
            plot_data <- plot_data |>
                dplyr::group_by(direction) |>
                dplyr::slice_max(weight, n = top_n) |>
                dplyr::ungroup()
        }

        plot_data <- plot_data |>
            dplyr::mutate(exposure = forcats::fct_reorder(exposure, weight))

        p <- ggplot2::ggplot(
            plot_data,
            ggplot2::aes(
                x = weight,
                y = exposure
            )
        ) +
            ggplot2::geom_col(fill = "blue4") +
            ggplot2::facet_wrap(~direction,
                scales = "free_y"
            ) +
            ggplot2::labs(
                title = "WQS Weights by Direction",
                x = "Weight",
                y = "Exposure"
            )
    } else {
        # Single direction
        plot_data <- weights |>
            dplyr::arrange(dplyr::desc(weight)) |>
            dplyr::mutate(exposure = forcats::fct_reorder(exposure, weight))

        if (!is.null(top_n) && top_n < n_exp) {
            plot_data <- plot_data |>
                dplyr::slice_max(weight, n = top_n)
        }

        p <- ggplot2::ggplot(
            plot_data,
            ggplot2::aes(
                x = weight,
                y = exposure
            )
        ) +
            ggplot2::geom_col(fill = "blue4") +
            ggplot2::labs(
                title = "WQS Weights",
                x = "Weight",
                y = "Exposure"
            )
    }

    p <- p +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

    # Add threshold line
    if (threshold && n_exp > 1) {
        p <- p +
            ggplot2::geom_vline(
                xintercept = 1 / n_exp,
                linetype = "dashed",
                color = "magenta4"
            )
    }

    p
}

#' Plot WQS effect
#' @keywords internal
#' @noRd
.plot_wqs_effect <- function(
  results
) {
    coefs <- results$mixture_effect

    if (is.null(coefs) || nrow(coefs) == 0) {
        stop("No WQS coefficients found in results.")
    }

    plot_data <- coefs |>
        dplyr::mutate(
            term_label = dplyr::case_when(
                term == "wqs" ~ "WQS Index",
                term == "pwqs" ~ "Positive Index",
                term == "nwqs" ~ "Negative Index",
                TRUE ~ term
            ),
            significant = p_value < 0.05
        )

    ggplot2::ggplot(
        plot_data,
        ggplot2::aes(
            x = term_label,
            y = estimate,
            fill = significant
        )
    ) +
        ggplot2::geom_col(width = 0.6) +
        ggplot2::geom_errorbar(
            ggplot2::aes(
                ymin = estimate - 1.96 * std_error,
                ymax = estimate + 1.96 * std_error
            ),
            width = 0.2
        ) +
        ggplot2::geom_hline(
            yintercept = 0,
            linetype = "dashed"
        ) +
        ggplot2::scale_fill_manual(
            values = c("TRUE" = "blue4", "FALSE" = "grey60"),
            name = "p < 0.05"
        ) +
        ggplot2::labs(
            title = "WQS Index Effect",
            x = NULL,
            y = "Effect Estimate"
        ) +
        ggplot2::theme_bw() +
        ggplot2::coord_flip() +
        ggplot2::theme(legend.position = "bottom")
}

# ---- BKMR plots ----

#' Plot BKMR results
#' @keywords internal
#' @noRd
.plot_bkmr <- function(
  results,
  plot_type,
  ...
) {
    switch(plot_type,
        pips = .plot_bkmr_pips(results),
        univariate = .plot_bkmr_univariate(results, ...),
        overall = .plot_bkmr_overall(results, ...),
        interaction = .plot_bkmr_interaction(results, ...),
        stop("plot_type must be 'pips', 'univariate', 'overall', or 'interaction'")
    )
}

#' Plot BKMR PIPs
#' @keywords internal
#' @noRd
.plot_bkmr_pips <- function(results) {
    if (is.null(results$weights)) {
        stop("No PIPs available. Run BKMR with varsel = TRUE.")
    }

    plot_data <- results$weights |>
        dplyr::arrange(dplyr::desc(pip)) |>
        dplyr::mutate(
            exposure = forcats::fct_reorder(exposure, pip),
            important = pip > 0.5
        )

    ggplot2::ggplot(
        plot_data,
        ggplot2::aes(
            x = pip,
            y = exposure,
            fill = important
        )
    ) +
        ggplot2::geom_col() +
        ggplot2::geom_vline(
            xintercept = 0.5,
            linetype = "dashed",
            color = "midnightblue"
        ) +
        ggplot2::scale_fill_manual(
            values = c(
                "TRUE" = "magenta4",
                "FALSE" = "grey60"
            ),
            name = "PIP > 0.5"
        ) +
        ggplot2::scale_x_continuous(limits = c(0, 1)) +
        ggplot2::labs(
            title = "BKMR Posterior Inclusion Probabilities",
            x = "Posterior Inclusion Probability",
            y = "Exposure"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid.major.y = ggplot2::element_blank(),
            legend.position = "bottom"
        )
}

#' Plot BKMR univariate exposure-response
#' @keywords internal
#' @noRd
.plot_bkmr_univariate <- function(
  results,
  sel = NULL,
  ...
) {
    if (is.null(results$fit)) {
        stop("BKMR fit object required. Run run_mixture_analysis() with keep_fit = TRUE.")
    }

    .check_suggested("bkmr")
    fit <- results$fit

    pred_resp <- bkmr::PredictorResponseUnivar(fit,
        sel = sel,
        ...
    )

    ggplot2::ggplot(
        pred_resp,
        ggplot2::aes(
            x = z,
            y = est,
            ymin = est - 1.96 * se,
            ymax = est + 1.96 * se
        )
    ) +
        ggplot2::geom_ribbon(alpha = 0.3, fill = "blue4") +
        ggplot2::geom_line(color = "blue4", linewidth = 1) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::facet_wrap(~variable, scales = "free_x") +
        ggplot2::labs(
            title = "BKMR Univariate Exposure-Response",
            caption = "Effect of each exposure holding others at median",
            x = "Exposure Level (quantile)",
            y = "Estimated Effect"
        ) +
        ggplot2::theme_bw()
}

#' Plot BKMR overall mixture effect
#' @keywords internal
#' @noRd
.plot_bkmr_overall <- function(
  results,
  qs = seq(0.1, 0.9, by = 0.1),
  ...
) {
    if (is.null(results$fit)) {
        stop("BKMR fit object required. Run run_mixture_analysis() with keep_fit = TRUE.")
    }

    .check_suggested("bkmr")
    fit <- results$fit

    overall <- bkmr::OverallRiskSummaries(fit, qs = qs, ...)

    ggplot2::ggplot(
        overall,
        ggplot2::aes(
            x = quantile,
            y = est,
            ymin = est - 1.96 * sd,
            ymax = est + 1.96 * sd
        )
    ) +
        ggplot2::geom_ribbon(alpha = 0.3, fill = "blue4") +
        ggplot2::geom_line(color = "blue4", linewidth = 1) +
        ggplot2::geom_point(color = "blue4", size = 2) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::labs(
            title = "BKMR Overall Mixture Effect",
            caption = "Effect of increasing all exposures together vs. median",
            x = "Quantile of All Exposures",
            y = "Estimated Effect"
        ) +
        ggplot2::theme_bw()
}

#' Plot BKMR interactions
#' @keywords internal
#' @noRd
.plot_bkmr_interaction <- function(
  results,
  ...
) {
    if (is.null(results$fit)) {
        stop("BKMR fit object required. Run run_mixture_analysis() with keep_fit = TRUE.")
    }

    .check_suggested("bkmr")
    fit <- results$fit

    interactions <- bkmr::PredictorResponseBivar(fit, ...)

    ggplot2::ggplot(
        interactions,
        ggplot2::aes(
            x = z1,
            y = est,
            color = factor(z2)
        )
    ) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::facet_wrap(~ variable1 + variable2, scales = "free") +
        ggplot2::scale_color_viridis_d(name = "Level of\nVariable 2") +
        ggplot2::labs(
            title = "BKMR Bivariate Exposure-Response",
            caption = "Interaction between pairs of exposures",
            x = "Level of Variable 1",
            y = "Estimated Effect"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "right")
}
