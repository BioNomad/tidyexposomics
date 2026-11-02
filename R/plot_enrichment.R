#' Plot Enrichment Results from exposomicset
#'
#' Visualize enrichment results stored in a `MultiAssayExperiment` object.
#' Supports dotplots, heatmaps, cnetplots, networks,
#' and multi-panel summary plots.
#'
#' @param exposomicset A `MultiAssayExperiment` object with enrichment results
#' added via `run_enrichment()`.
#' @param feature_type Character; one of `"degs"`, `"degs_robust"`,
#' `"omics"`, `"factor_features"`,
#'   `"degs_cor"`, `"omics_cor"`, or `"factor_features_cor"`.
#'   Defines which enrichment results to use.
#' @param plot_type Type of plot to generate. One of `"dotplot"`,
#'  `"cnet"`, `"network"`, `"heatmap"`, or `"summary"`.
#'
#' @param top_n Integer; number of top `go_group`s to include
#' (used in `"dotplot"`). Default is `5`.
#' @param n_per_group Integer; number of terms per group to plot
#' (used in `"dotplot"`). Default is `5`.
#' @param add_top_genes Logical; if `TRUE`, appends top shared genes to
#'  dotplot facets. Default is `TRUE`.
#' @param top_n_genes Integer; number of top genes to show in each group
#' (used in `"dotplot"`). Default is `5`.
#'
#' @param fill_data Optional data frame used to color heatmap tiles and cnet
#'   gene nodes. Needs `exp_name`, `feature`, and a numeric column named by
#'   `fill_col`. Build it from any source (differential abundance,
#'   exposure-omics association, or your own scores). If `NULL`, tiles and
#'   nodes are grey.
#' @param fill_col Name of the numeric column in `fill_data`.
#'   Default is `"estimate"`.
#' @param fill_label Legend title. Defaults to `fill_col`.
#'
#' @param overlap_thresh Numeric; Jaccard threshold for edges in
#' the network plot. Default is `0.2`.
#' @param node_radius Numeric; node size in network plot. Default is `0.2`.
#' @param pie_colors Optional named vector of colors for pie charts
#'  (network and cnet).
#' @param label_top_n Integer; number of top nodes to label in network.
#'  Default is `NULL`.
#' @param label_colour Color of node labels in network. Default is `"black"`.
#' @param net_facet_by Column used to facet the network plot
#' (e.g., `"category"`). Default is `NULL`.
#'
#' @param max_terms Integer; max number of terms to include in the cnet plot.
#'  Default is `30`.
#' @param node_size Numeric; base node size for cnet plot. Default is `1`.
#' @param term_node_correction Scaling factor for term nodes in cnet plot.
#' Default is `0.2`.
#' @param gene_node_correction Scaling factor for gene nodes in cnet plot.
#' Default is `3`.
#'
#' @param go_groups Optional character vector of GO group names to subset
#' enrichment results (all plots). Required for `"heatmap"`.
#' @param layout_algo Graph layout algorithm to use in `"network"` and `"cnet"`
#'  plots. Default is `"fr"`.
#' @param edge_alpha Transparency of network/cnet plot edges.
#'  Default is `0.3`.
#' @param label_size Font size for labels in network and cnet plots.
#'  Default is `3`.
#' @param feature_col Column in the feature metadata holding the gene
#'  identifiers used in `run_enrichment()`. Must match the `feature_col`
#'  used there. Used to map `fill_data` onto enriched genes.
#'  Default is `"feature"`.
#'
#' @return A `ggplot` or `patchwork` object corresponding to the
#' requested plot type.
#'
#' @details
#' This function visualizes results from `run_enrichment()`
#'  using one of several plot types:
#' \itemize{
#'   \item `"dotplot"`: Enrichment terms grouped by GO group, colored by
#'   significance.
#'   \item `"heatmap"`: Term - gene matrix with optional fill and
#'   shared gene highlighting.
#'   \item `"network"`: Graph of term overlap based on shared genes,
#'   faceted by metadata if desired.
#'   \item `"cnet"`: Gene - term bipartite graph with gene values
#'   and term pie slices.
#'   \item `"summary"`: Multi-panel figure with GO group ridgeplots,
#'   gene counts, and Venn diagram.
#' }
#'
#' For `"heatmap"` and `"cnet"`, `fill_data` must have one row per
#' `exp_name` x `feature`. When a gene symbol maps to several features within
#' an experiment, the value with the largest magnitude is kept.
#'
#' @examples
#' # create example data
#' mae <- make_example_data(
#'     n_samples = 30,
#'     return_mae = TRUE
#' )
#'
#' # perform differential abundance analysis
#' mae <- run_differential_abundance(
#'     exposomicset = mae,
#'     formula = ~ smoker + sex,
#'     abundance_col = "counts",
#'     method = "limma_voom",
#'     action = "add"
#' )
#'
#' # perform enrichment analysis
#' mae <- run_enrichment(
#'     exposomicset = mae,
#'     feature_type = "degs",
#'     feature_col = "symbol",
#'     species = "goa_human",
#'     deg_logfc_threshold = log2(1),
#'     deg_pval_col = "P.Value",
#'     deg_pval_threshold = 0.2,
#'     action = "add"
#' )
#'
#' # create an enrichment plot
#' enr_plot <- plot_enrichment(
#'     exposomicset = mae,
#'     feature_type = "degs",
#'     plot_type = "network"
#' )
#'
#' @importFrom MultiAssayExperiment metadata
#' @importFrom purrr pluck
#' @export
plot_enrichment <- function(
  exposomicset,
  feature_type = c(
      "degs",
      "degs_robust",
      "omics",
      "factor_features",
      "degs_cor",
      "omics_cor",
      "factor_features_cor"
  ),
  plot_type = c(
      "dotplot",
      "cnet",
      "network",
      "heatmap",
      "summary"
  ),
  # Dotplot arguments
  top_n = 5,
  n_per_group = 5,
  add_top_genes = TRUE,
  top_n_genes = 5,
  # Fill arguments (heatmap and cnet)
  fill_data = NULL,
  fill_col = "estimate",
  fill_label = NULL,
  # Network arguments
  overlap_thresh = 0.2,
  node_radius = 0.2,
  pie_colors = NULL,
  label_top_n = NULL,
  label_colour = "black",
  net_facet_by = NULL,
  # Cnet arguments
  max_terms = 30,
  node_size = 1,
  term_node_correction = 0.2,
  gene_node_correction = 3,
  # Shared arguments
  go_groups = NULL,
  layout_algo = "fr",
  edge_alpha = 0.3,
  label_size = 3,
  feature_col = "feature"
) {
    # Check that the user ran the enrichment analysis
    if (!"enrichment" %in% names(MultiAssayExperiment::metadata(exposomicset))) {
        stop("Please run `run_enrichment` first.")
    }

    feature_type <- match.arg(feature_type)
    plot_type <- match.arg(plot_type)

    enr_res <- exposomicset@metadata |>
        purrr::pluck("enrichment") |>
        purrr::pluck(feature_type)

    fill <- if (plot_type %in% c("heatmap", "cnet")) {
        .resolve_fill_data(
            exposomicset = exposomicset,
            feature_col = feature_col,
            fill_data = fill_data,
            fill_col = fill_col,
            fill_label = fill_label
        )
    }

    switch(plot_type,
        "dotplot" = {
            .plot_dotplot_enrichment(
                enr_res,
                top_n = top_n,
                n_per_group = n_per_group,
                add_top_genes = add_top_genes,
                top_n_genes = top_n_genes,
                go_groups = go_groups
            )
        },
        "cnet" = {
            .plot_cnet_enrichment(
                enr_res,
                fill = fill,
                layout_algo = layout_algo,
                go_groups = go_groups,
                max_terms = max_terms,
                edge_alpha = edge_alpha,
                node_size = node_size,
                label_size = label_size,
                term_node_correction = term_node_correction,
                gene_node_correction = gene_node_correction
            )
        },
        "heatmap" = {
            .plot_heatmap_enrichment(
                enr_res,
                fill = fill,
                go_groups = go_groups
            )
        },
        "network" = {
            .plot_network_enrichment(
                enr_res = enr_res,
                feature_type = feature_type,
                go_groups = go_groups,
                overlap_thresh = overlap_thresh,
                layout_algo = layout_algo,
                node_radius = node_radius,
                edge_alpha = edge_alpha,
                pie_colors = pie_colors,
                label_top_n = label_top_n,
                label_size = label_size,
                label_colour = label_colour,
                net_facet_by = net_facet_by
            )
        },
        "summary" = {
            .plot_summary_enrichment(enr_res)
        }
    )
}

# --- Enrichment DotPlot --------
#' @noRd
#'
#' @title Internal: Dotplot of Enrichment Terms by GO Group
#'
#' @description
#' Generates a dotplot of enriched terms colored by -log10(p) and
#' sized by the number of genes, optionally grouped by `go_group` and
#'  annotated with top genes.
#'
#' @param enr_res A data frame of enrichment results.
#' @param top_n Number of top GO groups to plot based on a scoring heuristic.
#' @param n_per_group Number of terms per group to show.
#' @param add_top_genes Logical; whether to append top genes to GO group labels.
#' @param top_n_genes Number of top genes to show in label.
#' @param go_groups Optional subset of GO groups to display.
#'
#' @return A `ggplot` dotplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_point labs scale_color_gradient facet_grid theme element_text
#' @importFrom dplyr filter group_by mutate reframe arrange slice_head ungroup inner_join pull bind_cols
#' @importFrom tibble rownames_to_column
#' @importFrom purrr map
#' @importFrom ggpubr theme_pubr rotate_x_text
#'
#' @keywords internal
.plot_dotplot_enrichment <- function(
  enr_res,
  top_n,
  n_per_group,
  add_top_genes,
  top_n_genes,
  go_groups
) {
    # require(ggplot2)
    .check_suggested(pkg = "forcats")

    if (nrow(enr_res) < 1) {
        return(NULL)
    }

    if (!is.null(go_groups)) {
        enr_res <- enr_res |>
            dplyr::filter(go_group %in% go_groups)
    } else {
        enr_res <- enr_res |>
            dplyr::filter(go_group %in% c(
                enr_res |>
                    dplyr::group_by(go_group) |>
                    dplyr::mutate(score_tmp = (-log10(p_adjust) * n_with_sel)) |>
                    dplyr::reframe(score = mean(score_tmp)) |>
                    dplyr::arrange(desc(score)) |>
                    dplyr::slice_head(n = top_n) |>
                    dplyr::pull(go_group)
            ))
    }

    go_group_genes_df <- enr_res |>
        (\(df) split(df, df$go_group))() |>
        purrr::map(~ .x |>
            dplyr::pull(ids) |>
            (\(x) strsplit(x, ","))() |>
            unlist() |>
            table() |>
            sort() |>
            tail(n = top_n_genes) |>
            names() |>
            (\(genes) {
                # Split into groups of 5 and add line breaks
                gene_chunks <- split(genes, ceiling(seq_along(genes) / 5))
                paste(
                    vapply(
                        gene_chunks,
                        function(chunk) paste(chunk, collapse = ", "),
                        FUN.VALUE = character(1)
                    ),
                    collapse = "\n"
                )
            })()) |>
        bind_cols() |>
        t() |>
        as.data.frame() |>
        setNames("gene_col") |>
        tibble::rownames_to_column("go_group")

    if (add_top_genes) {
        enr_res <- enr_res |>
            dplyr::inner_join(
                go_group_genes_df,
                by = "go_group"
            ) |>
            dplyr::mutate(go_group = gsub("_", " ", go_group))
    } else {
        enr_res <- enr_res |>
            dplyr::mutate(go_group = gsub("_", " ", go_group))
    }

    if ("category" %in% colnames(enr_res)) {
        enr_res |>
            dplyr::mutate(go_group = paste(go_group, gene_col, sep = "\n")) |>
            dplyr::group_by(category, exp_name, go_group) |>
            dplyr::arrange(p_adjust) |>
            dplyr::slice_head(n = n_per_group) |>
            dplyr::ungroup() |>
            ggplot(aes(
                x = forcats::fct_reorder(exp_name, -n_with_sel),
                y = forcats::fct_reorder(term_name, -log10(p_adjust)),
                size = n_with_sel,
                color = -log10(p_adjust)
            )) +
            geom_point() +
            ggpubr::theme_pubr(
                legend = "right",
                base_size = 10
            ) +
            facet_grid(go_group ~ category,
                space = "free",
                scales = "free"
            ) +
            ggpubr::rotate_x_text(angle = 45) +
            scale_color_gradient(low = "thistle1", high = "#8E0152") +
            theme(
                strip.text.x = element_text(
                    face = "bold.italic",
                    angle = 90
                ),
                strip.text.y = element_text(
                    face = "bold.italic",
                    angle = 0
                )
            ) +
            labs(
                x = "",
                y = "",
                title = "",
                size = "Number of Genes",
                color = expression(-Log[10](P)),
            )
    } else {
        enr_res |>
            dplyr::mutate(go_group = paste(go_group, gene_col, sep = "\n")) |>
            dplyr::group_by(exp_name, go_group) |>
            dplyr::arrange(p_adjust) |>
            dplyr::slice_head(n = n_per_group) |>
            dplyr::ungroup() |>
            ggplot(aes(
                x = forcats::fct_reorder(exp_name, -n_with_sel),
                y = forcats::fct_reorder(term_name, -log10(p_adjust)),
                size = n_with_sel,
                color = -log10(p_adjust)
            )) +
            geom_point() +
            ggpubr::theme_pubr(
                legend = "right",
                base_size = 10
            ) +
            facet_grid(go_group ~ ., space = "free", scales = "free") +
            ggpubr::rotate_x_text(angle = 45) +
            scale_color_gradient(low = "thistle1", high = "#8E0152") +
            theme(
                strip.text.x = element_text(
                    face = "bold.italic",
                    angle = 90
                ),
                strip.text.y = element_text(
                    face = "bold.italic",
                    angle = 0
                )
            ) +
            labs(
                x = "",
                y = "",
                title = "",
                size = "Number of Genes",
                color = expression(-Log[10](P)),
            )
    }
}

# --- Enrichment Heatmap -----
#' @noRd
#'
#' @title Internal: Heatmap of Enrichment Terms by Gene
#'
#' @description
#' Tile heatmap of enrichment terms by associated genes, optionally filled
#' by a resolved fill table.
#'
#' @param enr_res Enrichment results.
#' @param fill Output of `.resolve_fill_data()`, or `NULL` for grey tiles.
#' @param go_groups GO groups to include in the heatmap.
#'
#' @return A `ggplot` heatmap object.
#'
#' @importFrom ggplot2 ggplot aes theme_bw facet_grid theme element_text element_blank labs geom_tile
#' @importFrom dplyr filter distinct count mutate left_join if_else
#' @importFrom tidyr separate_rows
#'
#' @keywords internal
.plot_heatmap_enrichment <- function(
  enr_res,
  fill,
  go_groups
) {
    .check_suggested("forcats")

    if (nrow(enr_res) < 1) {
        return(NULL)
    }

    if (is.null(go_groups)) {
        stop("Please specify a `go_group`")
    }

    enr_long <- enr_res |>
        dplyr::filter(go_group %in% go_groups) |>
        tidyr::separate_rows("ids", sep = ",\\s*")

    shared_genes <- enr_long |>
        dplyr::distinct(exp_name, ids) |>
        dplyr::count(ids, name = "n_exp") |>
        dplyr::mutate(shared = n_exp > 1)

    if (!is.null(fill)) {
        enr_long <- enr_long |>
            dplyr::left_join(fill$data, by = c("exp_name", "ids"))
    }

    enr_long <- enr_long |>
        dplyr::left_join(shared_genes, by = "ids") |>
        dplyr::mutate(ids = dplyr::if_else(shared, paste0("*", ids), ids))

    p <- if ("category" %in% colnames(enr_long)) {
        enr_long |>
            ggplot2::ggplot(ggplot2::aes(
                x = forcats::fct_reorder2(ids, term_name, category),
                y = forcats::fct_reorder2(term_name, exp_name, category)
            )) +
            ggplot2::facet_grid(
                category ~ exp_name,
                scales = "free",
                space = "free"
            )
    } else {
        enr_long |>
            ggplot2::ggplot(ggplot2::aes(
                x = forcats::fct_reorder(ids, term_name),
                y = forcats::fct_reorder(term_name, exp_name)
            )) +
            ggplot2::facet_grid(
                ~exp_name,
                scales = "free",
                space = "free"
            )
    }

    p <- p +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(
                angle = 90,
                face = "italic"
            ),
            axis.ticks.x = ggplot2::element_blank(),
            strip.text.x = ggplot2::element_text(face = "bold.italic"),
            strip.text.y = ggplot2::element_text(face = "bold.italic"),
            plot.caption = ggplot2::element_text(face = "italic"),
            plot.title = ggplot2::element_text(face = "bold.italic")
        ) +
        ggplot2::labs(
            caption = "* Appears in more than one omic",
            x = "",
            y = "",
            title = paste(
                "Enrichment Term Heatmap for",
                paste(gsub("_", " ", go_groups), collapse = ", ")
            )
        )

    if (is.null(fill)) {
        p +
            ggplot2::geom_tile(fill = "grey30")
    } else {
        p +
            ggplot2::geom_tile(ggplot2::aes(fill = fill_value)) +
            .scale_value(fill, "fill")
    }
}

# --- Enrichment Network Plot -------
#' @noRd
#'
#' @title Internal: Network Plot of Enrichment Term Overlaps
#'
#' @description
#' Builds a graph from overlapping genes between terms and renders a
#'  circular or force-directed layout colored by experiment.
#'
#' @param enr_res Enrichment results.
#' @param feature_type Type of feature used.
#' @param go_groups Optional GO groups to subset.
#' @param overlap_thresh Minimum Jaccard overlap to draw edges.
#' @param layout_algo Graph layout algorithm (e.g., "fr").
#' @param node_radius Radius of node pie arcs.
#' @param edge_alpha Edge transparency.
#' @param pie_colors Optional named vector of colors for experiments.
#' @param label_top_n Number of top nodes to label by degree centrality.
#' @param label_size Font size of node labels.
#' @param label_colour Color of node labels.
#' @param net_facet_by Optional column to facet the network by.
#'
#' @return A `ggraph` object.
#'
#' @importFrom igraph graph_from_data_frame cluster_louvain
#' @importFrom ggplot2 facet_wrap theme element_text labs scale_fill_manual theme_void
#' @importFrom ggrepel geom_label_repel
#' @importFrom dplyr filter mutate select distinct group_by ungroup left_join top_n pull
#' @importFrom purrr map
#' @importFrom tidyr separate_rows
#' @importFrom stringr str_trim
#' @importFrom tibble as_tibble
#'
#' @keywords internal
.plot_network_enrichment <- function(
  enr_res = enr_res,
  feature_type = feature_type,
  go_groups = go_groups,
  overlap_thresh = overlap_thresh,
  layout_algo = layout_algo,
  node_radius = node_radius,
  edge_alpha = edge_alpha,
  pie_colors = pie_colors,
  label_top_n = label_top_n,
  label_size = label_size,
  label_colour = label_colour,
  net_facet_by = net_facet_by
) {
    # require(igraph)
    # require(tidygraph)
    # require(ggraph)
    .check_suggested(pkg = "ggraph")
    .check_suggested(pkg = "tidygraph")

    if (nrow(enr_res) < 1) {
        return(NULL)
    }

    # Get enrichment table
    enr <- enr_res

    if (!is.null(go_groups)) {
        enr <- enr |>
            dplyr::filter(go_group %in% go_groups)
    }

    # Build pairwise overlaps by term_name
    term_overlaps <- split(enr, enr$term_name) |>
        purrr::map(~ .x |>
            tidyr::separate_rows("ids", sep = ",") |>
            dplyr::pull(ids) |>
            stringr::str_trim() |>
            unique()) |>
        .get_pairwise_overlaps()

    # Build deduplicated edge list
    edge_df <- term_overlaps |>
        dplyr::mutate(
            source_std = pmin(source, target),
            target_std = pmax(source, target)
        ) |>
        dplyr::select(
            source = source_std,
            target = target_std,
            overlap
        ) |>
        dplyr::filter(overlap > overlap_thresh) |>
        dplyr::distinct()

    # Ensure unique node list by term_name
    node_tbl <- enr |>
        dplyr::select(term_name) |>
        dplyr::distinct() |>
        dplyr::mutate(name = term_name) # igraph expects "name"

    # Only keep nodes involved in edges
    keep_nodes <- unique(c(edge_df$source, edge_df$target))
    node_tbl <- node_tbl |> dplyr::filter(name %in% keep_nodes)

    # Build igraph object
    g <- igraph::graph_from_data_frame(
        edge_df,
        directed = FALSE,
        vertices = node_tbl
    ) |>
        tidygraph::as_tbl_graph()

    # Layout
    layout <- ggraph::create_layout(g, layout = layout_algo)

    # Compute pie slice positions for omics types
    node_omics <- enr |>
        dplyr::select(term_name, exp_name) |>
        dplyr::distinct() |>
        dplyr::group_by(term_name) |>
        dplyr::mutate(
            n = dplyr::n(),
            angle = 2 * pi / n,
            start = (row_number() - 1) * angle,
            end = row_number() * angle
        ) |>
        dplyr::ungroup()

    # Join layout with slice data
    arc_data <- layout |>
        as_tibble() |>
        dplyr::left_join(node_omics, by = c("name" = "term_name"))

    # Stop if join failed
    if (!"exp_name" %in% names(arc_data)) {
        stop("Join failed: 'exp_name' not found in arc_data.")
    }

    # Default pie colors
    if (is.null(pie_colors)) {
        omics_types <- unique(arc_data$exp_name)
        n_colors <- max(3, length(omics_types))
        col_options <- c(
            "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
            "#FFD92F", "#E5C494", "#B3B3B3"
        )
        # Repeat colors if necessary
        pie_colors <- setNames(
            rep(col_options, length.out = length(omics_types)),
            omics_types
        )
    }

    # Final plot
    p <- ggraph::ggraph(layout) +
        ggraph::geom_edge_link(aes(width = overlap),
            alpha = edge_alpha,
            color = "gray50"
        ) +
        ggraph::scale_edge_width(range = c(0.2, 2), guide = "none") +
        ggraph::geom_node_arc_bar(
            data = arc_data,
            aes(
                x0 = x, y0 = y, r0 = 0, r = node_radius,
                start = start, end = end,
                fill = exp_name
            ),
            color = "transparent"
        ) +
        scale_fill_manual(values = tidy_exp_pal) +
        theme_void() +
        labs(fill = "Experiment")

    if (!is.null(net_facet_by)) {
        # Add category info to node positions (layout)
        category_map <- enr |>
            dplyr::select(term_name, !!rlang::sym(net_facet_by)) |>
            dplyr::distinct()

        layout <- layout |>
            dplyr::left_join(category_map, by = c("name" = "term_name"))

        arc_data <- arc_data |>
            dplyr::left_join(category_map, by = c("name" = "term_name"))

        if (!net_facet_by %in% names(layout)) {
            stop(sprintf("Could not join `net_facet_by` by %s.", net_facet_by))
        }

        p <- p +
            facet_wrap(as.formula(paste("~", net_facet_by))) +
            theme(strip.text.x = element_text(
                face = "bold.italic",
                size = 12
            ))
    }

    if (!is.null(label_top_n)) {
        community <- igraph::cluster_louvain(g)

        # Compute layout + node metrics
        node_metrics <- tidygraph::activate(g, nodes) |>
            dplyr::mutate(
                degree = tidygraph::centrality_degree(),
                cluster = community$membership
            ) |>
            dplyr::as_tibble()

        #
        top_nodes <- node_metrics |>
            dplyr::group_by(cluster) |>
            dplyr::top_n(label_top_n, degree) |>
            dplyr::ungroup() |>
            dplyr::pull(name)

        layout_labels <- layout |>
            dplyr::filter(name %in% top_nodes)

        p <- p +
            ggrepel::geom_label_repel(
                data = layout_labels,
                aes(x = x, y = y, label = name),
                size = label_size,
                color = label_colour,
                force = 2, # stronger repel
                max.overlaps = Inf,
                segment.color = "grey30",
                segment.size = 0.2,
                box.padding = 0.3,
                point.padding = 0.2
            )
    }

    p
}

# --- Enrichment Cnetplot ------
#' @noRd
#'
#' @title Internal: Gene - Term Network (Cnetplot)
#'
#' @description
#' Bipartite graph between genes and enriched terms, with gene nodes colored
#' by a resolved fill table and term nodes drawn as pies by experiment.
#'
#' @param enr_res Enrichment results.
#' @param fill Output of `.resolve_fill_data()`, or `NULL` for grey genes.
#' @param layout_algo Graph layout algorithm. Default is `"fr"`.
#' @param go_groups Optional subset of GO groups.
#' @param max_terms Maximum number of terms when `go_groups` is `NULL`.
#' @param edge_alpha Transparency of edges.
#' @param node_size Base size for graph nodes.
#' @param label_size Size of text labels.
#' @param term_node_correction Scaling factor for term nodes.
#' @param gene_node_correction Scaling factor for gene nodes.
#'
#' @return A `ggraph` plot object.
#'
#' @importFrom igraph graph_from_data_frame
#' @importFrom ggplot2 aes theme_void labs theme element_text
#' @importFrom dplyr filter mutate select distinct left_join count slice_head pull bind_rows group_by ungroup add_count if_else
#' @importFrom tidyr separate_rows
#' @importFrom tibble as_tibble
#'
#' @keywords internal
.plot_cnet_enrichment <- function(
  enr_res,
  fill,
  layout_algo = "fr",
  go_groups = NULL,
  max_terms = 30,
  edge_alpha = 0.5,
  node_size = 1,
  label_size = 3,
  term_node_correction = 0.2,
  gene_node_correction = 3
) {
    if (nrow(enr_res) < 1) {
        return(NULL)
    }
    .check_suggested(pkg = "ggnewscale")
    .check_suggested(pkg = "ggraph")
    .check_suggested(pkg = "tidygraph")

    # Explode to gene by term by exp rows
    long_enr <- enr_res |>
        tidyr::separate_rows(ids, sep = ",\\s*") |>
        dplyr::select(
            term_name,
            go_group,
            ids,
            exp_name
        ) |>
        dplyr::distinct()

    # Filter GO groups or select top terms
    if (!is.null(go_groups)) {
        long_enr <- long_enr |>
            dplyr::filter(go_group %in% go_groups)
    } else {
        long_enr <- long_enr |>
            dplyr::filter(
                term_name %in% (
                    long_enr |>
                        dplyr::count(term_name, sort = TRUE) |>
                        dplyr::slice_head(n = max_terms) |>
                        dplyr::pull(term_name)
                )
            )
    }

    # Attach fill values (NA when no fill source) and node names
    long_enr <- long_enr |>
        (\(df) {
            if (is.null(fill)) {
                dplyr::mutate(df, fill_value = NA_real_)
            } else {
                dplyr::left_join(df, fill$data, by = c("exp_name", "ids"))
            }
        })() |>
        dplyr::mutate(
            gene_node = paste0("gene:", ids, "::", exp_name),
            term_node = paste0("term:", term_name)
        )

    # Build edges and nodes
    edges <- long_enr |>
        dplyr::select(from = gene_node, to = term_node)

    gene_meta <- long_enr |>
        dplyr::select(name = gene_node, fill_value) |>
        dplyr::distinct() |>
        dplyr::mutate(type = "gene")

    term_meta <- long_enr |>
        dplyr::select(name = term_node, exp_name) |>
        dplyr::distinct() |>
        dplyr::mutate(type = "term")

    # Graph layout
    layout <- igraph::graph_from_data_frame(
        edges,
        vertices = dplyr::bind_rows(term_meta, gene_meta) |>
            dplyr::distinct(name, .keep_all = TRUE),
        directed = FALSE
    ) |>
        ggraph::create_layout(layout = layout_algo)

    # Pie slices for term nodes
    term_arc_data <- long_enr |>
        dplyr::distinct(name = term_node, exp_name) |>
        dplyr::group_by(name) |>
        dplyr::mutate(
            n = dplyr::n(),
            angle = 2 * pi / n,
            start = (dplyr::row_number() - 1) * angle,
            end = dplyr::row_number() * angle
        ) |>
        dplyr::ungroup() |>
        dplyr::left_join(
            layout |>
                tibble::as_tibble() |>
                dplyr::select(name, x, y),
            by = "name"
        )

    # Gene labels with conditional exp_name
    gene_labels <- long_enr |>
        dplyr::distinct(ids, exp_name) |>
        dplyr::add_count(ids, name = "n_exps") |>
        dplyr::mutate(
            label = dplyr::if_else(
                n_exps > 1,
                paste0(ids, " (", exp_name, ")"),
                ids
            ),
            name = paste0("gene:", ids, "::", exp_name)
        ) |>
        dplyr::select(name, label)

    layout <- layout |>
        dplyr::left_join(gene_labels, by = "name") |>
        dplyr::mutate(
            label = ifelse(
                is.na(label),
                gsub("^.*?:", "", name),
                label
            )
        )

    gene_nodes <- layout |>
        dplyr::filter(type == "gene")

    # Edges and term pies, then reset colour scale before genes
    p <- ggraph::ggraph(layout) +
        ggraph::geom_edge_link(
            alpha = edge_alpha,
            color = "grey70"
        ) +
        ggraph::geom_node_arc_bar(
            data = term_arc_data,
            ggplot2::aes(
                x0 = x,
                y0 = y,
                r0 = 0,
                r = node_size * term_node_correction,
                start = start,
                end = end,
                fill = exp_name
            ),
            color = "transparent"
        ) +
        scale_fill_tidy_exp(name = "Experiment") +
        ggnewscale::new_scale_color()

    # Gene nodes: coloured by fill values, or grey when no fill supplied
    p <- if (is.null(fill)) {
        p +
            ggraph::geom_node_point(
                data = gene_nodes,
                ggplot2::aes(x = x, y = y),
                color = "grey60",
                size = node_size * gene_node_correction
            )
    } else {
        p +
            ggraph::geom_node_point(
                data = gene_nodes,
                ggplot2::aes(
                    x = x,
                    y = y,
                    color = fill_value
                ),
                size = node_size * gene_node_correction
            ) +
            .scale_value(fill, "colour")
    }

    p +
        ggraph::geom_node_label(
            data = layout,
            ggplot2::aes(
                x = x,
                y = y,
                label = label,
                fontface = ifelse(type == "gene", "italic", "bold")
            ),
            size = label_size,
            repel = TRUE
        ) +
        ggplot2::theme_void() +
        ggplot2::labs(title = "Gene-Term Network") +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold.italic")
        )
}


# --- Enrichment Summary -------
#' @noRd
#'
#' @title Internal: Multi-panel Summary of Enrichment Results
#'
#' @description
#' Assembles barplots, ridgeplots, and a Venn diagram to summarize enrichment
#' results across GO groups and experiments.
#'
#' @param enr_res Enrichment result data frame.
#'
#' @return A `patchwork` object combining multiple ggplots.
#'
#' @importFrom ggplot2 ggplot aes geom_col coord_flip labs theme_minimal theme element_blank element_text element_line element_rect element_text scale_y_continuous scale_fill_gradient scale_size_continuous scale_fill_manual
#' @importFrom dplyr mutate distinct count filter left_join arrange pull select group_by ungroup
#' @importFrom tidyr separate_rows
#' @importFrom purrr map
#'
#' @keywords internal
.plot_summary_enrichment <- function(
  enr_res
) {
    if (nrow(enr_res) < 1) {
        return(NULL)
    }
    .check_suggested("ggridges")
    .check_suggested("ggvenn")
    .check_suggested("patchwork")

    # Clean up enr_res so that go_groups have no underscore
    enr_res <- enr_res |>
        dplyr::mutate(go_group = gsub("_", " ", go_group))

    # --- Barplot: # of enriched terms per exp_name
    bar_exp <- enr_res |>
        dplyr::distinct(
            term_name,
            exp_name
        ) |>
        dplyr::count(exp_name) |>
        ggplot(aes(
            x = reorder(exp_name, n),
            y = n,
            fill = exp_name
        )) +
        geom_col() +
        coord_flip() +
        scale_fill_tidy_exp() +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
        labs(
            y = "No. Enriched Terms",
            x = ""
        ) +
        theme_minimal() +
        theme(legend.position = "none")

    # --- Barplot: # of enriched terms per category
    if ("category" %in% colnames(enr_res)) {
        bar_cat <- enr_res |>
            dplyr::distinct(term_name, category) |>
            dplyr::count(category) |>
            ggplot(aes(
                x = reorder(category, n),
                y = n,
                fill = category
            )) +
            geom_col() +
            coord_flip() +
            scale_fill_tidy_exp(rev = TRUE) +
            labs(
                y = "No. Enriched Terms",
                x = "Category"
            ) +
            theme_minimal() +
            theme(legend.position = "none")
    }
    # --- Combination Plot ----
    if ("go_group" %in% colnames(enr_res)) {
        # Count how many terms per group
        group_term_counts <- enr_res |>
            dplyr::filter(!is.na(go_group)) |>
            dplyr::count(go_group, name = "n_terms") |>
            dplyr::filter(n_terms > 0)

        valid_groups <- group_term_counts$go_group

        # Order go_groups by number of terms
        group_levels <- group_term_counts |>
            dplyr::arrange(dplyr::desc(n_terms)) |>
            dplyr::pull(go_group)

        # Count genes per valid group (for dot size + ridge fill)
        go_counts <- enr_res |>
            dplyr::filter(go_group %in% valid_groups) |>
            tidyr::separate_rows(ids, sep = ",\\s*") |>
            dplyr::distinct(go_group, ids) |>
            dplyr::count(go_group, name = "gene_count")

        # Prepare enrichment table with metadata
        enr_go <- enr_res |>
            dplyr::filter(go_group %in% valid_groups) |>
            dplyr::left_join(go_counts, by = "go_group") |>
            dplyr::mutate(go_group = factor(go_group, levels = group_levels))

        # --- (Optional) Barplot: number of unique categories per group
        if ("category" %in% colnames(enr_res)) {
            bar_cat_group <- enr_res |>
                dplyr::filter(go_group %in% valid_groups) |>
                dplyr::distinct(go_group, category) |>
                dplyr::count(go_group, category, name = "n_cat") |>
                dplyr::mutate(go_group = factor(go_group, levels = group_levels))

            bar_category <- bar_cat_group |>
                ggplot(aes(x = n_cat, y = go_group, fill = category)) +
                geom_col(position = "stack") +
                scale_fill_tidy_exp(rev = TRUE) +
                scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
                labs(
                    x = "No. Categories",
                    y = NULL,
                    fill = "Category"
                ) +
                theme_minimal() +
                theme(
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    panel.grid = element_blank(),
                    legend.position = "right"
                )
        }


        # --- Dotplot of total gene count per GO group
        dot_gene <- go_counts |>
            dplyr::mutate(go_group = factor(go_group, levels = group_levels)) |>
            ggplot(aes(x = 1, y = go_group, size = gene_count)) +
            geom_point(color = "midnightblue", stroke = 0.3) +
            scale_size_continuous(range = c(2, 8)) +
            labs(x = NULL, y = NULL, size = "No. Genes") +
            coord_cartesian(clip = "off") +
            theme_minimal() +
            theme(
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank(),
                panel.grid = element_blank(),
                legend.position = "right"
            )

        # --- Ridgeplot of -log10(padj), filled by gene count
        ridge_p <- enr_go |>
            ggplot(aes(x = -log10(padj), y = go_group, fill = gene_count)) +
            ggridges::geom_density_ridges(scale = 1.2, alpha = 0.6) +
            scale_fill_gradient(low = "thistle", high = "midnightblue") +
            labs(x = expression("-Log"[10] * "P"), y = NULL, fill = "No. Genes") +
            ggridges::theme_ridges() +
            theme(
                legend.position = "none",
                axis.title.x = element_text(size = 12)
            )

        # --- Barplot: number of terms per group by exp_name
        group_exp_terms <- enr_res |>
            dplyr::filter(go_group %in% valid_groups) |>
            dplyr::distinct(go_group, exp_name, term_name) |>
            dplyr::count(go_group, exp_name, name = "n_terms") |>
            dplyr::mutate(go_group = factor(go_group, levels = group_levels))

        bar_expname <- group_exp_terms |>
            ggplot(aes(x = n_terms, y = go_group, fill = exp_name)) +
            geom_col(position = "stack") +
            scale_fill_tidy_exp() +
            scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
            labs(x = "No. Terms", y = NULL, fill = "Experiment") +
            theme_minimal() +
            theme(
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                legend.position = "right"
            )

        # --- Combine panels
        if (exists("bar_category")) {
            ridge <- bar_category + dot_gene + ridge_p + bar_expname +
                patchwork::plot_layout(widths = c(1.5, 1, 3.5, 2), guides = "collect")
        } else {
            ridge <- dot_gene + ridge_p + bar_expname +
                patchwork::plot_layout(widths = c(1, 3.5, 2), guides = "collect")
        }
    }


    # --- Venn diagram of term overlap
    venn_data <- enr_res |>
        dplyr::select(exp_name, term_name) |>
        dplyr::distinct() |>
        split(~exp_name) |>
        purrr::map(~ .x$term_name)

    venn <- ggvenn::ggvenn(
        venn_data,
        fill_color = tidy_exp_pal,
        stroke_color = "transparent",
        set_name_size = 4,
        show_percentage = TRUE
    ) +
        theme(
            plot.margin = margin(0, 0, 0, 0),
            plot.background = element_blank(),
            panel.background = element_blank()
        )

    # --- Layout assembly
    layout_left <- ridge
    layout_right <- (bar_exp / venn) +
        patchwork::plot_layout(heights = c(1, 2.5))

    final_plot <- (layout_left | layout_right) +
        patchwork::plot_annotation(title = "Enrichment Summary") &
        theme(plot.title = element_text(face = "bold.italic"))

    return(final_plot)
}

# --- Resolve Fill Values for Enrichment Plots ------------
#' @noRd
#'
#' @title Internal: Resolve a standardized fill table for enrichment plots
#'
#' @description
#' Maps a user-supplied table of per-feature values onto the gene identifiers
#' used in the enrichment results. Returns a list with `data`
#' (columns `exp_name`, `ids`, `fill_value`), a legend `label`, and a
#' `diverging` flag, or `NULL` when no fill data was supplied.
#'
#' @param exposomicset A `MultiAssayExperiment` object.
#' @param feature_col Column in the feature metadata holding the gene
#'   identifiers used in `run_enrichment()`.
#' @param fill_data Data frame with `exp_name`, `feature`, and `fill_col`.
#' @param fill_col Name of the numeric column in `fill_data`.
#' @param fill_label Legend title.
#'
#' @importFrom dplyr select inner_join summarise all_of
#' @importFrom rlang %||%
#' @keywords internal
.resolve_fill_data <- function(
  exposomicset,
  feature_col,
  fill_data,
  fill_col,
  fill_label
) {
    if (is.null(fill_data)) {
        return(NULL)
    }

    required <- c("exp_name", "feature", fill_col)
    if (!is.data.frame(fill_data) || !all(required %in% names(fill_data))) {
        stop(
            "`fill_data` must be a data frame with columns: ",
            paste(required, collapse = ", ")
        )
    }

    if (!is.numeric(fill_data[[fill_col]])) {
        stop("`", fill_col, "` must be numeric.")
    }

    fill_df <- fill_data |>
        dplyr::select(
            exp_name,
            feature,
            fill_value = dplyr::all_of(fill_col)
        ) |>
        dplyr::inner_join(
            pivot_feature(exposomicset) |>
                dplyr::select(
                    exp_name = .exp_name,
                    feature = .feature,
                    ids = dplyr::all_of(feature_col)
                ),
            by = c("exp_name", "feature")
        ) |>
        dplyr::summarise(
            fill_value = fill_value[which.max(abs(fill_value))],
            .by = c(exp_name, ids)
        )

    list(
        data = fill_df,
        label = fill_label %||% fill_col,
        diverging = any(fill_df$fill_value < 0, na.rm = TRUE)
    )
}

# --- Fill / Colour Scale for Resolved Fill ------------
#' @noRd
#' @keywords internal
.scale_value <- function(
  fill,
  aesthetic = c("fill", "colour")
) {
    aesthetic <- match.arg(aesthetic)

    scale_fn <- if (fill$diverging) {
        switch(aesthetic,
            fill = ggplot2::scale_fill_gradient2,
            colour = ggplot2::scale_colour_gradient2
        )
    } else {
        switch(aesthetic,
            fill = ggplot2::scale_fill_gradient,
            colour = ggplot2::scale_colour_gradient
        )
    }

    scale_args <- if (fill$diverging) {
        list(
            low = "blue4",
            mid = "white",
            high = "red4",
            midpoint = 0
        )
    } else {
        list(
            low = "thistle",
            high = "midnightblue"
        )
    }

    do.call(
        scale_fn,
        c(list(name = fill$label, na.value = "grey60"), scale_args)
    )
}
