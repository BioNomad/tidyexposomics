#' Plot Enrichment Results from ExpOmicSet
#'
#' Visualize enrichment results stored in a `MultiAssayExperiment` object.
#' Supports dotplots, heatmaps, cnetplots, networks,
#' and multi-panel summary plots.
#'
#' @param expomicset A `MultiAssayExperiment` object with enrichment results
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
#' @param heatmap_fill Logical; whether to fill tiles by logFC in the heatmap.
#' Default is `TRUE`.
#' @param logfc_thresh Numeric; log2 fold change threshold for filtering
#' (heatmap only). Default is `log2(1)`.
#' @param pval_col Column name of the p-value used for filtering in `"degs"`
#'  heatmap. Default is `"P.Value"`.
#' @param pval_thresh Threshold for `pval_col` (heatmap only).
#' Default is `0.05`.
#' @param score_metric Column for stability score
#' (used in `"degs_robust"` heatmap). Default is `"stability_score"`.
#' @param score_thresh Numeric; threshold for `score_metric`
#' (heatmap only). Default is `NULL`.
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
#' enrichment results (all plots).
#' @param layout_algo Graph layout algorithm to use in `"network"` and `"cnet"`
#'  plots. Default is `"fr"`.
#' @param edge_alpha Transparency of network/cnet plot edges.
#'  Default is `0.3`.
#' @param label_size Font size for labels in network and cnet plots.
#'  Default is `3`.
#' @param feature_col Column name used to join gene-level metadata.
#' Default is `"feature"`.
#' @param logfc_col Column name used for log2 fold change values.
#' Default is `"logFC"`.
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
#'   \item `"heatmap"`: Term–gene matrix with optional logFC fill and
#'   shared gene highlighting.
#'   \item `"network"`: Graph of term overlap based on shared genes,
#'   faceted by metadata if desired.
#'   \item `"cnet"`: Gene–term bipartite graph with gene logFC values
#'   and term pie slices.
#'   \item `"summary"`: Multi-panel figure with GO group ridgeplots,
#'   gene counts, and Venn diagram.
#' }
#'
#' @examples
#' \dontrun{
#' # Dotplot of top GO groups
#' plot_enrichment(expomicset, feature_type = "degs",
#'  plot_type = "dotplot")
#'
#' # Heatmap for selected groups
#' plot_enrichment(expomicset,
#' feature_type = "degs_robust",
#'  plot_type = "heatmap",
#'  go_groups = c("Group_1", "Group_2"))
#'
#' # Gene-term cnetplot
#' plot_enrichment(expomicset,
#' feature_type = "degs",
#'  plot_type = "cnet")
#'
#' # Summary multi-panel
#' plot_enrichment(expomicset,
#' feature_type = "degs",
#' plot_type = "summary")
#' }
#'
#' @export
plot_enrichment <- function(
    expomicset,
    feature_type = c("degs",
                     "degs_robust",
                     "omics",
                     "factor_features",
                     "degs_cor",
                     "omics_cor",
                     "factor_features_cor"),
    plot_type = c("dotplot",
                  "cnet",
                  "network",
                  "heatmap",
                  "summary"),

    # Dotplot arguments
    top_n=5,
    n_per_group=5,
    add_top_genes=TRUE,
    top_n_genes=5,

    # Heatmap arguments
    heatmap_fill = TRUE,
    logfc_thresh = log2(1),
    pval_col = "P.Value",
    pval_thresh = 0.05,
    score_metric = "stability_score",
    score_thresh = NULL,

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
    term_node_correction=0.2,
    gene_node_correction=3,

    # Shared arguments
    go_groups = NULL,
    layout_algo = "fr",
    edge_alpha = 0.3,
    label_size = 3,
    feature_col = "feature",
    logfc_col = "logFC"


){

  # Check that the user ran the enrichment analysis
  if(!"enrichment" %in% names(MultiAssayExperiment::metadata(expomicset))){
    stop("Please run `run_enrichment` first.")
  }

  # match arguments
  feature_type <-  match.arg(feature_type)
  plot_type <- match.arg(plot_type)

  enr_res <- expomicset@metadata |>
    purrr::pluck("enrichment") |>
    purrr::pluck(feature_type)

  plot <- switch(
    plot_type,
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
    "cnet"= {
      .plot_cnet_enrichment(
        enr_res,
        expomicset = expomicset,
        feature_type = feature_type,
        layout_algo = layout_algo,
        go_groups = go_groups,
        max_terms = max_terms,
        edge_alpha = edge_alpha,
        node_size = node_size,
        label_size = label_size,
        term_node_correction=term_node_correction,
        gene_node_correction=gene_node_correction,
        feature_col = feature_col,
        logfc_col = logfc_col
      )
    },
    "heatmap"= {
      .plot_heatmap_enrichment(
        enr_res = enr_res,
        expomicset = expomicset,
        feature_type = feature_type,
        score_metric = score_metric,
        score_thresh = score_thresh,
        go_groups = go_groups,
        heatmap_fill = heatmap_fill,
        feature_col = feature_col,
        logfc_col = logfc_col,
        logfc_thresh = logfc_thresh,
        pval_col = pval_col,
        pval_thresh = pval_thresh
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
      .plot_summary_enrichment(
        enr_res
      )
    }
  )

  return(plot)

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
#' @importFrom ggplot2 ggplot aes geom_point labs scale_color_gradient
#' facet_grid theme element_text
#' @importFrom dplyr filter group_by mutate reframe arrange slice_head
#' ungroup inner_join pull
#' @importFrom tibble rownames_to_column
#' @importFrom purrr map
#' @importFrom forcats fct_reorder
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
){
  #require(ggplot2)

  if(!is.null(go_groups)){
    enr_res <- enr_res |>
      dplyr::filter(go_group %in% go_groups)
  }else{
    enr_res <- enr_res |>
      dplyr::filter(go_group %in% c(
        enr_res |>
          dplyr::group_by(go_group) |>
          dplyr::mutate(score_tmp= (-log10(p_adjust) * n_with_sel)) |>
          dplyr::reframe(score=mean(score_tmp)) |>
          dplyr::arrange(desc(score)) |>
          dplyr::slice_head(n=top_n) |>
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
                 })()
    ) |>
    as.data.frame() |>
    t() |>
    as.data.frame() |>
    setNames("gene_col") |>
    tibble::rownames_to_column("go_group")

  # go_group_genes_df <- enr_res |>
  #   (\(df) split(df, df$go_group))() |>
  #   purrr::map(~ .x |>
  #                dplyr::pull(ids) |>
  #                (\(x) strsplit(x, ","))() |>
  #                unlist() |>
  #                table() |>
  #                sort() |>
  #                tail(n = top_n_genes) |>
  #                names() |>
  #                (\(genes) {
  #                  # Split into groups of 5 and add line breaks
  #                  gene_chunks <- split(genes,
  #                                       ceiling(seq_along(genes) / 5))
  #                  paste(sapply(gene_chunks,
  #                               function(chunk) paste(
  #                                 chunk, collapse = ", ")),
  #                        collapse = "\n")
  #                })()
  #   ) |>
  #   as.data.frame() |>
  #   t() |>
  #   as.data.frame() |>
  #   setNames("gene_col") |>
  #   tibble::rownames_to_column("go_group")

  if(add_top_genes){
    enr_res <- enr_res |>
      dplyr::inner_join(
        go_group_genes_df,
        by = "go_group"
      ) |>
      dplyr::mutate(go_group = gsub("_"," ", go_group))
  } else{
    enr_res <- enr_res |>
      dplyr::mutate(go_group = gsub("_"," ", go_group))
  }

  if ("category" %in% colnames(enr_res)){
    enr_res |>
      dplyr::mutate(go_group = paste(go_group,gene_col,sep="\n")) |>
      dplyr::group_by(category,exp_name,go_group) |>
      dplyr::arrange(p_adjust) |>
      dplyr::slice_head(n=n_per_group) |>
      dplyr::ungroup() |>
      ggplot(aes(
        x = forcats::fct_reorder(exp_name, -n_with_sel),
        y = forcats::fct_reorder(term_name, -log10(p_adjust)),
        size = n_with_sel,
        color = -log10(p_adjust)
      )) +
      geom_point() +
      ggpubr::theme_pubr(legend = "right",
                         base_size = 10) +
      facet_grid(go_group~category,
                 space ="free",
                 scales = "free") +
      ggpubr::rotate_x_text(angle=45)+
      scale_color_gradient(low="thistle1",high="#8E0152") +
      theme(strip.text.x = element_text(
        face = "bold.italic",
        angle=90),
        strip.text.y = element_text(
          face="bold.italic",
          angle=0))+
      labs(
        x="",
        y="",
        title="",
        size="Number of Genes",
        color=expression(-Log[10](P)),
      )
  } else{
    enr_res |>
      dplyr::mutate(go_group = paste(go_group,gene_col,sep="\n")) |>
      dplyr::group_by(exp_name,go_group) |>
      dplyr::arrange(p_adjust) |>
      dplyr::slice_head(n=n_per_group) |>
      dplyr::ungroup() |>
      ggplot(aes(
        x = forcats::fct_reorder(exp_name, -n_with_sel),
        y = forcats::fct_reorder(term_name, -log10(p_adjust)),
        size = n_with_sel,
        color = -log10(p_adjust)
      )) +
      geom_point() +
      ggpubr::theme_pubr(legend = "right",
                         base_size = 10) +
      facet_grid(go_group~., space ="free",scales = "free") +
      ggpubr::rotate_x_text(angle=45)+
      scale_color_gradient(low="thistle1",high="#8E0152") +
      theme(strip.text.x = element_text(
        face = "bold.italic",
        angle=90),
        strip.text.y = element_text(
          face="bold.italic",
          angle=0))+
      labs(
        x="",
        y="",
        title="",
        size="Number of Genes",
        color=expression(-Log[10](P)),
      )
  }
}

# --- Enrichment Heatmap -----
#' @noRd
#'
#' @title Internal: Heatmap of Enrichment Term × Gene LogFC Matrix
#'
#' @description
#' Generates a tile-based heatmap of enrichment terms by associated genes,
#'  optionally colored by log2 fold change.
#'
#' @param enr_res Enrichment results with exploded gene-term relationships.
#' @param expomicset The original `MultiAssayExperiment` object.
#' @param feature_type Feature source (e.g., "degs", "degs_robust").
#' @param score_metric Column name for stability filtering (for `degs_robust`).
#' @param score_thresh Numeric threshold for stability score filtering.
#' @param go_groups GO groups to include in the heatmap.
#' @param heatmap_fill Logical; whether to fill tiles by logFC.
#' @param feature_col Column for matching features.
#' @param logfc_col Column for log2 fold change values.
#' @param logfc_thresh Threshold for absolute logFC.
#' @param pval_col Column for p-values.
#' @param pval_thresh Significance threshold.
#'
#' @return A `ggplot` heatmap object.
#'
#' @importFrom ggplot2 ggplot aes theme_bw facet_grid theme
#' element_text labs geom_tile scale_fill_gradient2
#' @importFrom dplyr filter distinct count mutate left_join inner_join select
#' @importFrom purrr pluck
#' @importFrom tidyr separate_rows
#' @importFrom forcats fct_reorder fct_reorder2
#'
#' @keywords internal
.plot_heatmap_enrichment <- function(
    enr_res = enr_res,
    expomicset = expomicset,
    feature_type = feature_type,
    score_metric = score_metric,
    score_thresh = score_thresh,
    go_groups = go_groups,
    heatmap_fill = heatmap_fill,
    feature_col = feature_col,
    logfc_col = logfc_col,
    logfc_thresh = logfc_thresh,
    pval_col = pval_col,
    pval_thresh = pval_thresh
){
  enr_res <- enr_res |>
    tidyr::separate_rows("ids",sep=", ")

  if(is.null(go_groups)){
    stop("Please specify a `go_group`")
  } else {
    enr_res <- enr_res |>
      dplyr::filter(go_group %in% go_groups)
  }

  shared_genes <- enr_res |>
    tidyr::separate_rows(ids, sep = ", ") |>
    dplyr::distinct(exp_name, ids) |>
    dplyr::count(ids, name = "n_exp") |>
    dplyr::mutate(shared = n_exp > 1)

  if (heatmap_fill){

    if (feature_type == "degs_robust"){

      if (is.null (score_thresh)){
        score_thresh <- expomicset@metadata |>
          purrr::pluck("differential_analysis") |>
          purrr::pluck("sensitivity_analysis") |>
          purrr::pluck("score_thresh")

      } else {
        score_thresh <- score_thresh
      }

      # grab feature information from sensitivity analysis
      f_stable_df <- expomicset@metadata |>
        purrr::pluck("differential_analysis") |>
        purrr::pluck("sensitivity_analysis") |>
        purrr::pluck("feature_stability") |>
        dplyr::filter(!!dplyr::sym(score_metric) > score_thresh) |>
        dplyr::select(c(feature,exp_name))

      enr_res <- enr_res |>
        inner_join(
          expomicset@metadata |>
            purrr::pluck("differential_analysis") |>
            purrr::pluck("differential_abundance") |>
            dplyr::select(c(exp_name,
                            feature,
                            feature_col,
                            logfc_col)),
          by=c("exp_name"="exp_name",
               "ids" = feature_col),
          relationship = "many-to-many"
        ) |>
        inner_join(
          f_stable_df,
          by = c(
            "exp_name"="exp_name",
            "feature" = "feature"
          ))

    } else if (feature_type == "degs") {
      enr_res <- enr_res |>
        inner_join(
          expomicset@metadata |>
            purrr::pluck("differential_analysis") |>
            purrr::pluck("differential_abundance") |>
            dplyr::filter(!!dplyr::sym(pval_col) < pval_thresh,
                          abs(!!dplyr::sym(logfc_col)) > logfc_thresh) |>
            dplyr::select(c(exp_name,
                            feature,
                            feature_col,
                            logfc_col)),
          by=c("exp_name"="exp_name",
               "ids" = feature_col),
          relationship = "many-to-many"
        )
    } else{
      enr_res <- enr_res |>
        inner_join(
          expomicset@metadata |>
            purrr::pluck("differential_analysis") |>
            purrr::pluck("differential_abundance") |>
            dplyr::select(c(exp_name,
                            feature,
                            feature_col,
                            logfc_col)),
          by=c("exp_name"="exp_name",
               "ids" = feature_col),
          relationship = "many-to-many"
        )
    }
  }

  # Join shared status
  enr_res <- enr_res |>
    dplyr::left_join(shared_genes, by = "ids") |>
    dplyr::mutate(
      ids = ifelse(shared, paste0("*",ids), ids)
    )

  go_groups_clean <- gsub("_"," ",go_groups)

  if ("category" %in% colnames(enr_res)){
    p <- enr_res |>
      ggplot(aes(
        x = forcats::fct_reorder2(ids, term_name,category),
        y = forcats::fct_reorder2(term_name, exp_name,category)
      )) +
      theme_bw() +
      facet_grid(category~exp_name,
                 scales = "free",
                 space = "free") +
      theme(
        axis.text.x = element_text(angle = 90, face = "italic"),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(face = "bold.italic"),
        strip.text.y = element_text(face = "bold.italic"),
        plot.caption = element_text(face = "italic"),
        plot.title = element_text(face = "bold.italic")
      ) +
      labs(
        caption = "* Appears in more than one omic",
        x = "",
        y = "",
        title = paste("Enrichment Term Heatmap for ",
                      paste(go_groups_clean, collapse = ","))
      )
  } else {
    p <- enr_res |>
      ggplot(aes(
        x = forcats::fct_reorder(ids, term_name),
        y = forcats::fct_reorder(term_name, exp_name)
      )) +
      theme_bw() +
      facet_grid(~exp_name,
                 scales = "free",
                 space = "free") +
      theme(
        axis.text.x = element_text(angle = 90, face = "italic"),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(face = "bold.italic"),
        plot.caption = element_text(face = "italic"),
        plot.title = element_text(face = "bold.italic")
      ) +
      labs(
        caption = "* Appears in more than one omic",
        x = "",
        y = "",
        title = paste("Enrichment Term Heatmap for ",
                      paste(go_groups_clean, collapse = ","))
      )
  }

  if (heatmap_fill) {
    p <- p +
      geom_tile(aes(fill = !!sym(logfc_col)))+
      scale_fill_gradient2(
        name = expression("Log"[2]*"FC"),
        low = "blue4",
        mid = "white",
        high = "red4",
        midpoint = 0
      )
  } else {
    p <- p +
      geom_tile(fill = "grey30")
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
#' @importFrom tidygraph as_tbl_graph activate centrality_degree
#' @importFrom ggraph create_layout ggraph geom_edge_link scale_edge_width
#' geom_node_arc_bar
#' @importFrom ggplot2 facet_wrap theme element_text labs scale_fill_manual
#'  theme_void
#' @importFrom ggrepel geom_label_repel
#' @importFrom dplyr filter mutate select distinct group_by ungroup
#' left_join top_n pull
#' @importFrom purrr map
#' @importFrom tidyr separate_rows
#' @importFrom stringr str_trim
#' @importFrom tibble as_tibble
#' @importFrom RColorBrewer brewer.pal
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
  #require(igraph)
  #require(tidygraph)
  #require(ggraph)

  # Get enrichment table
  enr <- enr_res

  if(!is.null(go_groups)){
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
    dplyr::select(source = source_std,
                  target = target_std,
                  overlap) |>
    dplyr::filter(overlap>overlap_thresh) |>
    dplyr::distinct()

  # Ensure unique node list by term_name
  node_tbl <- enr |>
    dplyr::select(term_name) |>
    dplyr::distinct() |>
    dplyr::mutate(name = term_name)  # igraph expects "name"

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
    stop("Join failed: 'exp_name' not found in arc_data.
         Check that term_name is consistent.")
  }

  # Default pie colors
  if (is.null(pie_colors)) {
    omics_types <- unique(arc_data$exp_name)
    n_colors <- max(3, length(omics_types))
    pie_colors <- setNames(
      RColorBrewer::brewer.pal(
        n = n_colors,
        "Set2")[seq_along(omics_types)], omics_types)
  }

  # Final plot
  p <- ggraph(layout) +
    geom_edge_link(aes(width = overlap),
                   alpha = edge_alpha,
                   color = "gray50") +
    scale_edge_width(range = c(0.2, 2),guide = "none")+
    geom_node_arc_bar(
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
      stop(glue::glue(
        "Could not join `net_facet_by = '{net_facet_by}'` column."))
    }

    p <- p +
      facet_wrap(as.formula(paste("~", net_facet_by)))+
      theme(strip.text.x = element_text(
        face="bold.italic",
        size = 12))
  }

  if (!is.null(label_top_n)){

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
        force = 2,        # stronger repel
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
#' @title Internal: Gene–Term Network (Cnetplot)
#'
#' @description
#' Builds a bipartite graph between genes and enriched terms,
#' with logFC values for genes and pie slices for multi-omic terms.
#'
#' @param enr_res Enrichment results.
#' @param expomicset A `MultiAssayExperiment` object.
#' @param feature_type Source of features ("degs", "omics", etc.).
#' @param layout_algo Graph layout algorithm. Default is `"fr"`.
#' @param go_groups Optional subset of GO groups.
#' @param max_terms Maximum number of terms to include.
#' @param edge_alpha Transparency of edges.
#' @param node_size Base size for graph nodes.
#' @param label_size Size of text labels.
#' @param term_node_correction Scaling factor for term nodes.
#' @param gene_node_correction Scaling factor for gene nodes.
#' @param feature_col Feature column used for matching.
#' @param logfc_col Column name with log fold change values.
#'
#' @return A `ggraph` plot object.
#'
#' @importFrom ggraph ggraph create_layout geom_edge_link
#' geom_node_arc_bar geom_node_point geom_node_label
#' @importFrom igraph graph_from_data_frame
#' @importFrom tidygraph as_tbl_graph
#' @importFrom ggplot2 scale_color_gradient2 theme_void labs element_text
#' @importFrom dplyr filter mutate select distinct left_join
#' count slice_head pull bind_rows group_by ungroup add_count if_else
#' @importFrom tidyr separate_rows
#' @importFrom tibble as_tibble
#' @importFrom purrr pluck
#' @importFrom ggnewscale new_scale_color
#'
#' @keywords internal
.plot_cnet_enrichment <- function(
    enr_res,
    expomicset,
    feature_type,
    layout_algo = "fr",
    go_groups = NULL,
    max_terms = 30,
    edge_alpha = 0.5,
    node_size = 1,
    label_size = 3,
    term_node_correction=0.2,
    gene_node_correction=3,
    feature_col = "feature_clean",
    logfc_col = "logFC"
) {
  #require(ggraph)
  #require(igraph)
  #require(tidygraph)
  #require(ggnewscale)

  # Explode to gene–term–exp rows
  long_enr <- enr_res |>
    tidyr::separate_rows(ids,
                         sep = ",\\s*") |>
    dplyr::select(term_name,
                  go_group,
                  ids,
                  exp_name) |>
    dplyr::distinct()

  # Filter GO groups or select top terms
  if (!is.null(go_groups)) {
    long_enr <- long_enr |>
      dplyr::filter(go_group %in% go_groups)
  } else {
    top_terms <- long_enr |>
      dplyr::count(term_name, sort = TRUE) |>
      dplyr::slice_head(n = max_terms) |>
      dplyr::pull(term_name)
    long_enr <- long_enr |>
      dplyr::filter(term_name %in% top_terms)
  }

  # Join logFC from DA results
  da_df <- expomicset@metadata |>
    purrr::pluck("differential_analysis",
                 "differential_abundance") |>
    dplyr::select(exp_name,
                  feature,
                  !!sym(feature_col),
                  !!sym(logfc_col)) |>
    dplyr::distinct()

  long_enr <- long_enr |>
    dplyr::left_join(
      da_df,
      by = c("exp_name",
             "ids" = feature_col),
      relationship = "many-to-many"
    )

  # Assign unique node names
  long_enr <- long_enr |>
    dplyr::mutate(
      gene_node = paste0("gene:", ids, "::", exp_name),
      term_node = paste0("term:", term_name)
    )

  # Build edges
  edges <- long_enr |>
    dplyr::select(from = gene_node, to = term_node)

  # Node metadata
  gene_meta <- long_enr |>
    dplyr::select(name = gene_node,
                  !!sym(logfc_col)) |>
    dplyr::distinct() |>
    dplyr::mutate(type = "gene")

  term_meta <- long_enr |>
    dplyr::select(name = term_node,
                  exp_name) |>
    dplyr::distinct() |>
    dplyr::mutate(type = "term")

  nodes <- dplyr::bind_rows(term_meta,
                            gene_meta) |>
    dplyr::distinct(name, .keep_all = TRUE)

  # Graph layout
  g <- igraph::graph_from_data_frame(
    edges,
    vertices = nodes,
    directed = FALSE)

  layout <- ggraph::create_layout(g,
                                  layout = layout_algo)

  # Pie slices for term nodes
  term_arc_data <- long_enr |>
    dplyr::mutate(name = paste0("term:", term_name)) |>
    dplyr::select(name, exp_name) |>
    dplyr::distinct() |>
    dplyr::group_by(name) |>
    dplyr::mutate(
      n = n(),
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
      label = if_else(n_exps > 1, paste0(ids, " (", exp_name, ")"), ids),
      name = paste0("gene:", ids, "::", exp_name)
    ) |>
    dplyr::select(name, label)

  # Merge labels into layout
  layout <- layout |>
    dplyr::left_join(gene_labels,
                     by = "name") |>
    dplyr::mutate(
      label = ifelse(is.na(label),
                     gsub("^.*?:", "", name),
                     label))

  # Build plot
  p <- ggraph(layout) +
    geom_edge_link(alpha = edge_alpha, color = "grey70")

  # Term nodes as pie chart by exp_name
  p <- p +
    geom_node_arc_bar(
      data = term_arc_data,
      aes(
        x0 = x, y0 = y, r0 = 0, r = node_size * term_node_correction,
        start = start, end = end, fill = exp_name
      ),
      color = "transparent"
    ) +
    scale_fill_tidy_exp(name = "Experiment")

  # Reset scale before gene coloring
  p <- p +
    ggnewscale::new_scale_color()

  # Gene nodes colored by logFC
  p <- p +
    geom_node_point(
      data = layout |>
        dplyr::filter(type == "gene"),
      aes(x = x,
          y = y,
          color = !!sym(logfc_col)),
      size = node_size * gene_node_correction
    ) +
    scale_color_gradient2(
      name = expression(Log[2]*"FC"),
      low = "blue4",
      mid = "white",
      high = "red4",
      midpoint = 0,
      na.value = "grey60"
    )

  # Add labels
  p <- p +
    geom_node_label(
      data = layout,
      aes(
        x = x, y = y,
        label = label,
        fontface = ifelse(type == "gene", "italic", "bold")
      ),
      size = label_size,
      repel = TRUE
    ) +
    theme_void()

  p <-  p +
    labs(
      title = "Gene–Term Network"
    )+
    theme(plot.title = element_text(face="bold.italic"))

  return(p)
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
#' @importFrom ggplot2 ggplot aes geom_col coord_flip labs theme_minimal
#' theme element_blank element_text element_line element_rect
#' element_text scale_y_continuous scale_fill_gradient
#'  scale_size_continuous scale_fill_manual
#' @importFrom dplyr mutate distinct count filter left_join
#' arrange pull select group_by ungroup
#' @importFrom tidyr separate_rows
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom purrr map
#' @importFrom ggridges geom_density_ridges theme_ridges
#' @importFrom ggvenn ggvenn
#'
#' @keywords internal
.plot_summary_enrichment <- function(
    enr_res
) {
  #require(ggplot2)

  # Clean up enr_res so that go_groups have no underscore
  enr_res <- enr_res |>
    dplyr::mutate(go_group = gsub("_", " ", go_group))

  # --- Barplot: # of enriched terms per exp_name
  bar_exp <- enr_res |>
    dplyr::distinct(term_name,
                    exp_name) |>
    dplyr::count(exp_name) |>
    ggplot(aes(x = reorder(exp_name, n),
               y = n,
               fill = exp_name)) +
    geom_col() +
    coord_flip() +
    scale_fill_tidy_exp() +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    labs(y = "No. Enriched Terms",
         x = "") +
    theme_minimal() +
    theme(legend.position = "none")

  # --- Barplot: # of enriched terms per category
  if ("category" %in% colnames(enr_res)) {
    bar_cat <- enr_res |>
      dplyr::distinct(term_name, category) |>
      dplyr::count(category) |>
      ggplot(aes(x = reorder(category, n),
                 y = n,
                 fill = category)) +
      geom_col() +
      coord_flip() +
      scale_fill_tidy_exp(rev = TRUE) +
      labs(y = "No. Enriched Terms",
           x = "Category") +
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
        labs(x = "No. Categories",
             y = NULL,
             fill = "Category") +
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
      labs(x = expression("-Log"[10]*"P"), y = NULL, fill = "No. Genes") +
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
    split(~ exp_name) |>
    purrr::map(~ .x$term_name)

  venn <- ggvenn::ggvenn(
    venn_data,
    fill_color = tidy_exp_pal,
    stroke_color = "transparent",
    set_name_size = 4,
    show_percentage = TRUE
  )+
    theme(
      plot.margin = margin(0, 0, 0, 0),
      plot.background = element_blank(),
      panel.background = element_blank()

    )

  # --- Layout assembly
  layout_left <- ridge
  layout_right <- (bar_exp / venn) +
    patchwork::plot_layout(heights = c(1, 2.5))

  final_plot <- (layout_left | layout_right)+
    patchwork::plot_annotation(title = "Enrichment Summary") &
    theme(plot.title = element_text(face="bold.italic"))

  return(final_plot)
}



