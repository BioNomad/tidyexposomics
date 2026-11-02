#' Generate Example Data for Testing
#'
#' This helper function generates a reproducible dummy dataset
#' containing exposures, mRNA data, and proteomics data. It can
#' optionally return the data as a \code{MultiAssayExperiment} using
#' \code{\link{create_exposomicset}}.
#'
#' @param n_samples Integer. Number of samples to simulate (default: 12).
#' @param n_proteins Integer. Number of proteins to simulate (default: 80).
#' @param use_batch Logical. If \code{TRUE},
#' include a "batch" variable in the exposure data (default: FALSE).
#' @param return_mae Logical. If \code{TRUE},
#' return a \code{MultiAssayExperiment} created using
#' \code{create_exposomicset()} (default: FALSE).
#'
#' @return Either:
#' \itemize{
#'   \item A named list containing \code{codebook},
#'    \code{exposure}, \code{omics}, and \code{row_data},
#'    if \code{return_mae = FALSE}.
#'   \item A \code{MultiAssayExperiment}, if \code{return_mae = TRUE}.
#' }
#'
#' @examples
#' # Return as a list
#' dummy <- make_example_data()
#'
#' # Return as a MultiAssayExperiment
#' mae <- make_example_data(return_mae = TRUE)
#'
#' @export
make_example_data <- function(
  n_samples = 12,
  n_proteins = 80,
  use_batch = FALSE,
  return_mae = FALSE
) {
    sample_ids <- paste0("S", seq_len(n_samples))

    # Gene sets with some known function
    xeno_genes <- c(
        "AADAC", "BAK1", "SLCO1B1", "GATA3", "POR", "HTR2B",
        "CYP26A1", "SLC28A3", "EMX2", "KCNH2",
        "RORC", "RBM22", "CCND1", "TLR3", "LCN2", "CRKL",
        "NUDT15", "XPO1", "ROGDI", "FZD1",
        "BCL2", "NTRK1", "PTPRM", "MIR34B",
        "FUT1", "CDH1", "ABCD3", "FMO4", "SLCO1B3", "CYP2U1"
    )

    antigen_genes <- c(
        "ERAP1", "PSMB8", "IFI30", "CD1D", "IGHE", "MARCHF8",
        "THBS1", "FCER2", "FCGR2B", "LILRB2",
        "CTSH", "HLA-DRB5", "LGMN", "CD209",
        "HLA-DRB1", "HLA-E", "TAPBP", "CST7", "CD1E", "RAB10",
        "HLA-G", "MARCHF1", "HLA-DQB2", "HLA-DRA", "TRAF6",
        "PIKFYVE", "CCL21", "HLA-DPB1", "RAB35", "TAPBPL"
    )

    # Randomly sampled gene sets from msig
    other_genes <- c(
        "MYCT1", "NCDN", "MGAT4C", "ATF1", "REEP6", "NEURL4", "HOATZ",
        "RBM39", "HLA-DRB5", "EFNA1",
        "MRPS18A", "MTO1", "MIRLET7A1", "MCTS2", "ZMAT3", "CHRM1",
        "TRAV22", "PBLD", "NOTCH2", "RHBDL1",
        "KIF27", "CXCL10", "ASB8", "CDKN1C", "ZPLD1", "MTHFD2", "ANKIB1",
        "POLD2", "TOM1L2", "SYVN1",
        "PPP3CA", "PRAMEF27", "ARHGEF3", "TBC1D24", "EMG1", "CCL21",
        "TUSC3", "SNORD111", "CEBPA", "TMEM63C",
        "SLC7A8", "ZNF738", "NEU3", "FSCN2", "BUB1", "MYEF2", "GMIP",
        "BAIAP2L1", "SEPSECS", "EPB41L1",
        "RIMS3", "IFT70A", "GLIS2", "IGLV3-12", "CGRRF1", "MESP2",
        "TCEA2", "TIMP2", "HECTD2", "CRTAC1",
        "C1R", "BRMS1L", "IZUMO1", "FDX1", "RPL14", "SLC35B2", "KRT16",
        "DPH2", "OCSTAMP", "NME5"
    )

    gene_symbols <- unique(c(antigen_genes, xeno_genes, other_genes))
    n_genes <- length(gene_symbols)
    mrna_ids <- paste0("feat_", seq_len(n_genes))

    # Simulate exposures
    exposure <- data.frame(
        age = sample(30:70, n_samples, replace = TRUE),
        sex = sample(c("M", "F"), n_samples, replace = TRUE),
        bmi = round(rnorm(n_samples, 25, 3), 1),
        smoker = sample(c("yes", "no"), n_samples, replace = TRUE),
        alcohol = sample(c("none", "low", "moderate", "high"),
            n_samples,
            replace = TRUE
        ),
        exposure_pm25 = runif(n_samples, 5, 30),
        exposure_no2 = runif(n_samples, 10, 60),
        row.names = sample_ids
    )

    if (use_batch) {
        exposure$batch <- sample(c("A", "B"), n_samples, replace = TRUE)
    }

    # Codebook with categories
    codebook <- data.frame(
        variable = colnames(exposure),
        description = c(
            "Age of subject", "Sex at birth", "Body mass index",
            "Current smoking status", "Alcohol intake",
            "Air pollution: PM2.5", "Air pollution: NO2",
            if (use_batch) "Batch label" else NULL
        ),
        category = c(
            "Demographics", "Demographics", "Demographics",
            "Lifestyle", "Lifestyle",
            "Environment", "Environment",
            if (use_batch) "Technical" else NULL
        ),
        stringsAsFactors = FALSE
    )

    # Simulate latent structure
    lv1 <- scale(rnorm(n_samples)) # antigen signal
    lv2 <- scale(rnorm(n_samples)) # xeno signal

    # Simulate mRNA matrix - continuous, positive values
    mrna_mat <- matrix(
        rlnorm(n_genes * n_samples,
            meanlog = 7.5,
            sdlog = 1
        ),
        nrow = n_genes,
        dimnames = list(
            mrna_ids,
            sample_ids
        )
    )

    # Add latent structure and exposure signal
    mrna_mat[mrna_ids[gene_symbols %in% antigen_genes], ] <- mrna_mat[
        mrna_ids[gene_symbols %in% antigen_genes],
    ] +
        matrix(
            rep(lv1,
                each = sum(gene_symbols %in% antigen_genes)
            ),
            nrow = sum(gene_symbols %in% antigen_genes)
        )

    mrna_mat[mrna_ids[gene_symbols %in% xeno_genes], ] <- mrna_mat[
        mrna_ids[gene_symbols %in% xeno_genes],
    ] +
        matrix(
            rep(lv2,
                each = sum(gene_symbols %in% xeno_genes)
            ),
            nrow = sum(gene_symbols %in% xeno_genes)
        )

    mrna_mat[mrna_ids[gene_symbols %in% antigen_genes], ] <- mrna_mat[
        mrna_ids[gene_symbols %in% antigen_genes],
    ] +
        matrix(
            rep(scale(exposure$exposure_pm25),
                each = sum(gene_symbols %in% antigen_genes)
            ),
            nrow = sum(gene_symbols %in% antigen_genes)
        )

    # Create rowData for mRNA
    row_data_mrna <- S4Vectors::DataFrame(
        symbol = gene_symbols,
        row.names = mrna_ids
    )

    # Simulate proteomics with real gene names
    matched_genes <- sample(gene_symbols, n_proteins, replace = TRUE)
    prot_ids <- paste0("prot_", seq_len(n_proteins))

    prot_mat <- mrna_mat[mrna_ids[match(matched_genes, gene_symbols)], ] +
        matrix(
            rlnorm(n_proteins * n_samples,
                meanlog = 7.5,
                sdlog = 0.5
            ),
            nrow = n_proteins
        )

    rownames(prot_mat) <- prot_ids
    colnames(prot_mat) <- sample_ids

    # Create rowData for proteomics
    row_data_prot <- S4Vectors::DataFrame(
        symbol = matched_genes,
        row.names = prot_ids
    )

    # Combine everything
    row_data <- list(
        mRNA = row_data_mrna,
        proteomics = row_data_prot
    )

    omics <- list(
        mRNA = mrna_mat,
        proteomics = prot_mat
    )

    if (return_mae) {
        res_list <- list(
            codebook = codebook,
            exposure = exposure,
            omics = omics,
            row_data = row_data
        )
        res <- create_exposomicset(
            codebook = res_list$codebook,
            exposure = res_list$exposure,
            omics = res_list$omics,
            row_data = res_list$row_data
        )
    } else {
        res <- list(
            codebook = codebook,
            exposure = exposure,
            omics = omics,
            row_data = row_data
        )
    }

    return(res)
}


#' Make an example single-cell dataset with multiple exposures
#'
#' Creates a sample `SingleCellExperiment` object with
#' sample-level exposures (two continuous, one binary, plus age), cell types,
#' and negative binomial counts. A subset of genes truly responds to each
#' exposure, and some exposures act only in specific cell types. The true
#' effect sizes are stored in `rowData()` so results can be benchmarked.
#'
#' Exposures are constant within a sample, so the effective sample size is the
#' number of samples, not the number of cells.
#'
#' @param n_samples Number of samples (donors).
#' @param cells_per_sample Number of cells simulated per sample.
#' @param n_genes Number of genes.
#' @param n_celltypes Number of cell types (named `CT1`, `CT2`, ...).
#' @param prop_responsive Proportion of genes with a true effect for each
#'   exposure.
#' @param effect_sd Standard deviation of the true (log-scale) exposure effects.
#' @param prop_celltype_responsive Probability that a cell type responds to a
#'   given exposure. `CT1` always responds to every exposure.
#' @param size Negative binomial size (dispersion) parameter. Smaller values
#'   give more overdispersion.
#'
#' @return A `SingleCellExperiment` with:
#' \describe{
#'   \item{`assays$counts`}{Genes x cells count matrix.}
#'   \item{`colData`}{`sample`, `pm25`, `bpa`, `smoker`, `age`, `celltype`.}
#'   \item{`rowData`}{`beta_pm25`, `beta_bpa`, `beta_smoker`: true effects,
#'     where nonzero means a truly responsive gene.}
#'   \item{`metadata$responsive`}{Cell type x exposure matrix of which cell
#'     types respond to which exposures.}
#' }
#'
#' @examples
#' sce <- make_example_sc()
#'
#' @export
make_example_sc <- function(n_samples = 12,
                            cells_per_sample = 200,
                            n_genes = 1000,
                            n_celltypes = 4,
                            prop_responsive = 0.08,
                            effect_sd = 0.6,
                            prop_celltype_responsive = 0.6,
                            size = 2) {
    n_cells <- n_samples * cells_per_sample

    # sample-level exposures, repeated across each sample's cells
    sample_df <- data.frame(
        sample = sprintf("S%02d", seq_len(n_samples)),
        pm25   = stats::rlnorm(n_samples, 2.5, 0.4),
        bpa    = stats::rlnorm(n_samples, 0, 0.6),
        smoker = factor(stats::rbinom(n_samples, 1, 0.4), levels = 0:1, labels = c("no", "yes")),
        age    = round(stats::runif(n_samples, 25, 70))
    )

    cell_df <- sample_df[rep(seq_len(n_samples), each = cells_per_sample), ] |>
        transform(
            celltype = factor(
                sample(paste0("CT", seq_len(n_celltypes)), n_cells, replace = TRUE)
            )
        )
    rownames(cell_df) <- sprintf("cell%05d", seq_len(n_cells))

    ct_idx <- as.integer(cell_df$celltype)
    sample_idx <- rep(seq_len(n_samples), each = cells_per_sample)

    # standardized exposure matrix: cells x exposures
    Z <- cbind(
        pm25   = scale(cell_df$pm25)[, 1],
        bpa    = scale(cell_df$bpa)[, 1],
        smoker = as.numeric(cell_df$smoker == "yes")
    )

    # true exposure effects: genes x exposures
    beta <- matrix(
        stats::rbinom(n_genes * 3, 1, prop_responsive) * stats::rnorm(n_genes * 3, 0, effect_sd),
        nrow = n_genes,
        dimnames = list(NULL, paste0("beta_", colnames(Z)))
    )

    # which cell types respond to which exposure (CT1 always responds)
    responsive <- matrix(
        stats::rbinom(n_celltypes * 3, 1, prop_celltype_responsive),
        nrow = n_celltypes,
        dimnames = list(paste0("CT", seq_len(n_celltypes)), colnames(Z))
    )
    responsive[1, ] <- 1

    # log-mean: gene baseline + cell type + exposure (cell-type specific) + sample + cell
    log_mu <- stats::rnorm(n_genes, -1, 1.2) +
        matrix(
            stats::rbinom(n_genes * n_celltypes, 1, 0.15) * stats::rnorm(n_genes * n_celltypes, 0, 1.5),
            nrow = n_genes
        )[, ct_idx] +
        beta %*% t(Z * responsive[ct_idx, ]) +
        matrix(stats::rnorm(n_genes * n_samples, 0, 0.2), nrow = n_genes)[, sample_idx] +
        rep(stats::rnorm(n_cells, 0, 0.3), each = n_genes)

    genes <- sprintf("gene%04d", seq_len(n_genes))

    SingleCellExperiment::SingleCellExperiment(
        assays = list(
            counts = matrix(
                stats::rnbinom(n_genes * n_cells, mu = exp(log_mu), size = size),
                nrow = n_genes,
                dimnames = list(genes, rownames(cell_df))
            )
        ),
        colData = cell_df,
        rowData = data.frame(beta, row.names = genes),
        metadata = list(responsive = responsive)
    )
}


#' Make an example codebook for the simulated single-cell data
#'
#' Returns a data dictionary describing the `colData()` variables produced by
#' `make_example_sc()`. Each row is one variable, annotated with a label, its
#' role in an analysis, its exposure domain, the level at which it varies, and
#' its unit.
#'
#' @return A data frame with one row per variable and columns:
#' \describe{
#'   \item{`variable`}{Column name in `colData()`.}
#'   \item{`label`}{Human-readable description.}
#'   \item{`role`}{One of `"id"`, `"exposure"`, `"covariate"`, or
#'     `"cell_annotation"`.}
#'   \item{`domain`}{Category of the variable (e.g. `"air pollution"`).}
#'   \item{`level`}{Whether the variable varies by `"sample"` or by `"cell"`.}
#'   \item{`unit`}{Unit of measurement, or `NA` if not applicable. Units are
#'     illustrative, since the simulator does not assign real ones.}
#' }
#'
#'
#' @examples
#' make_example_sce_codebook()
#'
#' make_example_sce_codebook() |>
#'     dplyr::filter(role == "exposure")
#'
#' @export
make_example_sce_codebook <- function() {
    data.frame(
        variable = c("sample", "pm25", "bpa", "smoker", "age", "celltype"),
        label = c(
            "Sample ID",
            "Fine particulate matter (PM2.5)",
            "Bisphenol A (BPA)",
            "Current smoking status",
            "Age at collection",
            "Cell type annotation"
        ),
        role = c("id", "exposure", "exposure", "exposure", "covariate", "cell_annotation"),
        domain = c("identifier", "air pollution", "chemical", "lifestyle", "demographic", "cell type"),
        level = c("sample", "sample", "sample", "sample", "sample", "cell"),
        unit = c(NA, "ug/m3", "ng/mL", NA, "years", NA)
    )
}
