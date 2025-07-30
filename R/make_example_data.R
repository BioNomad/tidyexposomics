#' Generate Example Data for Testing
#'
#' This helper function generates a reproducible dummy dataset
#' containing exposures, mRNA data, and proteomics data. It can
#' optionally return the data as a \code{MultiAssayExperiment} using
#' \code{\link{create_expomicset}}.
#'
#' @param n_samples Integer. Number of samples to simulate (default: 12).
#' @param n_proteins Integer. Number of proteins to simulate (default: 80).
#' @param use_batch Logical. If \code{TRUE},
#' include a "batch" variable in the exposure data (default: FALSE).
#' @param return_mae Logical. If \code{TRUE},
#' return a \code{MultiAssayExperiment} created using
#' \code{create_expomicset()} (default: FALSE).
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
  xeno_genes <- c('AADAC','BAK1','SLCO1B1','GATA3','POR','HTR2B',
                  'CYP26A1','SLC28A3','EMX2','KCNH2',
                  'RORC','RBM22','CCND1','TLR3','LCN2','CRKL',
                  'NUDT15','XPO1','ROGDI','FZD1',
                  'BCL2','NTRK1','PTPRM','MIR34B',
                  'FUT1','CDH1','ABCD3','FMO4','SLCO1B3','CYP2U1')

  antigen_genes <- c('ERAP1','PSMB8','IFI30','CD1D','IGHE','MARCHF8',
                     'THBS1','FCER2','FCGR2B','LILRB2',
                     'CTSH','HLA-DRB5','LGMN','CD209',
                     'HLA-DRB1','HLA-E','TAPBP','CST7','CD1E','RAB10',
                     'HLA-G','MARCHF1','HLA-DQB2','HLA-DRA','TRAF6',
                     'PIKFYVE','CCL21','HLA-DPB1','RAB35','TAPBPL')

  # Randomly sampled gene sets from msig
  other_genes <- c('MYCT1','NCDN','MGAT4C','ATF1','REEP6','NEURL4','HOATZ',
                   'RBM39','HLA-DRB5','EFNA1',
                   'MRPS18A','MTO1','MIRLET7A1','MCTS2','ZMAT3','CHRM1',
                   'TRAV22','PBLD','NOTCH2','RHBDL1',
                   'KIF27','CXCL10','ASB8','CDKN1C','ZPLD1','MTHFD2','ANKIB1',
                   'POLD2','TOM1L2','SYVN1',
                   'PPP3CA','PRAMEF27','ARHGEF3','TBC1D24','EMG1','CCL21',
                   'TUSC3','SNORD111','CEBPA','TMEM63C',
                   'SLC7A8','ZNF738','NEU3','FSCN2','BUB1','MYEF2','GMIP',
                   'BAIAP2L1','SEPSECS','EPB41L1',
                   'RIMS3','IFT70A','GLIS2','IGLV3-12','CGRRF1','MESP2',
                   'TCEA2','TIMP2','HECTD2','CRTAC1',
                   'C1R','BRMS1L','IZUMO1','FDX1','RPL14','SLC35B2','KRT16',
                   'DPH2','OCSTAMP','NME5')

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
                     replace = TRUE),
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
  lv1 <- scale(rnorm(n_samples))  # antigen signal
  lv2 <- scale(rnorm(n_samples))  # xeno signal

  # Simulate mRNA matrix - continuous, positive values
  mrna_mat <- matrix(rlnorm(n_genes * n_samples,
                            meanlog = 7.5,
                            sdlog = 1),
                     nrow = n_genes,
                     dimnames = list(mrna_ids,
                                     sample_ids))

  # Add latent structure and exposure signal
  mrna_mat[mrna_ids[gene_symbols %in% antigen_genes], ] <- mrna_mat[
    mrna_ids[gene_symbols %in% antigen_genes], ] +
    matrix(rep(lv1,
               each = sum(gene_symbols %in% antigen_genes)),
           nrow = sum(gene_symbols %in% antigen_genes))

  mrna_mat[mrna_ids[gene_symbols %in% xeno_genes], ] <- mrna_mat[
    mrna_ids[gene_symbols %in% xeno_genes], ] +
    matrix(rep(lv2,
               each = sum(gene_symbols %in% xeno_genes)),
           nrow = sum(gene_symbols %in% xeno_genes))

  mrna_mat[mrna_ids[gene_symbols %in% antigen_genes], ] <- mrna_mat[
    mrna_ids[gene_symbols %in% antigen_genes], ] +
    matrix(rep(scale(exposure$exposure_pm25),
               each = sum(gene_symbols %in% antigen_genes)),
           nrow = sum(gene_symbols %in% antigen_genes))

  # Create rowData for mRNA
  row_data_mrna <- S4Vectors::DataFrame(
    symbol = gene_symbols,
    row.names = mrna_ids
  )

  # Simulate proteomics with real gene names
  matched_genes <- sample(gene_symbols, n_proteins, replace = TRUE)
  prot_ids <- paste0("prot_", seq_len(n_proteins))

  prot_mat <- mrna_mat[mrna_ids[match(matched_genes, gene_symbols)], ] +
    matrix(rlnorm(n_proteins * n_samples,
                  meanlog = 7.5,
                  sdlog = 0.5),
           nrow = n_proteins)

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

  if(return_mae){
    res_list <- list(
      codebook = codebook,
      exposure = exposure,
      omics = omics,
      row_data = row_data
    )
    res <- create_expomicset(
      codebook = res_list$codebook,
      exposure = res_list$exposure,
      omics = res_list$omics,
      row_data = res_list$row_data
    )
  } else{
    res <- list(
      codebook = codebook,
      exposure = exposure,
      omics = omics,
      row_data = row_data
    )
  }

  return(res)

}
