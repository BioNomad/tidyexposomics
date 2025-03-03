library(GEOquery)
library(tidyverse)


meta <- pData(b$GSE40732_series_matrix.txt.gz)

meta_filt <- meta[meta$`asthma:ch1`=="TRUE",]

exp=exprs(b$GSE40732_series_matrix.txt.gz)
exp_filt=exp[,colnames(exp) %in% rownames(meta_filt)]



meth=exprs(a$GSE40576_series_matrix.txt.gz)
meth_filt=meth[,colnames(meth) %in% rownames(meta_filt)]

# ensure meth_filt and exp_filt match sample order of meta_filt
meth_filt=meth_filt[,match(rownames(meta_filt),colnames(meth_filt))]
# ensure exp_filt and meth_filt match sample order of meta_filt
exp_filt=exp_filt[,match(rownames(meta_filt),colnames(exp_filt))]

stopifnot(all(rownames(meta_filt) == colnames(meth_filt)))
stopifnot(all(rownames(meta_filt) == colnames(exp_filt)))


# grab methyl feature data
meth_features <- fData(a$GSE40576_series_matrix.txt.gz)

# grab exp feature data
exp_features <- fData(b$GSE40732_series_matrix.txt.gz)

# ensure that the features are in the same order
meth_features <- meth_features[match(rownames(meth_features), rownames(meth_filt)),]
exp_features <- exp_features[match(rownames(exp_features), rownames(exp_filt)),]

# ensure that the features match
stopifnot(all(rownames(meth_features) == rownames(meth_filt)))
stopifnot(all(rownames(exp_features) == rownames(exp_filt)))


# make expresssion se
exp_se <- SummarizedExperiment(
  assays = list(counts = as.matrix(exp_filt)),
  rowData = exp_features,
  colData = meta_filt
)

tnf=exp_se |> tidybulk() |> filter(.feature=="NM_000594") |> dplyr::select(.sample,counts)

generate_correlated_variable <- function(X, cor_target = 0.7, n = length(X), noise_sd = 5) {
  if (length(X) != n) {
    stop("Length of X must match specified number of observations (n).")
  }
  
  # Generate noise
  noise <- rnorm(n, mean = 0, sd = noise_sd)
  
  # Standardize X
  X_scaled <- scale(X)
  
  # Generate Y with specified correlation
  Y <- cor_target * X_scaled + sqrt(1 - cor_target^2) * scale(noise)
  
  # Rescale Y to have a similar range as X
  Y <- as.numeric(Y) * sd(X) + mean(X)
  
  return(Y)
}


tnf$fev1_fvc <- generate_correlated_variable(tnf$counts, cor_target = 0.5, n = length(tnf$counts), noise_sd = 5)

tnf$pm25 <- generate_correlated_variable(tnf$counts, cor_target = 0.3, n = length(tnf$counts), noise_sd = 5)

tnf$no2 <- generate_correlated_variable(tnf$counts, cor_target = 0.25, n = length(tnf$counts), noise_sd = 5)

tnf$ige_total <- generate_correlated_variable(tnf$counts, cor_target = 0.6, n = length(tnf$counts), noise_sd = 5)

tnf

