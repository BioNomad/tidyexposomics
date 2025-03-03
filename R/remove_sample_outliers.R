remove_sample_outliers <- function(
    expomicset,
    outliers=NULL
){
  
  if(!"pca" %in% names(expomicset@metadata)){
    stop("PCA not performed on the data. Please run 'pca_analysis' first.")
  }
  if(is.null(outliers)){
    outliers <-  expomicset@metadata$pca$outliers
  } else{
    outliers <- outliers
  }
  message("Removing outliers: ", paste(outliers, collapse=", "))
  expomicset <- expomicset[,!colnames(expomicset) %in% outliers]
  return(expomicset)
}
