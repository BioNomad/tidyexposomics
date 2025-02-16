remove_sample_outliers <- function(
    expOmicSet,
    outliers=NULL
){
  
  if(!"pca" %in% colnames(expOmicSet@metadata)){
    stop("PCA not performed on the data. Please run 'pca_analysis' first.")
  }
  if(is.null(outliers)){
    outliers <-  expOmicSet@metadata$pca$outliers
  } else{
    outliers <- outliers
  }
  message("Removing outliers: ", paste(outliers, collapse=", "))
  expOmicSet <- expOmicSet[,!colnames(expOmicSet) %in% outliers]
  return(expOmicSet)
}
