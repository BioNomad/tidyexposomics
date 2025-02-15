random_heatmap <- function(nrow,
                           ncol,
                           low="white",
                           high="blue"){
  # code to generate a small heatmap sample
  library(ggplot2)
  library(reshape2)
  
  # create a matrix of random numbers
  set.seed(123)
  m <- matrix(runif(nrow*ncol), nrow=nrow,ncol = ncol)
  
  # create a heatmap
  ggplot(melt(m), aes(Var1,
                      Var2,
                      fill = value)) +
    geom_tile() +
    scale_fill_gradient(
      low = low,
      high = high)+
    theme_void()+
    theme(legend.position = "none")
}
