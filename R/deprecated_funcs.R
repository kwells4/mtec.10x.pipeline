#' Plots a UMAP and/or tSNE and/or a PCA
#' 
#' This function allows you to plot a tsne and/or a pca given a seurat
#' object and either a quality metric, cluster, or a gene. Can only plot genes,
#' anything from meta.data (quality) or cluster. Returns
#' a list with ggplot objects. Returns objects under slots pca, tsne, umap
#' depending on which of the three you ran.
#' @param mtec a Seurat object
#' @param col_by the name of a gene or column from meta data with which to
#' colour the plot. Can also be "cluster"
#' @param UMAP OPTIONAL if the plot should be a umap. Defaults to TRUE
#' @param tSNE OPTIONAL if the plot should be a tSNE. Defaults to FALSE
#' @param PCA OPTIONAL if the plot should be a PCA. Defaults to FALSE
#' @param color OPTIONAL the colors used to color the plots. For gene values, provide a list
#' of two colors. For discrete values, provide a list as long as the object you want to color
#' example: if you are coloring clusters and have 5 clusters, provide 5 colors. Default is
#' set1 from ggplot2 for discrete and red/blue for continuous.
#' @param save_plot OPTIONAL the path of the save file if the file should be saved.
#' Must end in .pdf or .png. Default is NULL with no saved plot but the plot is printed.
#' @param show_legend OPTIONAL if a legend for the plot should be shown. Default
#' is TRUE
#' @keywords tSNE, PCA, UMAP
#' @import ggplot2
#' @import RColorBrewer
#' @export
#' @examples 
#' tSNE_PCA(mTEC.10x.data::mtec_trace, "Aire")
#' tSNE_PCA(mTEC.10x.data::mtec_trace, "nUMI", tSNE = FALSE, PCA = TRUE)

tSNE_PCA <- function(mtec, col_by, UMAP = TRUE, tSNE = FALSE,
	                 PCA = FALSE, color = NULL, save_plot = NULL,
                   show_legend = TRUE) {
  .Deprecated("plotDimRed")
  if (UMAP){
  	print("umap")
  	plotDimRed(mtec = mtec, col_by = col_by, plot_type = "umap", color = color,
  		       save_plot = save_plot, show_legend = show_legend)
  }
  if (tSNE){
  	plotDimRed(mtec = mtec, col_by = col_by, plot_type = "tsne", color = color,
  		       save_plot = save_plot, show_legend = show_legend)
  }
  if (PCA){
  	plotDimRed(mtec = mtec, col_by = col_by, plot_type = "pca", color = color,
  		       save_plot = save_plot, show_legend = show_legend)
  }
}
