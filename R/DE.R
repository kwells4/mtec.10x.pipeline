#' Finds significant markers (DE genes) between clusters
#' 
#' This function allows you to find significant markers between
#' two clusters (or pairwise between any clusters you choose).
#' Uses seurat FindMarkers to find DE genes.
#' Requries a seurat object where FindClusters has already been
#' run. If no list is given, will perform pairwise analysis of all
#' clusters. If list is given, will perform pairwise analysis of list
#' given. Can optionally plot a volcano plot with any number of top
#' genes labeled on the plot. Returns a new seurat object containing
#' DE genes in the assay slot (in a list called DE)
#' 
#' @inheritParams DE_two_clusters
#' @param seurat_object A seurat object
#' @param cluster_list OPTIONAL A list of cluster names to perform DE
#' between. Defaults to NULL (all clusters in object)
#' @param adj_p_val OPTIONAL The pvalue to call for significance
#' Defaults to 0.05
#' @param logFC OPTIONAL The log full change to call significant.
#' Defaults to 1 (a 2 fold change)
#' @param plot_volcano OPTIONAL if a volcano plot should be made.
#' Defaults to FALSE.
#' @param num_genes OPTIONAL the number of genes to plot on the 
#' volcano plot if you choose to plot it. Defaults to 5.
#' @param test_use OPTIONAL which DE test to use
#' @param \dots paramaters to pass to DE_two_clusters
#' @keywords DE, volcano plot
#' @export
#' @examples
#' \dontrun{
#' significant_markers(mTEC.10x.data::mtec_trace)
#' significant_markers(mTEC.10x.data::mtec_trace, cluster_list = c(5, 6),
#' plot_volcano = TRUE)
#' }

significant_markers <- function(seurat_object,
                                cluster_list = NULL, ...){
  
  # Do DE of all clusters if no list is given
  if (is.null(cluster_list)){
    clusters <- seurat_object@ident
    cluster_vals <- levels(factor(clusters))
  } else {
    cluster_vals <- cluster_list
  }
  
  # Determine all pair-wise combinations of clusters given
  combinations <- utils::combn(cluster_vals, 2)

  # Run DE between all pairwise clusters and return a list of all DE genes
  de_list_all <- sapply(1:ncol(combinations), function(x)
    DE_two_clusters(seurat_object, combinations[1, x], combinations[2, x], ...))
  
  # Add the DE lists to the seurat object (append to an existing list or create a
    # new slot)
  if (is.null(seurat_object@assay$DE)){
    seurat_object@assay$DE <- de_list_all
  } else {
    previously_covered <- names(de_list_all) %in% names(seurat_object@assay$DE)
    previously_covered_DE <- names(de_list_all)[previously_covered]
    if(length(previously_covered_DE) > 0) {
      print(paste0("replacing ", previously_covered_DE, " in DE slot"))
      seurat_object@assay$DE <- modifyList(seurat_object@assay$DE, de_list_all)
    } else {
      seurat_object@assay$DE <- c(seurat_object@assay$DE, de_list_all)
    }
  }

  return(seurat_object)
}

#' Finds significant markers (DE genes) between clusters
#' 
#' This function allows you to find significant markers between
#' two clusters (or pairwise between any clusters you choose).
#' Uses seurat FindMarkers to find DE genes.
#' Requries a seurat object where FindClusters has already been
#' run. If no list is given, will perform pairwise analysis of all
#' clusters. If list is given, will perform pairwise analysis of list
#' given. Can optionally plot a volcano plot with any number of top
#' genes labeled on the plot. Returns a new seurat object containing
#' DE genes in the assay slot (in a list called DE)
#' 
#' @inheritParams plot_volcano
#' @param seurat_object A seurat object
#' @param ident_1 the identity of the first group to compare
#' @param ident_1 the identity of the second group to compare
#' @param adj_p_val OPTIONAL The pvalue to call for significance
#' Defaults to 0.05
#' @param logFC OPTIONAL The log full change to call significant.
#' Defaults to 1 (a 2 fold change)
#' @param plot_volcano OPTIONAL if a volcano plot should be made.
#' Defaults to FALSE.
#' @param test_use OPTIONAL which DE test to use
#' @param \dots paramaters to pass to plot_volcano
#' @keywords DE, volcano plot
#' @export
#' @examples
#' \dontrun{
#' significant_markers(mTEC.10x.data::mtec_trace)
#' significant_markers(mTEC.10x.data::mtec_trace, cluster_list = c(5, 6),
#' plot_volcano = TRUE)
#' }

DE_two_clusters <- function(seurat_object,
                            ident_1, ident_2,
                            adj_p_val    = 0.05,
                            logFC        = 1.0,
                            plot_volcano = FALSE,
                            test_use     = "wilcox",
                            ...){
    
  # Find markers between two clusters
  cluster_markers <- find_markers(seurat_object, ident_1, ident_2, adj_p_val, logFC,
                                    test_use)
    
  # Name the comparison
  slot_name <- paste0(as.character(ident_1), "v", as.character(ident_2))
    
  if (plot_volcano){
    title <- slot_name
      
    # Get a list of significant genes and make a volcano plot
    DE_sig <- plot_volcano(cluster_markers, title, num_genes)
      
  } else {
      
    # Make a dataframe consisting of significant and high fold change genes
    DE_sig <- cluster_markers[cluster_markers$group ==
                                "Significant&FoldChange", ]
  }
  DE_sig$genes <- rownames(DE_sig)
    
  # Add newest DE to list
  de_list <- list(DE_sig)
  names(de_list) <- slot_name
  return(de_list)
}

#' Finds gene markers between two clusters.
#' 
#' This function finds markers (DE genes) and returns a 
#' data frame with an additional column containing significance level
#' (significant/fold change). Meant to be called by significant markers,
#' but can work independently. Uses FindMarkers from seurat to call
#' DE genes
#' @param seurat_object A seurat object
#' @param ident_1 the name of the first cluster you wish to compare
#' @param ident_2 the name of the second cluster you wish to compare
#' @param adj_p_val OPTIONAL the p_value used to call significance.
#' Defaults to 0.05
#' @param logFC OPTIONAL the log fold change used to call significance.
#' Defaults to 1.0.
#' @param test_use OPTIONAL which DE test to use, default is wilcox
#' @keywords find markers
#' @export

find_markers <- function(seurat_object, ident_1, ident_2, adj_p_val = 0.05,
                         logFC = 1.0, test_use = "wilcox"){
  # Run FindMarkers from seurat on the two clusters that are called. Use
  # A random seed for consistent results.
  cluster_markers <- Seurat::FindMarkers(object          = seurat_object,
                                         ident.1         = ident_1,
                                         ident.2         = ident_2,
                                         random.seed     = 0,
                                         logfc.threshold = 0.01,
                                         test.use        = test_use)
    
  # Add a "group" column to dataframe. Make default "NotSignificant"
  cluster_markers["group"] <- "NotSignificant"
  
  # Find rows that have a low enough adjusted p_value to be significant
  # But which has a low log fold change. Change these values to "significant"
  # in the group column
  cluster_markers[which(cluster_markers["p_val_adj"] < adj_p_val &
                        abs(cluster_markers["avg_logFC"]) <=
                          logFC), "group"] <- "Significant"
  
  # Find rows that have a high p_value but large log fold change. Call these
  # "FoldChange"
  cluster_markers[which(cluster_markers["p_val_adj"] >= adj_p_val &
                        abs(cluster_markers["avg_logFC"]) >= 
                          logFC), "group"] <- "FoldChange"
  
  # Find rows with small p_value and high fold change. Call these
  # "Significant&FoldChange"
  cluster_markers[which(cluster_markers["p_val_adj"] < adj_p_val &
                        abs(cluster_markers["avg_logFC"]) > 
                          logFC), "group"] <- "Significant&FoldChange"
  
  return(cluster_markers)
}

#' Plots a volcano plot

#' This function plots a volcano plot and labels highest fold change
#' genes (both up and down). Works best when called by significant_markers,
#' but will work if given a full gene list.
#' @param marker_genes a data frame containing a list of genes and
#' a colum containing significance information. This df is made by the 
#' find_markers function
#' @param title the title of the plot
#' @param num_genes OPTIONAL the number of genes to label. Defalts to
#' 5.
#' @keywords volcano plot
#' @import ggplot2
#' @import ggrepel
#' @export

plot_volcano <- function(marker_genes, title, num_genes = 5){
  # Make a dataframe consisting of only the significant and high fold
  # change genes
  marker_genes_sig <- marker_genes[marker_genes$group ==
                                     "Significant&FoldChange", ]
  
  # Order these top hits by fold change and p_value, keep only the top 5
  top_peaks <- marker_genes_sig[with(marker_genes_sig,
                                     order(avg_logFC, 
                                           p_val_adj)),] [1:num_genes, ]
  
  
  # Order again by p_value and reverse log_fold change to find genes upregulated
  # in second cluster. Keep only the top 5 and add it to the first 5
  top_peaks <- rbind(top_peaks,
                     marker_genes_sig[with(marker_genes_sig,
                                           order(-avg_logFC,
                                                 p_val_adj)), ][1:num_genes, ])
  
  print(top_peaks)
  top_peaks$GENE_ID <- rownames(top_peaks)
  p <- ggplot2::ggplot(marker_genes,
                       ggplot2::aes_(x = ~avg_logFC, y = ~-log10(p_val_adj))) +
       ggplot2::geom_point(ggplot2::aes_(color = ~group)) +
       ggplot2::scale_color_manual(values = c("green", "grey", "black", "red")) +
       ggplot2::theme_classic(base_size = 16) +
       ggplot2::ggtitle(title) + 
       ggrepel::geom_text_repel(
         data = top_peaks,
         ggplot2::aes_(label = ~GENE_ID),
         size = 5,
         box.padding = unit(0.35, "lines"),
         point.padding = unit(0.75, "lines"),
         arrow = arrow(length = unit(0.01, "npc"))
        )
  print(p)
  return(marker_genes_sig)
}


#' Plot a heatmap of specific genes
#' This function plots a gene list and will label genes of interest if given.
#' Can only plot DE genes put in the DE slot after running significant markers.
#' @param seurat_object A seurat object
#' @param cell_color OPTIONAL a list of colors to color the cell labels at the
#' top of the plot. Defaults to NULL and uses color brewer set1
#' @param subset_list OPTIONAL a list of genes to subset the data. Defaults
#' to NULL
#' @param color_list OPTIONAL a list of genes to color gene names on the plot.
#' Defaults to NULL
#' @param order_cells OPTIONAL if cells should be ordered based on cluster number.
#' Defaults to TRUE
#' @param seed OPTIONAL the seed to use for reproducability. Defaults to 0.
#' @keywords heatmap
#' @export
#' @examples
#' \dontrun{
#' plot_heatmap(mTEC.10x.data::mtec_trace)
#' }
#' plot_heatmap(mTEC.10x.data::mtec_trace, subset_list = mTEC.10x.data::TFs,
#' color_list = c("Aire", "Hmgb2"), order_cells = TRUE)

plot_heatmap <- function(seurat_object, cell_color = NULL, subset_list = NULL,
  color_list = NULL, color_list2 = NULL, order_cells = TRUE,
  seed = 0){
  count_data <- seurat_object@data
  
  # Make a complete df of DE genes (combining all slots in list)
  DE_df <- data.table::rbindlist(seurat_object@assay$DE)
  DE_genes <- unique(DE_df$genes)
  
  # Subset the list if desired (ie by a list of specific genes)
  if (!is.null(subset_list)){
    DE_genes <- DE_genes[DE_genes %in% subset_list]
    count_data <- count_data[rownames(count_data) %in% DE_genes, ]
  }
  count_data <- as.matrix(count_data)
  
  # Center values to plot on heatmap
  count_data_heatmap <- t(scale(t(count_data), scale = FALSE))
  cluster <- as.data.frame(seurat_object@ident)
  names(cluster) <- "cluster_val"
  
  # Order cells by cluster
  if (order_cells){
    cluster <- cluster[order(cluster$cluster_val), , drop=FALSE]
    count_data_heatmap <- count_data_heatmap[, match(rownames(cluster),
                                                   colnames(count_data_heatmap))]
  }

  colors <- as.numeric(cluster$cluster_val)
  if (!is.null(cell_color)) {
    col1 <- cell_color
  } else {
    col1 <- RColorBrewer::brewer.pal(length(levels(cluster$cluster_val)), "Set1")
  }
  
  cols <- rep("black", nrow(count_data_heatmap))
  
  # Color some text red if desired
  if (!is.null(color_list)){
    cols[row.names(count_data_heatmap) %in% color_list] <- "red"
  }
  if (!is.null(color_list2)){
    cols[row.names(count_data_heatmap) %in% color_list2] <- "blue"
  }

  # Seed for reporducibility
  set.seed(seed)
  gplots::heatmap.2(count_data_heatmap,
                    density.info  = "none",
                    labCol        = FALSE,
                    Colv          = !order_cells,
                    colRow        = cols,
                    ColSideColors = col1[colors],
                    trace         = "none",
                    col           = grDevices::colorRampPalette(c("blue", "yellow")),
                    dendrogram    = "row")


}
