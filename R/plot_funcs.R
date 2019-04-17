#' Plots a UMAP and/or tSNE and/or a PCA
#' 
#' This function allows you to plot a tsne and/or a pca given a seurat
#' object and either a quality metric, cluster, or a gene. Can only plot genes,
#' anything from meta.data (quality) or cluster. Returns
#' a list with ggplot objects. Returns objects under slots pca, tsne, umap
#' depending on which of the three you ran.
#' @inheritParams discretePlots
#' @inheritParams continuousPlots
#' @param mtec a Seurat object
#' @param col_by the name of a gene or column from meta data with which to
#' colour the plot. Can also be "cluster"
#' @param plot_type OPTIONAL what type of diminsional reduction to plot. Must be 
#' a diminisional reduction listed in the dr slot of the seurat object. Default is
#' UMAP.
#' @param dims_use OPTIONAL
#' @return a ggplot object that is the dimensional reduction of your data
#' colored by the parameter of your choice.
#' @keywords tSNE, PCA, UMAP, diminsional reduction plot
#' @import ggplot2
#' @import RColorBrewer
#' @export
#' @examples 
#' plotDimRed(mTEC.10x.data::mtec_trace, "Aire")
#' plotDimRed(mTEC.10x.data::mtec_trace, "nUMI", tSNE = FALSE, PCA = TRUE)

plotDimRed <- function(mtec, col_by, plot_type = "umap", dims_use = NULL, ...) {
	# Determine where in Seurat object to find variable to color by
	if (col_by %in% rownames(mtec@data)){
		mtec_data <- mtec@data
		col_by_data <- as.data.frame(mtec_data[col_by, ])
	}else if (col_by %in% colnames(mtec@meta.data)){
		mtec_meta_data <- mtec@meta.data
		col_by_data <- as.data.frame(mtec_meta_data[, col_by, drop = FALSE])
	}else if (col_by == "cluster" | col_by == "Cluster"){
		col_by_data <- as.data.frame(mtec@ident)
	}else {
		stop("col_by must be a gene, metric from meta data or 'cluster'")
	}

	# Make the name in the data frame the same regardless of what it was originally
	names(col_by_data) <- "colour_metric"

  if (is.null(dims_use)){
    dims_use <- c(1,2)
  }
  # Make a data frame based on the cell embeddings from the plot type of choice
  if (plot_type %in% names(mtec@dr)){
    plot_coord <- mtec@dr[[plot_type]]@cell.embeddings
    plot_names <- colnames(plot_coord)
    ndims <- length(plot_names)
    plot_cols <- lapply(dims_use, function(x){
      if (x > ndims) {
        stop("dims_use must be equal to or less than number of dimensions")
      } else {
        plot_col <- plot_names[x]
        return(plot_col)
      }
    })
    plot_cols <- unlist(plot_cols)
    plot_coord <- plot_coord[colnames(plot_coord) %in% plot_cols, ]
    axis_names <- colnames(plot_coord)
    colnames(plot_coord) <- c("dim1", "dim2")
    plot_df <- merge(plot_coord, col_by_data, by = "row.names")
  } else {
    stop("plot type must be a dimensional reduction in dr slot of Seurat object")
  }
	# Plot as discrete
	if (!is.numeric(col_by_data$colour_metric)){
		return_plot <- discretePlots(plot_df, axis_names = axis_names,
                                 col_by = col_by, ...)

	# Plot as continuous
	}else{
		return_plot <- continuousPlots(plot_df, axis_names = axis_names,
                                   col_by = col_by, ...)
	}
  return(return_plot)
}

#' Plots discrete colours on a umap, tsne, and or a PCA
#' 
#' This function allows you to plot discrete colours on a tsne and or a pca
#' given a seurat object and dataframe. Must have already run FindClusters by
#' seurat. Called by plotDimRed and works best through that function. Returns
#' a ggplot object.
#' @param plot_df a data_frame containing values used for plotting. This is df
#' is made by the plotDimRed function. x and y values should be in columns labeled
#' dim1 and dim2. discrete descriptions of how to color are given in the column
#' colour_metric.
#' @param col_by what to use to colour the plot. Can be anything from meta.data
#' or clusters. This will be the name of the color bar.
#' @param axis_names OPTIONAL What to name the axes. This will be the default
#' seurat name of the dimensional reduction if called by plotDimRed. Default is
#' c("dim1", "dim2")
#' @param color OPTIONAL what color should be used to color the plots. Must be the same
#' number of items as col_by. Defaults to set1 from RColorBrewer
#' @param save_plot OPTIONAL the path of the save file if the file should be saved.
#' Must end in .pdf or .png. Default is NULL with no saved plot.
#' @param show_legend OPTIONAL if a legend for the plot should be shown. Default
#' is TRUE
#' @return a ggplot object that is the dimensional reduction of your data
#' colored by the parameter of your choice.
#' @keywords tSNE, PCA, UMAP, clusters, dimensional reduction plot
#' @import ggplot2
#' @import RColorBrewer
#' @export

discretePlots <- function(plot_df, col_by, axis_names = c("dim1", "dim2"),
                          color = NULL, save_plot = NULL, show_legend = TRUE){
 	base_plot <- ggplot2::ggplot(data = plot_df,
                               ggplot2::aes_(~dim1, ~dim2))

    # Add colors based on metric chosen
  base_plot <- base_plot +
      ggplot2::geom_point(ggplot2::aes_(colour = ~colour_metric),
                          show.legend = show_legend) 

  nColors <- length(levels(factor(plot_df$colour_metric)))

  # Color based on RColorBrewer if own palette isn't chosen
  if (is.null(color)) {
    base_plot <- base_plot + ggplot2::scale_color_manual(
      values = grDevices::colorRampPalette(
          RColorBrewer::brewer.pal(9, "Set1"))(nColors), name = col_by)
  } else {
    base_plot <- base_plot + ggplot2::scale_color_manual(values = color, name = col_by)
  }

  if (!(is.null(save_plot))){
    ggplot2::ggsave(save_plot, plot = base_plot)
  }
  return(base_plot)
}


#' Plots continuous colours on a tsne and or a PCA
#' 
#' This function allows you to plot continuous colours on a tsne and or a pca
#' given a seurat object and dataframe. Must have already run FindClusters by
#' seurat. CCalled by plotDimRed and works best through that function. Returns
#' a ggplot object.
#' @param plot_df a data_frame containing values used for plotting. This is df
#' is made by the plotDimRed function. x and y values should be in columns labeled
#' dim1 and dim2. discrete descriptions of how to color are given in the column
#' colour_metric.
#' @param col_by what to use to colour the plot. Can be anything from meta.data
#' or clusters. This will be the name of the color bar.
#' @param axis_names OPTIONAL What to name the axes. This will be the default
#' seurat name of the dimensional reduction if called by plotDimRed. Default is
#' c("dim1", "dim2")
#' @param color OPTIONAL what color should be used to color the plots. Must be a list
#' of two colors, the first color as low and the second color as high. Default is low
#' is blue and high is red.
#' @param save_plot OPTIONAL the path of the save file if the file should be saved.
#' Must end in .pdf or .png. Default is NULL with no saved plot.
#' @param show_legend OPTIONAL if a legend for the plot should be shown. Default
#' is TRUE
#' @return a ggplot object that is the dimensional reduction of your data
#' colored by the parameter of your choice.
#' @keywords tSNE, PCA, clusters, dimensional reduction plot
#' @import ggplot2
#' @import RColorBrewer
#' @export

continuousPlots <- function(plot_df, col_by, axis_names = c("dim1", "dim2"),
                            color = NULL, save_plot = NULL, show_legend = TRUE){
	if (is.null(color)) {
	  low <- "#00AFBB"
	  high <- "#FC4E07"
	} else {
	  low <- color[1]
	  high <- color[2]
	}

	base_plot <- ggplot2::ggplot(data = plot_df, ggplot2::aes_(~dim1, ~dim2))

	base_plot <- base_plot +
	  ggplot2::geom_point(ggplot2::aes_(colour = ~colour_metric),
                        show.legend = show_legend) +
 		ggplot2::scale_color_gradient(low = low,
 		                              high = high, name = col_by) 
 		 
  if (!(is.null(save_plot))){
    ggplot2::ggsave(save_plot, plot = base_plot)
  }
  return(base_plot)
}

#' Plots three plots in one image
#' 
#' This function allows you to plot three plots on one image that are violin
#' or jitter plots given a seurat object. Must have already run FindClusters by
#' seurat.
#' @param seurat_object a Seurat object
#' @param geneset a list of three genes to plot
#' @param cell_cycle OPTIONAL if to color plots based on cell cycle scores.
#' Defaults to FALSE
#' @param plot_jitter OPTIONAL if the plot should be a jitter plot. Defaults to TRUE
#' @param plot_violin OPTIONAL if the plot should be a violin. Defaults to FALSE
#' @param jitter_and_violin OPTIONAL if the plot should be a jitter on top of a violin.
#' Defaults to FALSE
#' @param color OPTIONAL what colors should be used to color the plots. Defaults to 
#' RColorBrewer set 9.
#' @param sep_by OPTIONAL what to use to separate on the x axis. Defaults to cluster
#' @param save_plot OPTIONAL the path of the save file if the file should be saved.
#' Must end in .pdf or .png. Default is NULL with no saved plot but the plot is printed.
#' @param nrow OPTIONAL number of rows for output plots default is number of plots
#' @param ncol OPTIONAL number of columns for output plots default is 1
#' @param group_color OPTIONAL if a plot of group colors should be output. Default
#' is true
#' @return a ggplot object containing the jitter or violin plots arranged with
#' grid extra
#' @keywords violin, jitter, trio
#' @import gridExtra
#' @export

trio_plots <- function(seurat_object, geneset, cell_cycle = FALSE,
                       plot_jitter = TRUE, plot_violin = FALSE,
                       jitter_and_violin = FALSE, color = NULL,
                       sep_by = "cluster", save_plot = NULL,
                       nrow = NULL, ncol = NULL, group_color = TRUE){
  gene_list_stage <- c()
  if (plot_jitter) {
    if (group_color) {
      # Make a jitter plot based on expression of each gene given in the gene
      # set color by stage
      for (gene in geneset) {
        gene_stage <- jitter_plot(seurat_object, gene, sep_by,
                                color = color)
    
        # Add this plot object to a list
        gene_list_stage[[gene]] <- gene_stage
      }
    
      # Make a plot consisting of all plots made above
      full_plot <- gridExtra::grid.arrange(grobs = gene_list_stage,
                                           nrow = length(geneset))
    }
    # Make jitter plots colored by cell cycle stage
    if(cell_cycle){
      gene_list_cycle <- c()
      for (gene in geneset) {
        gene_cycle <- jitter_plot(seurat_object, gene, "stage", "cycle_phase",
                                 color = c("black", "red", "purple"))
    
        gene_list_cycle[[gene]] <- gene_cycle
      }
    
      # Arrange all plots into one figure
      full_plot <- gridExtra::grid.arrange(grobs = gene_list_cycle,
                                           nrow = length(geneset))
    }
  }
  if (plot_violin || jitter_and_violin) {
    for (gene in geneset) {
      gene_stage <- violin_plot(seurat_object, gene, sep_by,
                                color = color,
                                plot_jitter = jitter_and_violin)
      
      # Add this plot object to a list
      gene_list_stage[[gene]] <- gene_stage
    }
    
    # Make a plot consisting of all plots made above
    if (is.null(nrow)){
      nrow <- length(geneset)
    }
    if (is.null(ncol)){
      ncol <- 1
    }
    full_plot <- gridExtra::grid.arrange(grobs = gene_list_stage,
                                         nrow = nrow, ncol = ncol)
    
  }
  if (!(is.null(save_plot))){
    ggplot2::ggsave(save_plot, plot = full_plot)
  }
  return(full_plot)
}

#' Plots a jitter plot
#' 
#' This function allows you to plot a jitter plot given a seurat object. Must have
#' already run FindClusters by seurat.
#' @param seurat_object a Seurat object
#' @param y_val a gene (or any conitnuous value from metadata) to plot
#' @param x_val cluster (or any discrete value from metadata) to plot
#' @param col_by OPTIONAL what should be used to color the plots. Defaults to cluster
#' @param color OPTIONAL what colors should be used to color the plots. Defaults to 
#' RColorBrewer set 9.
#' @keywords jitter
#' @import ggplot2
#' @import RColorBrewer
#' @export

jitter_plot <- function(seurat_object, y_val, x_val,
                        col_by = NULL, color = NULL) {
  plot_data <- make_plot_df(seurat_object, y_val, x_val,
                            col_by, color)
  # Determine the number of different colors needed.
  nColors <- length(unique(plot_data$col_by))
  
  # Plot the main plot
  plot_base <- ggplot2::ggplot(data = plot_data, ggplot2::aes_(~x_value, 
                                                               ~y_value,
                                                               color = ~col_by)) +
    #ggplot2::theme_classic() + 
    ggplot2::ylab(y_val) + ggplot2::xlab(x_val) +
    ggplot2::geom_point() + ggplot2::geom_jitter(shape = 16) +
    ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_blank())
  
  
  
  if (is.null(color)) {
    plot_base <- plot_base +
      ggplot2::scale_color_manual(values =
                                    (grDevices::colorRampPalette(
                                      RColorBrewer::brewer.pal(9, 
                                                               "Set1")))(nColors), name = col_by)
  } else {
    plot_base <- plot_base + ggplot2::scale_color_manual(values = color, 
                                                         name = col_by)  
  }
  
  return(plot_base)
}  


#' Plots a violin plot
#' 
#' This function allows you to plot a violin plot given a seurat object. Must have
#' already run FindClusters by seurat.
#' @param seurat_object a Seurat object
#' @param y_val a gene (or any conitnuous value from metadata) to plot
#' @param x_val cluster (or any discrete value from metadata) to plot
#' @param col_by OPTIONAL what should be used to color the plots. Defaults to cluster
#' @param color OPTIONAL what colors should be used to color the plots. Defaults to 
#' RColorBrewer set 9.
#' @param plot_jitter OPTIONAL if a jitter plot should be plotted on top of the violin
#' plot. Defaults to FALSE
#' @keywords violin
#' @import ggplot2
#' @import RColorBrewer
#' @export

violin_plot <- function(seurat_object, y_val, x_val,
                        col_by = NULL, color = NULL,
                        plot_jitter = FALSE) {
  plot_data <- make_plot_df(seurat_object, y_val, x_val,
                            col_by, color)
  # Determine the number of different colors needed.
  nColors <- length(unique(plot_data$col_by))
  
  # Plot the main plot
  plot_base <- ggplot2::ggplot(data = plot_data, ggplot2::aes_(~x_value, 
                                                               ~y_value,
                                                               fill = ~col_by)) +
    #ggplot2::theme_classic() + 
    ggplot2::ylab(y_val) + ggplot2::xlab(x_val) +
    ggplot2::geom_violin(scale = "width") +
    ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_blank())
  
  if (plot_jitter) {
    plot_base <- plot_base + ggplot2::geom_jitter(shape = 16)
  }
  
  if (is.null(color)) {
    plot_base <- plot_base +
      ggplot2::scale_fill_manual(values =
                                    (grDevices::colorRampPalette(
                                      RColorBrewer::brewer.pal(9, 
                                                               "Set1")))(nColors), name = col_by)
  } else {
    plot_base <- plot_base + ggplot2::scale_fill_manual(values = color, 
                                                         name = col_by)  
  }
  
  return(plot_base)
}


#' Plots a data frame that is easily placed into ggplot2
#' 
#' This function allows you to take a seurat object and use many pieces from it to make
#' a data frame that easily interacts with ggplot2
#' @param seurat_object a Seurat object
#' @param y_val a gene (or any conitnuous value from metadata) to plot
#' @param x_val cluster (or any discrete value from metadata) to plot
#' @param col_by OPTIONAL what should be used to color the plots. Defaults to cluster
#' @keywords ggplot2, data frame
#' @export

make_plot_df <- function(seurat_object, y_val, x_val,
                          col_by = NULL) {
  # Add y_value to a data frame used for plotting. This value can be a gene
  # or a value from meta data like nGene
  if (y_val %in% rownames(seurat_object@data)) {
    # Add column containing normalized UMI counts
    # Should be able to fix this. Take out as and do y_val = then don't need names()
    plot_data <- as.data.frame(seurat_object@data[y_val, ])
    #plot_data <- data.frame("y_value" = seurat_object@data[y_val, ])
  }
  else if (y_val %in% colnames(seurat_object@meta.data)) {
    # Add column contining meta data info.
    plot_data <- as.data.frame(seurat_object@meta.data[, y_val, drop = FALSE])
    #plot_data <- data.frame("y_value" = seurat_object@meta.data[, y_val,
    #                                                           drop = FALSE])
  }
  else {
    stop("y_val must be a gene, metric from meta data")
  }
  # Name the column
  names(plot_data) <- "y_value"
  
  # Add a column contining the x_value. This should be something discrete
  # Like timepoint or cluster
  if (x_val %in% colnames(seurat_object@meta.data)) {
    # Should be able to fix this. Take out as and do x_val = then don't need names()
    x_plot_data <- as.data.frame(seurat_object@meta.data[, x_val, drop = FALSE])
    #x_plot_data <- data.frame("x_value" = seurat_object@meta.data[, x_val,
     #                                                             drop = FALSE])
    plot_data <- merge(plot_data, x_plot_data, by = "row.names")
    rownames(plot_data) <- plot_data$Row.names
    plot_data$Row.names <- NULL
  } else if (x_val == "cluster") {
    x_plot_data <- as.data.frame(seurat_object@ident)
    #x_plot_data <- data.frame("x_value" = seurat_object@ident)
    plot_data <- merge(plot_data, x_plot_data, by = "row.names")
    rownames(plot_data) <- plot_data$Row.names
    plot_data$Row.names <- NULL
  } else {
    stop("x_val must be a metric from meta data or 'cluster'")
  }
  
  # Name the appropriate column of the plotting data
  names(plot_data)[2] <- "x_value"
  plot_data <- plot_data[match(rownames(seurat_object@meta.data),
                               rownames(plot_data)), ]
  
  # Determine how to color the plot. Default is the x_value but can be any
  # discrete value.
  if (is.null(col_by)) {
    plot_data$col_by <- plot_data$x_value
  } else if (col_by %in% colnames(seurat_object@meta.data)) {
    plot_data$col_by <- seurat_object@meta.data[ , col_by]
  } else if (col_by == "cluster") {
    plot_data$col_by <- seurat_object@ident
  } else {
    stop("x_val must be a metric from meta data or 'cluster'")
  }
  return(plot_data)
}

#' Makes a dataframe from multiple seurat objects for one variable
#' 
#' This function allows you to take multiple seurat objects and combine
#' them for one variable (like cell type) to make a plot of population percents
#' Could be improved to make the entire df and call the plotting function
#' @param seurat_object a Seurat object
#' @param sample_name the name of one sample to add to the plot
#' @param stage_df_all the data frame containing all samples
#' @keywords population percents
#' @export

populations_dfs <- function(seurat_object, sample_name, stage_df_all){
  stage_df <- data.frame(table(seurat_object@meta.data$stage))
  names(stage_df) <- c("stage", "count")
  stage_df$percent <- stage_df$count / sum(stage_df$count) * 100
  stage_df$sample <- sample_name
  if(is.null(stage_df_all)){
    stage_df_all <- stage_df
  } else {
    stage_df_all <- rbind(stage_df_all, stage_df)
  }
  return(stage_df_all)
}

#' Makes a plot of population percents
#' 
#' This function allows you to make a population percent plot from many different
#' samples.
#' @param stage_df_all a data frame contining the info to plot
#' @param color how to color the populations
#' @keywords population percents
#' @import ggplot2
#' @export

population_plots <- function(stage_df_all, color){
  
  plot_base <- ggplot2::ggplot(data = stage_df_all, ggplot2::aes_(x = ~sample,
                                                                  y = ~percent,
                                                                  fill = ~stage)) +
    #ggplot2::theme_classic() + 
    ggplot2::xlab("frequency")  +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = color, name = "stage")
  
  return(plot_base)

}


#' Plots a density plot
#' 
#' This function allows you to plot density plot given a seurat
#' object and either a quality metric or a gene.
#' @param mtec a Seurat object
#' @param y_val the name of a gene or quality metric with which to make the y
#' axis of the plot
#' @param boxplot OPTIONAL if a boxplot should be made. Defaults to TRUE
#' @param vioplot OPTIONAL if a violin plot should be made. Defaults to FALSE
#' @param density OPTIONAL if a density plot should be made. Defaults to FALSE
#' @import ggplot2
#' @import RColorBrewer
#' @export
#' @examples 
#' Exp_plots(mTEC.10x.data::mtec_trace, "Aire")
#' Exp_plots(mTEC.10x.data::mtec_trace, "nUMI")


Exp_plots <- function(mtec, y_val,
						 boxplot = TRUE, vioplot = FALSE, density = FALSE){
  if (y_val %in% rownames(mtec@data)){
    mtec_data <- mtec@data
    y_val_data <- as.data.frame(mtec_data[y_val, ])
  }else if (y_val %in% colnames(mtec@meta.data)){
    mtec_meta_data <- mtec@meta.data
    y_val_data <- as.data.frame(mtec_meta_data[, y_val, drop = FALSE])
  }else {
    stop("col_by must be a gene, metric from meta data or 'cluster'")
  }
	names(y_val_data) <- "y_value"

  y_val_data$cluster <- mtec@ident

	plot_base <- ggplot2::ggplot(data = y_val_data,
	                             ggplot2::aes_(~cluster, ~y_value)) +
	  #ggplot2::theme_classic() +
		ggplot2::ylab(y_val) +
		ggplot2::scale_colour_brewer(palette = "Set1", name = "cluster")	

	if (boxplot) {
		
		p_boxplot <- ggplot2::ggplot(data = y_val_data,
		                             ggplot2::aes_(~cluster, ~y_value)) +
		  ggplot2::geom_boxplot(ggplot2::aes_(colour = ~cluster)) +
			ggplot2::geom_jitter(shape = 16, position = position_jitter(0.2),
			                     ggplot2::aes_(color = ~cluster)) +
			#ggplot2::theme_classic() +
			ggplot2::ylab(y_val) +
			ggplot2::scale_colour_brewer(palette = "Set1", name = "cluster")

		print(p_boxplot)
	}
	if (vioplot) {
		p_vioplot <- ggplot2::ggplot(data = y_val_data,
		                             ggplot2::aes_(~cluster, ~y_value)) +
		  ggplot2::geom_violin(ggplot2::aes_(colour = ~cluster)) +
			ggplot2::geom_jitter(shape = 16, position = position_jitter(0.2),
			                     ggplot2::aes_(color = ~cluster)) +
			#ggplot2::theme_classic() +
			ggplot2::ylab(y_val) +
			ggplot2::scale_colour_brewer(palette = "Set1", name = "cluster")

		print(p_vioplot)
	}
	if (density){
		p_density <- ggplot2::ggplot(y_val_data, aes_(~y_value,
		                                              colour = ~cluster)) +
		  ggplot2::geom_density() +
		  ggplot2::scale_colour_brewer(palette = "Set1",
		                               ggplot2::aes_(name = ~y_value)) + 
		 	#ggplot2::theme_classic() +
		  ggplot2::xlab(y_val)
		
		print(p_density)
	}
	if (!boxplot & !vioplot & !density){
		stop("Must select one plot type")
	}
}


#' Plots PC loadings
#' 
#' This function allows you to plot PCA loadings for any number of genes for any
#' PC given a Seurat object
#' @param mtec a Seurat object
#' @param PC_val OPTIONAL which PC loadings to plot. Defaults to "PC1"
#' @param nLoadings OPTIONAL the number of loadings to plot. Defaults to 10
#' @param right_pull OPTIONAL plot genes pulling PC to the right. Defaults to
#' TRUE.
#' @param left_pull OPTIONAL plot genes pulling PC to the left. Defaults to
#' TRUE.
#' @param colour_by OPTIONAL a list of genes used for colouring. Defaults to
#' NULL
#' @param full_loadings OPTIONAL if all loadings should be used. If false, only
#' hvgs used
#' to create the PCA will be plotted. Defaults to FALSE.
#' @keywords tSNE, PCA, clusters
#' @import ggplot2
#' @import RColorBrewer
#' @import dplyr
#' @export
#' @examples 
#' PCA_loadings(mTEC.10x.data::mtec_trace)
#' PCA_loadings(mTEC.10x.data::mtec_trace, PC_val = "PC2", nLoadings = 20,
#'				colour_by = mTEC.10x.data::tuft_markers, left_pull = FALSE)

PCA_loadings <- function(mtec, PC_val = "PC1", nLoadings = 10,
                         left_pull = TRUE, right_pull = TRUE,
                         colour_by = NULL, full_loadings = FALSE){
	if (full_loadings){
		pc_loadings <- data.frame(mtec@dr$pca@gene.loadings.full)
	} else {
		pc_loadings <- data.frame(mtec@dr$pca@gene.loadings)
	}
	pc_loadings$gene <- rownames(pc_loadings)
	pc_loadings$PC <- pc_loadings[[PC_val]]
	if (left_pull && right_pull){
    	pc_loadings <- rbind(pc_loadings %>% top_n(nLoadings, pc_loadings$PC),
                               pc_loadings %>% top_n(-nLoadings, pc_loadings$PC))
	} else if (right_pull && !left_pull) {
    	pc_loadings <- pc_loadings %>% top_n(nLoadings, pc_loadings$PC)
	} else if (left_pull && !right_pull){
	  pc_loadings <- pc_loadings %>% top_n(-nLoadings, pc_loadings$PC)
	} else {
	  stop("either right_pull or left_pull must be TRUE")
	}
	pc_loadings <- pc_loadings[order(pc_loadings$PC), ]
	cols <- rep("black", nrow(pc_loadings))
	if (!is.null(colour_by)) {
    	cols[pc_loadings$gene %in% colour_by] <- "red"
	}
	print(ggplot2::ggplot(pc_loadings, ggplot2::aes_(x = ~PC, y = ~gene)) +
          ggplot2::geom_point(colour = "blue") +
          ggplot2::scale_y_discrete(limits = pc_loadings$gene) +
          ggplot2::labs(x = PC_val) +
          ggplot2::theme(axis.text.y = ggplot2::element_text(colour = cols)))
}
