#' Runs monocle on cells
#' 
#' This function allows you to run monocle on your samples. It also will add the
#' monocle output into the seurat object. While this can be run on a local, I
#' recommend only running this on a cluster and setting cores to at least 10.
#' If you only have a local, just use the mtec_trace object as it already
#' includes the data from monocle.
#' @param mtec a Seurat object
#' @param quality_plots OPTIONAL if quality plots (including density plot of
#' expression values, PC variance explained, and tSNE coloured by monocle and
#' seurat clusters) should be plotted. Defaults to FALSE.
#' @param cores OPTIONAL the number of cores to use when running
#' differentialGeneTest. very slow when run on one core. Defaults to 1.
#' @param seed OPTIONAL seed for reproducibility. Defaults to 0
#' @keywords monocle
#' @export
#' @examples
#' \dontrun{
#' run_monocle(mTEC.10x.data::mtec_trace, cores = 10)
#' }

run_monocle <- function(mtec, quality_plots = FALSE, cores = 1, seed = 0){
  
  if (!requireNamespace("monocle", quietly = TRUE)){
		stop("Package \"monocle\" needed for this function to work. Please install it.",
			call. = FALSE)
  }
  if (!requireNamespace("reshape", quietly = TRUE)){
		stop("Package \"reshape\" needed for this function to work. Please install it.",
			call. = FALSE)
  }
  if (!requireNamespace("BiocGenerics", quietly = TRUE)){
    stop("Package \"BiocGenerics\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

	# for seurat objects, a newCellDataSet is made with the negbinomial.size() as 
	# the expressionFamily. It also takes the raw data. If the raw.data and meta.data
	# names have been moved, reorder the raw.data to fit the meta.data
	mtec_dev <- monocle::importCDS(mtec)
	mtec_dev <- BiocGenerics::estimateSizeFactors(mtec_dev)
	mtec_dev <- BiocGenerics::estimateDispersions(mtec_dev)
	mtec_dev <- monocle::detectGenes(mtec_dev, min_expr = 0.1)
	gene_expressed <- Biobase::fData(mtec_dev)$num_cells_expressed >=3
	expressed_genes <- row.names(subset(Biobase::fData(mtec_dev),
	                                    subset = gene_expressed))
	
	Biobase::fData(mtec_dev)$use_for_ordering<-
	  Biobase::fData(mtec_dev)$num_cells_expressed > 0.05 * ncol(mtec_dev)
	
	set.seed(seed)
	mtec_dev <- monocle::reduceDimension(mtec_dev,
	                                     max_compoments   = 2,
	                                     norm_method      = "log",
	                                     num_dim          = 10,
	                                     reduction_method = "tSNE",
	                                     verbose          = T)

	set.seed(seed) 
	mtec_dev <- monocle::clusterCells(mtec_dev, verbose = F)

	set.seed(seed)
	clustering_DEG_genes <- monocle::differentialGeneTest(mtec_dev[expressed_genes, ],
	                                             fullModelFormulaStr = "~Cluster",
	                                             cores = cores)

	mtec_dev_ordering_genes <- row.names(clustering_DEG_genes)[order(
										clustering_DEG_genes$qval)][1:1000]

	mtec_dev <- monocle::setOrderingFilter(mtec_dev,
	                                       ordering_genes =
	                                         mtec_dev_ordering_genes)

	set.seed(seed)
	mtec_dev <- monocle::reduceDimension(mtec_dev, method = "DDRTree")

	set.seed(seed)
	mtec_dev <- monocle::orderCells(mtec_dev)

	if (quality_plots){
		L <- log(Biobase::exprs(mtec_dev[expressed_genes,]))
		melted_dens_df <- reshape::melt(Matrix::t(scale(Matrix::t(L))))
		
		density <- ggplot2::ggplot(melted_dens_df, aes_(~value)) +
		  ggplot2::geom_density() +
		  ggplot2::stat_function(fun = stats::dnorm, size = 0.5, color = "red") +
		  ggplot2::xlab("Standardized log(FPKM)") +
		  ggplot2::ylab("Density")
		
		print(density)

		print(monocle::plot_pc_variance_explained(mtec_dev, return_all = F))

		print(monocle::plot_cell_clusters(mtec_dev, color_by = 'Cluster'))

		#Colour by seurat cluster
		print(monocle::plot_cell_clusters(mtec_dev, color_by = 'res.0.6') +
			ggplot2::scale_color_brewer(palette = "Set1"))
	}
	return(mtec_dev)
}

#' Determines significant genes from a monocle object
#' 
#' This function finds genes that drive pseudotime. While this can be run on a
#' local, I recommend only running this on a cluster and setting cores to at
#' least 10. If you only have a local, just use the mtec_trace object as it
#' already includes the data from monocle.
#' @param mtec_dev a CellDataSet object
#' @param cores OPTIONAL the number of cores to use when running
#' differentialGeneTest. very slow when run on one core. Defaults to 1.
#' @keywords monocle
#' @export
#' @examples 
#' \dontrun{
#' monocle_genes(mTEC.10x.data::mtec_trace_monocle, cores = 10)
#' }

monocle_genes <- function(mtec_dev, cores = 1){
  
  if (!requireNamespace("monocle", quietly = TRUE)){
    stop("Package \"monocle\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
	set.seed(0)
	diff_test_res_all <- monocle::differentialGeneTest(mtec_dev,
	                                                   fullModelFormulaStr =
	                                                     "~sm.ns(Pseudotime)",
	                                                   cores = cores)

	diff_test_res_all$pval_adj <- diff_test_res_all$pval * nrow(diff_test_res_all)

	sig_genes <- diff_test_res_all[diff_test_res_all$pval_adj < 0.05,]

	sig_gene_names <- row.names(sig_genes)


	set.seed(0)
	all_genes_heatmap <- monocle::plot_pseudotime_heatmap(mtec_dev[sig_gene_names, ],
	                                                      cores = cores,
		                                                    show_rownames = FALSE,
		                                                    return_heatmap = TRUE)

	all_genes_df <- as.data.frame(stats::cutree(all_genes_heatmap$tree_row, k = 6))
	colnames(all_genes_df) <- "Cluster"
	all_genes_df$Gene <- rownames(all_genes_df)

	return(all_genes_df)
}

#' Runs diffusion pseudotime (dpt) on cells
#' 
#' This function allows you to run dpt on your samples. It also will add the
#' dpt output into the seurat object. While this can be run on a local, I
#' recommend only running this on a cluster. If you only have a just use the
#' mtec_trace object as it already includes the data from dpt.
#' @param mtec a Seurat object
#' @param seed OPTIONAL seed for reporducibility. Defaults to 0.
#' @keywords dpt
#' @export
#' @examples 
#' \dontrun{
#' run_dpt(mTEC.10x.data::mtec_trace)
#' }


run_dpt <- function(mtec, seed = 0){
  if (!requireNamespace("destiny", quietly = TRUE)){
    stop("Package \"destiny\" needed for this function to work. Please install 
         it.", call. = FALSE)
  }
	meta_data <- mtec@meta.data
	meta_data$cluster <- mtec@ident
	mtec_data <- mtec@data
	mtec_data <- as.matrix(mtec_data)
	pd <- methods::new('AnnotatedDataFrame', data = meta_data)
	mtecExpreSet <- Biobase::ExpressionSet(assayData=mtec_data, phenoData=pd)

	set.seed(seed)
	dm <- destiny::DiffusionMap(mtecExpreSet)

	set.seed(seed)
	dpt <- destiny::DPT(dm)

	return(dpt)
}

#' Adds the pseudotime information to a Seurat object
#' 
#' This function adds information from a pseudotime object to the seurat meta
#' data. Some parts of this function were borrowed from
#' maehrlab/thymusatlastools
#' @param mtec a Seurat object
#' @param pt_obj the pseudotime object
#' @param method OPTIONAL the method of pseudotime used. Can be "monocle" or
#' "dpt". Defaults to "monocle"
#' @param root_state OPTIONAL if you want to change the root state for
#' pseudotime. defaults to NULL.
#' @keywords monocle
#' @export
#' @examples 
#' add_pseudotime(mTEC.10x.data::mtec_trace, mTEC.10x.data::mtec_trace_monocle)
#' add_pseudotime(mTEC.10x.data::mtec_trace, mTEC.10x.data::mtec_trace_dpt,
#' method = "dpt")

add_pseudotime <- function(mtec, pt_obj, method = "monocle", root_state = NULL) {
  if (method == "monocle"){
    if (!requireNamespace("monocle", quietly = TRUE)){
      stop("Package \"monocle\" needed for this function to work. Please install 
           it.", call. = FALSE)
    }
		if (!is.null(root_state)){
			pt_obj <- monocle::orderCells(pt_obj, root_state = root_state)
		}
		pt_info <- data.frame(
		  monocle_dm_1 = monocle::reducedDimS(pt_obj)[1, ],
		  monocle_dim_2 = monocle::reducedDimS(pt_obj)[2, ],
		  monocle_branch = as.character(Biobase::pData(pt_obj)$State),
		  monocle_pt = Biobase::pData(pt_obj)$Pseudotime)
	}else if (method == "dpt"){
	  if (!requireNamespace("destiny", quietly = TRUE)){
	    stop("Package \"destiny\" needed for this function to work. Please install
	         it.", call. = FALSE)
	  }
	  pt_matrix <- destiny::as.matrix(pt_obj)
		pt_info <- data.frame(
			dpt_dm_1 = pt_obj@dm@eigenvectors[, 1],
			dpt_dm_2 = pt_obj@dm@eigenvectors[, 2],
			dpt_dm_3 = pt_obj@dm@eigenvectors[, 3],
			dpt_pt = pt_matrix[, 1],
			dpt_branch = as.character(pt_obj@branch[, 1])
			)
		rownames(pt_info) <- rownames(mtec@meta.data)
	}else {
		warning("Only transfer from monocle and dpt has been implemented. 
			Returning seruat object untouched")
		return(mtec)
	}
	cells <- rownames(mtec@meta.data)
	if (nrow(pt_info) == length(cells)){
		pt_info <- pt_info[match(cells, rownames(pt_info)),]
		mtec <- Seurat::AddMetaData(mtec, pt_info)
	} else {
		meta_data <- mtec@meta.data
		new_meta_data <- merge(meta_data, pt_info, by = "row.names", all.x = TRUE)
		rownames(new_meta_data) <- new_meta_data$Row.names
		new_meta_data$Row.names <- NULL
		mtec@meta.data <- new_meta_data
	}
	return(mtec)
}
