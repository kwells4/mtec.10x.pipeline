#' Adds number of reads to seurat object
#' 
#' This function allows you to add number of reads to seurat metadata
#' @param mtec a seurat object
#' @param reads_df the data frame containing number of reads per cell
#' @keywords meta.data, nReads
#'
#' @import Seurat
#' @export

add_nReads <- function(mtec, reads_df){

	# only keep cells in seurat object
	rownames(reads_df) <- regmatches(rownames(reads_df),
	                                 regexpr("[ATCG]+",
	                                         rownames(reads_df)))
	colnames(reads_df) <- "nReads"

	# add number of reads to seurat object
	mtec <- Seurat::AddMetaData(object = mtec, metadata = reads_df,
	                            col.name = "nReads")
	return(mtec)
}

#' Adds mito percent to seurat object
#' 
#' This function allows you to add mito percent to seurat meta data
#' @param mtec a seurat object
#' @param mito_pattern OPTIONAL the pattern of mito reads. Defaults to mouse
#' (^mt-)
#' @keywords meta.data, perc.mito
#'
#' @import Seurat
#' @export
#' @examples 
#' add_perc_mito(mTEC.10x.data::mtec_trace_empty)

add_perc_mito <- function(mtec, mito_pattern = "^mt-"){
	mito_genes <- grep(pattern = mito_pattern, x = rownames(x = mtec@data),
	                   value = TRUE)
	percent_mito <- Matrix::colSums(
	  mtec@raw.data[mito_genes, ])/Matrix::colSums(mtec@raw.data)
	mtec <- Seurat::AddMetaData(object = mtec, metadata = percent_mito,
	                            col.name = "percent_mito")
	return(mtec)
}

#' Initially processes cells in seurat object
#' 
#' This function allows you to filter cells, normalize data,
#' find variable genes,and scale data
#' @param mtec a seurat object
#' @param filter_by OPTIONAL what to use to filter cells. Must be in object
#' meta.data. Can be a list. Defaults to c("nGene", "percent.mito")
#' @param low_thresholds OPTIONAL low threshold for each factor used to filter
#' cells. Can be a list. Defaults to c(200, -Inf)
#' @param high_thresholds OPTIONAL high threshold for each factor used to filter
#' cells. Can be a list. Defaults to c(7500, 0.1)
#' @param scale_factor OPTIONAL how to scale the data to normalize. Defaults to
#' 10000.
#' @param filter_cells OPTIONAL if cells should be filtered. Defaults to TRUE
#' @param normalize_data OPTIONAL if values should be normalized. Defaults to
#' TRUE.
#' @param find_var_genes OPTIONAL if variable genes should be determined.
#' Defaults to TRUE.
#' @param scale_data OPTIONAL if the data should be scaled. Defaults to TRUE
#' @param PCA OPTIONAL if PCs should be found. Defaults to TRUE.
#' @keywords meta.data, perc.mito
#'
#' @import Seurat
#' @export
#' @examples 
#' process_cells(mTEC.10x.data::mtec_trace_empty, filter_by = c("nGene"),
#' high_thresholds = c(7500), low_thresholds = c(200))

process_cells <- function(mtec, filter_by = c("nGene", "percent_mito"),
						  low_thresholds = c(200, -Inf), high_thresholds = c(7500, 0.1),
						  scale_factor = 10000, filter_cells = TRUE,
						  normalize_data = TRUE, find_var_genes = TRUE,
						  scale_data = TRUE, PCA = TRUE){
  
	if (filter_cells){
		mtec <- Seurat::FilterCells(object = mtec, subset.names = filter_by,
							low.thresholds = low_thresholds,
							high.thresholds = high_thresholds)
	}
	
	if (normalize_data){
		mtec <- Seurat::NormalizeData(object = mtec,
		                              normalization.method= "LogNormalize",
		                              scale.factor = scale_factor)
	}

	if (find_var_genes){
		mtec <- Seurat::FindVariableGenes(object = mtec, mean.function = ExpMean,
								  dispersion.function = LogVMR,
								  x.low.cutoff = 0.0125, x.high.cutoff = 3,
								  y.cutoff = 0.5)
	}

	if (scale_data){
		mtec <- Seurat::ScaleData(object = mtec)
	}

	if (PCA){
		mtec <- Seurat::RunPCA(object = mtec, pc.genes = mtec@var.genes,
		                       do.print = TRUE, pcs.print = 1:5, 
		                       genes.print = 5, seed.use = 42)

		mtec <- Seurat::ProjectPCA(object = mtec, do.print = FALSE)
	}

	return(mtec)
}


#' Groups cells with clustering and/or tSNE
#' 
#' This function groups cells by clustering and or by tSNE given a Seurat
#' object
#' @param mtec a Seurat object
#' @param dims_use OPTIONAL which PC dimensions to use. Should be based on the
#' heatmap, jackstraw plot, and elbow plot. Defaults to 1:10
#' @param random_seed OPTIONAL the random seed to allow reproducibility.
#' Defaults to 0.
#' @param resolution OPTIONAL the resolution for the clustering. Best between
#' 0.6 and 1.2 look at Seurat tutorial for more info. Defaults to 0.6
#' @param cluster OPTIONAL if the cells should be clustered. Defaults to TRUE.
#' @param tSNE OPTIONAL if a tSNE should be calculated. Defaults to TRUE.
#' @param UMAP OPTIONAL if a umap should be calculated. Defaults to TRUE
#' @keywords tSNE, clusters, umap
#' @import Seurat
#' @export
#' @examples 
#' group_cells(mTEC.10x.data::mtec_trace)
#' group_cells(mTEC.10x.data::mtec_trace, dims_use = 1:5, tSNE = FALSE)


group_cells <- function(mtec, dims_use = 1:10, random_seed = 0,
						resolution = 0.6, cluster = TRUE, tSNE = TRUE,
						UMAP = TRUE){
	if (cluster){
		mtec <- Seurat::FindClusters(object = mtec, reduction.type = "pca",
							 dims.use = dims_use, resolution = resolution,
							 print.output = 0, save.SNN = TRUE,
							 random.seed = random_seed)
	}

	if (tSNE){
		mtec <- Seurat::RunTSNE(object = mtec, dims.use = dims_use, do.fast = TRUE,
						seed.use = random_seed)
	}

	if (UMAP) {
		mtec <- Seurat::RunUMAP(object = mtec, dims.use = dims_use, reduction.use = "pca",
								sed.usee = random_seed)
	}
	return(mtec)
}
