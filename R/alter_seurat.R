#' Changes cluster names of seurat object
#' 
#' This function allows you to change the names of clusters in a seurat object.
#' This could be to number ordered by developmental or cell type names
#' @param mtec a Seurat object
#' @param new_ids OPTIONAL A list of new ids. Must be the same length as current
#' ids. Defaults to c(1,2,6,0,4,3,5,7,8) as used for the initial mtec object
#' @keywords seurat
#' @import dplyr
#' @export
#' @examples 
#' change_clus_ids(mTEC.10x.data::mtec_trace)
#' change_clus_ids(mTEC.10x.data::mtec_trace,
#' new_ids = c(8, 7, 6, 5, 4, 3, 2, 1, 0))

change_clus_ids <- function(mtec, new_ids = c(1, 2, 6, 0, 4, 3, 5, 7, 8)){
  current_ids <- levels(mtec@ident)
  
  # map new values onto correct cells
  mtec@ident <- plyr::mapvalues(x = mtec@ident, from = current_ids,
                                to = new_ids)
  
  # Also replace the res.0.6 column (which contains cluster info)
  mtec@meta.data$res.0.6 <- mtec@ident
  if (typeof(new_ids) == "double"){
	  new_ids <- sort(new_ids)	
  }
  
  # Remake factors to be correct
  mtec@ident <- factor(mtec@ident, levels = unique(new_ids))
  mtec@meta.data$res.0.6 <- factor(mtec@meta.data$res.0.6, levels = unique(new_ids))
	return(mtec)
}

#' Subsets the seurat object to selected clusters
#' 
#' This function allows you to generate an empty seurat object with only
#' the cells in clusters of interest. Then the processing pipeline can
#' be repeated with the different group of cells.
#' @param mtec a Seurat object
#' @param clusters A list clusters to subset (ex c(1, 2, 3))
#' @keywords seurat, clusters
#' @export
#' @examples 
#' subset_seurat(mTEC.10x.data::mtec_trace, clusters = c(4, 5, 6))

subset_seurat <- function(mtec, clusters){
  mtec_data <- mtec@raw.data
  cluster_info <- data.frame(cluster = mtec@ident)
  cluster_info$cells <- row.names(cluster_info)
  
  # Only keeps cells in desired clusters
  cluster_cells <- cluster_info[cluster_info$cluster %in% clusters, ]
  cluster_cells <- row.names(cluster_cells)
  
  # Grab raw data info for the chosen cells
  mtec_data_clus <- mtec_data[, cluster_cells]
  
  # Make new Seurat object
  mtec_clus <- Seurat::CreateSeuratObject(raw.data  = mtec_data_clus,
                                          min.cells = 3,
                                          min.genes = 200, 
                                          project   = "10X_unbaised_mTEC_RFP_GFP")
	return(mtec_clus)
}

#' Adds stage information to seurat object
#'
#' This function allows you to add stage information to a seurat object. 
#' It depends on a list input by the user that links cluster number to 
#' stage.
#' @param mtec_obj a Seurat object
#' @param stage_list a list of stages and their corresponding cluster numbers.
#' @keywords seurat, clusters
#' @export
#' @examples
#' stage_list <- c("0" = "cTEC", "1" = "Immature", "2" = "Immature",
#'                "3" = "Immature", "4" = "Intermediate",
#'                "5" = "Mature", "6" = "Late_mature",
#'                "7" = "Tuft")
#' mtec_trace <- set_stage(mTEC.10x.data::mtec_trace, stage_list)

set_stage <- function(mtec_obj, stage_list) {
  mtec_obj@meta.data$stage <- "unknown"
 for (i in names(stage_list)) {
   mtec_obj@meta.data$stage[mtec_obj@ident == i] <- stage_list[i]
 }
  return(mtec_obj)
}