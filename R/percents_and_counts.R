#' Determines percent of cycling cells
#' 
#' This function determines the percent of cycling cells. The input must
#' be a seurat object with a meta_data column called "cycle_phase" and containing
#' values "G2M", "G1", "S". Returns a named list with the percent of cells that
#' are classified as G2M or S.
#' @param seurat_obj a Seurat object
#' @param subset_seurat OPTIONAL if the seurat object should be subset by some
#' parameter before determining percent of cycing. Default is FALSE.
#' @param subset_by OPTIONAL what meta data column should be used to subset
#' the data. Default is "exp"
#' @param subset_val OPTIONAL What value from the meta data column to use to subset
#' the seurat object. The value must be a value from the subset_by specified. Default
#' is "isoControlBeg"
#' @return a named list containing the percent of cycling cells and named by subset_val
#' If no data set was provided, the list will be named "all"
#' @keywords cell cycle, percent
#' @export
#' @examples 
#' percent_cycling_cells(mTEC.10x.data::mtec_trace)


percent_cycling_cells <- function(seurat_obj, subset_seurat = FALSE, subset_by = "exp",
                                  subset_val = "isoControlBeg"){
  if (subset_seurat) {
    # Subset the seurat_obj based on the subset_val
    cells_use <- rownames(seurat_obj@meta.data)[
      seurat_obj@meta.data[[subset_by]] == subset_val]
    new_seurat <- Seurat::SubsetData(seurat_obj, cells.use = cells_use)
  } else {
  	new_seurat <- seurat_obj
  	subset_val <- "all_cells"
  }

  # Determine the number of cycling cells
  cycling_cells <- table(new_seurat@meta.data$cycle_phase)
  if (!("S" %in% cycling_cells)) {
    cycling_cells["S"] = 0 
  }

  # Determine the percent of cells
  cycling_percent <- (cycling_cells["G2M"] +
    cycling_cells["S"])/nrow(new_seurat@meta.data)

  # Make the cycling percent a named list, named by the data set
  names(cycling_percent) <- subset_val

  # Return the percent as a named list
  return(cycling_percent)
}

#' Determines median umi count
#' 
#' This function determines the median umi value given a seurat object. It can
#' determine the umi value for an entier object or just one subset of it. For
#' this function, subset_by is best as batches as the output can tell you
#' differences in median umi count between samples. Helpful for doing downsampling
#' based on umi
#' @param seurat_obj a Seurat object
#' @param subset_seurat OPTIONAL if the seurat object should be subset by some
#' parameter before determining percent of cycing. Default is FALSE.
#' @param subset_by OPTIONAL what meta data column should be used to subset
#' the data. Default is "exp"
#' @param subset_val OPTIONAL What value from the meta data column to use to subset
#' the seurat object. The value must be a value from the subset_by specified. Default
#' is "isoControlBeg"
#' @return a value giving the number of unique molecules
#' @keywords umi
#' @export
#' @examples 
#' get_umi(mTEC.10x.data::mtec_trace)
#' get_umi(mTEC.10x.data::mtec_trace, subset_seurat = TRUE, subset_by = "stage",
#'	subset_val = "Aire_positive")

get_umi <- function(seurat_obj, subset_seurat = FALSE, subset_by = "exp",
  subset_val = "isoControlBeg"){
  if (subset_seurat){
    if (identical(names(seurat_obj@ident), rownames(seurat_obj@meta.data))){
      seurat_obj <- Seurat::SetAllIdent(seurat_obj, id = subset_by)
      seurat_obj <- Seurat::SubsetData(seurat_obj, ident.use = subset_val,
        subset.raw = TRUE)
    } else {
      stop("ident and meta.data slots not in the same order")
    }
  }
  cell_matrix <- as.matrix(seurat_obj@raw.data)
  umi <- median(colSums(cell_matrix))
  return(umi)
}

#' Returns information about gene sets
#' 
#' This function determines the number of genes or UMIs for genes in each gene set
#' expressed in each cell and also determines the percent of genes in each gene
#' set represented in the data set. This means the percent of genes from the gene
#' list that is expressed in at least one cell.
#' @param seurat_obj a Seurat object
#' @param gene_lists a named list of gene sets. Can be just one gene set or many.
#' @param downsample_UMI OPTIONAL if the raw data should be downsampled to a 
#' specific value. Default is FALSE
#' @param one_batch OPTIONAL if the seurat object should be subset. This is
#' especially helpful when comparing between samples within one seurat object.
#' Default it FALSE.
#' @param batch OPTIONAL if one_batch is true, what should be used to downsample?
#' Default is "exp", but must be a column of metadata
#' @param batch_name OPTIONAL the group to keep. This must be a value in the 
#' metadata column specified by batch. This is also the name of the returned
#' list. Default is "all_cells".
#' @param lowest_UMI OPTIONAL The value used to downsample the UMI. Default is NULL
#' @param count OPTIONAL what should be returned. Can be percents, counts of genes
#' or counts of UMIs. Can be a list or a character string with values "genes", "UMI",
#' and/or "percent"
#' @return a list with two parts. One part is the count of genes/UMI the second is
#' the percent of all genes expressed. Both parts contain lists of every gene set
#' named by the gene set
#' @keywords umi, gene, percent
#' @export

percents_and_counts <- function(seurat_obj, gene_lists, downsample_UMI = FALSE,
  one_batch = FALSE, batch = "exp", batch_name = "all_cells",
  lowest_UMI = NULL, count = "genes"){
  # If not looking at all sample in a seurat object, than subset to the desired
  # batch
  if (one_batch){
    if (identical(names(seurat_obj@ident), rownames(seurat_obj@meta.data))){
      seurat_obj <- Seurat::SetAllIdent(seurat_obj, id = batch)
      seurat_obj <- Seurat::SubsetData(seurat_obj, ident.use = batch_name,
        subset.raw = TRUE)
    } else {
      stop("ident and meta.data slots not in the same order")
    }
  }

  # Grab the raw.data slot from the seurat object
  cell_matrix <- as.matrix(seurat_obj@raw.data)

  # Downsample the UMI 
  if (downsample_UMI){
    # DropletUtils is required for this function
    if (!requireNamespace("DropletUtils", quietly = TRUE)){
      stop("Package \"DropletUtils\" needed for this function to work. Please install it.",
        call. = FALSE)
    }

    # Determine the median UMI for the dataset and determine factor based on the
    # given lowest UMI
    data_umi <- median(colSums(cell_matrix))
    if (is.null(lowest_UMI)) {
      stop("If downsampling, you must provide a value.
        This value will be the median number of UMI after downsampling")
    }
    factor <- lowest_UMI/data_umi
    
    # Use DropletUtils to downsample the raw matrix
    set.seed(0)
    cell_matrix <- DropletUtils::downsampleMatrix(cell_matrix, prop = factor)
  }
  
  return_list <- list()
  if ("genes" %in% count) {
    # Determine the number of genes in each gene set present in each
    # cell in the data set
    count_list <- lapply(names(gene_lists), function(x)
      gene_count_function(cell_matrix, gene_lists[[x]], x))
    count_df <- do.call(cbind, count_list)
    count_df$exp <- batch_name
    return_list$counts <- count_df
  }
  if ("UMI" %in% count) {
    # Determine the number of UMIs in each gene set present in each
    # cell in the data set
    umi_list <- lapply(names(gene_lists), function(x)
      umi_count_function(cell_matrix, gene_lists[[x]], x))
    umi_df <- do.call(cbind, umi_list)
    umi_df$exp <- batch_name
    return_list$umi <- umi_df
  }
  if ("percent" %in% count) {
    # Determine the percent of genes in each set expressed in ANY cell in the data set
    # ie percent of genes seen in at least one cell.
    gene_percent_list <- sapply(names(gene_lists), function(x)
      percent_list(cell_matrix, gene_lists[[x]], x))
    return_list$percents <- gene_percent_list
  }

  # Return both and name based on the batch 
  return_list <- list(return_list)
  names(return_list) <- batch_name
  return(return_list)
}

#' Returns the number of genes from a gene set in each cell
#' 
#' This function determines the number of genes for genes in a gene set
#' expressed in each cell. Called by percents_and_counts, but can
#' work on its own too
#' @param cell_matrix a matrix with cells as columns and genes as rows
#' @param gene_lists a list of genes
#' @param list name the name of the gene list
#' @return a named list of the number of genes expressed in every cell
#' @keywords gene count, gene set
#' @export

gene_count_function <- function(cell_matrix, gene_list, list_name){
  # Subset cell matrix to only be genes of interest
  gene_matrix <- cell_matrix[rownames(cell_matrix) %in% gene_list, ]
  # Counts number of genes expressed in each cell
  gene_count <- apply(gene_matrix, 2, function(x) sum(x > 0))
  gene_count <- data.frame(gene_count)
  # Return as a named list named by the gene list name
  names(gene_count) <- list_name
  return(gene_count)
}

#' Returns the number of UMIs from a gene set in each cell
#' 
#' This function determines the number of UMIs for genes in a gene set
#' expressed in each cell. Called by percents_and_counts, but can
#' work on its own too
#' @param cell_matrix a matrix with cells as columns and genes as rows
#' @param gene_lists a list of genes
#' @param list name the name of the gene list
#' @return a named list of the number of UMIs expressed in every cell
#' @keywords umi count, gene set
#' @export

umi_count_function <- function(cell_matrix, gene_list, list_name){
  # Subset cell matrix to only be genes of interest
  gene_matrix <- cell_matrix[rownames(cell_matrix) %in% gene_list, ]
  # Counts number of molecules expressed in each cell
  umi_count <- apply(gene_matrix, 2, function(x) sum(x))
  umi_count <- data.frame(umi_count)
  # Return as a named list named by the gene list name
  names(umi_count) <- list_name
  return(umi_count)
}

#' Returns the percent of genes expressed in a data set
#' 
#' This function determines the percent of genes in a gene set represented
#' in the data set. This means the percent of genes from the gene list that
#' is expressed in at least one cell. Called by percents_and_counts, but can
#' work on its own too
#' @param cell_matrix a matrix with cells as columns and genes as rows
#' @param gene_lists a list of genes
#' @param list name the name of the gene list
#' @return a named list of the percent of all genes expressed. The list is named by
#' the name of the gene set
#' @keywords percent, gene set
#' @export

percent_list <- function(cell_matrix, gene_list, gene_list_name){
  # Subset cell matrix to only be genes of interest
  gene_matrix <- cell_matrix[rownames(cell_matrix) %in% gene_list, ]
  # Counts number of cells expressing each gene.
  cell_count <- apply(gene_matrix, 1, function(x) sum(x > 0))
  # Counts number of genes expressed by at least 1 cell
  expr_genes <- cell_count[cell_count > 0]
  # Number of expressed genes
  n_expr_genes <- length(expr_genes)
  # percent of expressed genes
  percent <- n_expr_genes/length(gene_list)
  
  names(percent) <- gene_list_name
  # Returns the number
  return(percent)
}

#' Gets percents and counts from the output of percents_and_counts
#' 
#' This function is a relatively silly work around to get a list of
#' just the return from the percents or just the return from the counts
#' for all gene lists. Works best with sapply as shown in the example
#' @param percent_counts_list the return value from percents_and_counts
#' @param list_slot which slot in the list to pull
#' @param data_type OPTIONAL return "counts", "umi", or "percent". Must be one of the
#' three. Default is "counts"
#' @return a named list of the percent of all genes expressed for the gene set
#' listed in "list_slot"
#' @keywords percent, gene set, gene count, umi count
#' @export
#' @examples 
#' \dontrun{
#' percents <- sapply(names(percents_counts_all), function(x) 
#'  get_perc_count(percents_counts_all, x, percent = TRUE), USE.NAMES = TRUE)
#' }

get_perc_count <- function(percent_counts_list, list_slot, data_type = "counts"){
  percent_counts_one <- percent_counts_list[[list_slot]]
  if (data_type == "percent") {
    return_val <- percent_counts_one$percents
  } else if (data_type == "counts") {
    return_val <- percent_counts_one$counts
  } else if (data_type == "UMI") {
    return_val <- percents_counts_one$umi
  } else {
    stop("must bet either percent, count, or umi")
  }
  return(return_val)
}