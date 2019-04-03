#' Runs cell cycle analysis on cells
#' 
#' This function allows you to run cell cycle analysis on the population of
#' cells using cyclone from scran.
#' @param mtec a Seurat object
#' @param conversion_file a file to use to convert gene names to ensembl ids.
#' @param seed OPTIONAL seed for reproducibility. Defaults to 0
#' @param add_to_obj OPTIONAL if the data should be added to the Seurat object.
#' Defaults to TRUE
#' @keywords scran
#' @export
#' @examples 
#' \dontrun{
#' run_cyclone(mTEC.10x.data::mtec_trace)
#' run_cyclone(mTEC.10x.data::mtec_trace, seed = 10)
#' }

run_cyclone <- function(mtec, conversion_file, seed = 0,
                        add_to_obj = TRUE){

  # Throw an error if scran isn't installed
  if (!requireNamespace("scran", quietly = TRUE)) {
    stop("Package \"scran\" needed for this function to work. Please install it.",
         call. = FALSE)
    }
  mtec_data <- mtec@raw.data

  # Determines genes to use for mouse
  mm_pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds",
                                  package="scran"))

  colnames(conversion_file) <- c("Ens_id", "gene_id")

  expr_genes <- as.data.frame(rownames(mtec_data))

  colnames(expr_genes) <- "gene_id"

  mtec_genes <- merge(x = expr_genes, y = conversion_file,
                      all.x = TRUE, all.y = FALSE, by = "gene_id")

  mtec_genes <- mtec_genes[!duplicated(mtec_genes$gene_id), ]

  mtec_genes <- mtec_genes[order(match(mtec_genes$gene_id,
                                       expr_genes$gene_id)), ]

  rownames(mtec_data) <- mtec_genes$Ens_id

  # Set seed for reproducibility
  set.seed(0)
  
  # Run cyclone
  assigned <- scran::cyclone(as.matrix(mtec_data), pairs=mm_pairs)

  # Determine output of cyclone
  cell_cycle_stage <- as.data.frame(assigned$normalized.scores)
  rownames(cell_cycle_stage) <- colnames(mtec_data)
  cell_cycle_stage$cycle_phase <- assigned$phases

  # Add cyclone output to object if desired. Calls add_cell_cycle function
  if (add_to_obj){
    mtec <- add_cell_cycle(mtec, cell_cycle_stage)
		return(mtec)
  
  # Otherwise return the cyclone output
  } else {
    return(cell_cycle_stage)
  }	
}

#' Adds cell cycle info to the meta data
#' 
#' This function will add the cell cycle information created by run_cyclone to
#' the seurat object. Called directly by run_cyclone and works best that way.
#' @param mtec a Seurat object
#' @param cell_cycle_stage the cell cycle stage output file from cyclone. Can be
#' cell cycle stage determined with other methods. Requries that cell names are
#' row names of the data frame and that these names are identical to the names
#' in the seurat object.
#' @keywords seurat, cell cycle
#' @export

add_cell_cycle <- function(mtec, cell_cycle_stage){

  meta_data <- mtec@meta.data
  
  # Add cell cycle info to meta data
  meta_data <- merge(cell_cycle_stage, meta_data, by = "row.names",
                     all.x = FALSE, all.y = TRUE)
  rownames(meta_data) <- meta_data$Row.names
  meta_data$Row.names <- NULL

  # Make sure the order of the new meta.data is consistent with the package
  if(!identical(rownames(meta_data), colnames(mtec@data))) {
    print("reordering meta data rows")
    meta_data <- meta_data[match(colnames(mtec@data),
                                 rownames(meta_data)), ]
  }
  # Replace old meta data with new meta data containing cell cycle
  mtec@meta.data <- meta_data
  return(mtec)
}
