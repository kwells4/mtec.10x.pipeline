#' Plots qc information
#' 
#' This function allows you to plot a qc information given a seurat
#' object. 
#' @param mtec a Seurat object
#' @param qc_vals OPTIONAL a list of qc information to plot. Defaults to
#' c("nGene", "nUMI", "percent.mito").
#' @keywords qc
#' @import Seurat
#' @export
#' @examples 
#' qc_plot(mTEC.10x.data::mtec_trace, qc_vals = c("nGene", "nUMI"))

qc_plot <- function(mtec, qc_vals = c("nGene", "nUMI", "percent_mito")){
	Seurat::VlnPlot(object = mtec, features.plot = qc_vals, nCol = 3)
}

#' Plots PCA plots, loadings, heatmap, and jackstraw
#' 
#' This function allows you to plot a PCA and run all tests to determine
#' which PCs to use for the remainder of the analysis. More information on
#' how to determine significant PCs in the Seurat tutorial. Briefly, I used
#' two methods. 1. plotting the heatmap to see which PCs are driven by genes
#' that make sense and have clear separation. 2. using jackstraw to determine
#' which PCs. have an enrichment to genes with low p-values. These p-values
#' are determined by random permutation of a small subset of the data and
#' construct a null distribution against which genes can be tested. This 
#' method is slow.
#' @param mtec a Seurat object
#' @param visualize_pcs OPTIONAL a list of PCs to plot on a PCA plot and/or
#' a PC loading plot. If plotting a PCA, the list may only be length 2.
#' Defaults to c(1, 2)
#' @param test_pcs OPTIONAL a range of PCs to test with a heatmap or jackstraw
#' plot. Defaults to 1:12
#' @param pc_loadings OPTIONAL if PC loadings should be visualized. Defaults
#' to TRUE.
#' @param pca_plot OPTIONAL if a PCA plot should be plotted. Defaults to TRUE
#' @param heatmap OPTIONAL if a heatmap should be plotted for each of the test
#' PCs. Defaults to TRUE
#' @param jackstraw OPTIONAL if jackstraw should be run to determine PC
#' significance. May take a long time. Defaults to FALSE.
#' @param elbow OPTIONAL if an elbow plot should be created. Defaults to TRUE
#' @param cells.use OPTIONAL the number of cells to use to create the heatmap.
#' Defaults to 500.
#' @param use.full OPTIONAL if all genes should be loaded (or only the HVGs).
#' Defaults to FALSE
#' @param num_replicates OPTIONAL the number of replicates to make the jackstraw
#' plot. Defaults to 100.
#' @keywords PCA
#' @import Seurat
#' @export
#' @examples 
#' PC_plots(mTEC.10x.data::mtec_trace)
#' \dontrun{
#' PC_plots(mTEC.10x.data::mtec_trace, visualize_pcs = c(2, 3),
#' test_pcs = 1:10, jackstraw = TRUE)
#' }

PC_plots <- function(mtec, visualize_pcs = c(1, 2), test_pcs = 1:12,
                     pc_loadings = TRUE, pca_plot = TRUE, heatmap = TRUE,
                     jackstraw = FALSE, elbow = TRUE, cells.use = 500,
                     use.full = FALSE, num_replicates = 100){
	if (pc_loadings){
		Seurat::VizPCA(object = mtec, pcs.use = visualize_pcs)
	}

	if (pca_plot){
		if (length(visualize_pcs) == 2){
			Seurat:: PCAPlot(object = mtec, dim.1 = visualize_pcs[1],
			                 dim.2 = visualize_pcs[2])
		} else {
			stop("This PCA is 2 dimensions. Please only use two values in the
				visualize_pcs list")
		}
	}

	if (heatmap){
		Seurat::PCHeatmap(object = mtec, pc.use = test_pcs, cells.use = cells.use,
				  do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
	}

	if (jackstraw){
		mtec <- Seurat::JackStraw(object = mtec, num.replicate = num_replicates)
		Seurat::JackStrawPlot(object = mtec, PCs = test_pcs)
	}

	if (elbow){
		Seurat::PCElbowPlot(object = mtec)
	}
}
