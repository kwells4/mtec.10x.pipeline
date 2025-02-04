% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_funcs.R
\name{make_plot_df}
\alias{make_plot_df}
\title{Plots a data frame that is easily placed into ggplot2}
\usage{
make_plot_df(seurat_object, y_val, x_val, col_by = NULL)
}
\arguments{
\item{seurat_object}{a Seurat object}

\item{y_val}{a gene (or any conitnuous value from metadata) to plot}

\item{x_val}{cluster (or any discrete value from metadata) to plot}

\item{col_by}{OPTIONAL what should be used to color the plots. Defaults to cluster}
}
\description{
This function allows you to take a seurat object and use many pieces from it to make
a data frame that easily interacts with ggplot2
}
\keyword{data}
\keyword{frame}
\keyword{ggplot2,}
