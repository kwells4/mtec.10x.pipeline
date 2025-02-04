% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pseudotime.R
\name{run_dpt}
\alias{run_dpt}
\title{Runs diffusion pseudotime (dpt) on cells}
\usage{
run_dpt(mtec, seed = 0)
}
\arguments{
\item{mtec}{a Seurat object}

\item{seed}{OPTIONAL seed for reporducibility. Defaults to 0.}
}
\description{
This function allows you to run dpt on your samples. It also will add the
dpt output into the seurat object. While this can be run on a local, I
recommend only running this on a cluster. If you only have a just use the
mtec_trace object as it already includes the data from dpt.
}
\examples{
\dontrun{
run_dpt(mTEC.10x.data::mtec_trace)
}
}
\keyword{dpt}
