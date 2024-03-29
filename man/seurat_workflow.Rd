% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seurat_workflow.R
\name{seurat_workflow}
\alias{seurat_workflow}
\title{Runs a typical Seurat workflow on a Seurat object (up to
dimensionality reduction and clustering).}
\usage{
seurat_workflow(
  seurat_obj,
  num_variable_features,
  resolution_param = 0.8,
  visualization_method = "umap",
  num_dims = 10,
  algorithm = "louvain"
)
}
\arguments{
\item{seurat_obj}{A Seurat object that will be analyzed.}

\item{num_variable_features}{The number of variable features to use in the
analysis.}

\item{resolution_param}{The resolution parameter to use when clustering.}

\item{visualization_method}{Either "umap" or "tsne".}

\item{num_dims}{The number of principal components to use.}

\item{algorithm}{The clustering algorithm to use, either "louvain" or
"leiden".}
}
\value{
A Seurat object containing the relevant analysis results.
}
\description{
Given a Seurat object, returns a new Seurat that has been
normalized, had variable features identified,
scaled, had principal components computed, had clusters identified, and had
tSNE and UMAP embeddings determined.
}
