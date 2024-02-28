# callback (Calibrated Clustering via Knockoffs) <img src="man/figures/callback_logo.png" align="right" alt="" width="120"/>


[![R CMD check](https://github.com/lcrawlab/callback/actions/workflows/check-standard.yml/badge.svg)](https://github.com/lcrawlab/callback/actions/workflows/check-standard.yml)
[![Docker Image CI](https://github.com/lcrawlab/callback/actions/workflows/docker-image.yml/badge.svg)](https://github.com/lcrawlab/callback/actions/workflows/docker-image.yml)


## Introduction

Standard single-cell RNA-sequencing (scRNA-seq) pipelines nearly always include unsupervised clustering as a key step in identifying biologically distinct cell types. A follow-up step in these pipelines is to test for differential expression between the identified clusters. When algorithms over-cluster, downstream analyses will produce inflated P-values resulting in increased false discoveries. Here, we present `callback` (Calibrated Clustering via Knockoffs): a new method for protecting against over-clustering by controlling for the impact of double-dipping. Importantly, our approach can be applied to any clustering algorithm (implemented here are the Louvain, Leiden, K-means, and hierarchical clustering algorithms). `callback`provides state-of-the-art clustering performance and can rapidly analyze large-scale scRNA-seq studies, even on a personal laptop.


## Installation

You can install the lastest development version by using the [devtools](https://CRAN.R-project.org/package=devtools) library. To install this package with devtools, use this command:

```r
devtools::install_github("lcrawlab/callback")
```


## Tutorial

```r
library(Seurat)
library(SeuratData)

library(callback)


set.seed(123)

# load pbmc3k dataset
SeuratData::InstallData("pbmc3k")
data("pbmc3k")

pbmc3k <- UpdateSeuratObject(pbmc3k)

pbmc3k <- NormalizeData(pbmc3k)
pbmc3k <- FindVariableFeatures(pbmc3k)
pbmc3k <- ScaleData(pbmc3k)
pbmc3k <- RunPCA(pbmc3k)
pbmc3k <- FindNeighbors(pbmc3k)
pbmc3k <- RunUMAP(pbmc3k, dims = 1:10)

pbmc_default <- FindClusters(pbmc3k)
pbmc_callback <- FindClustersCallback(pbmc3k)

DimPlot(pbmc_default) + DimPlot(pbmc_callback)
```

## The Method


## Questions and Feedback
For questions or concerns with callback, please contact
[Alan DenAdel](mailto:alan_denadel@brown.edu).



