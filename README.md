# callback (Calibrated Clustering via Knockoffs) <img src="man/figures/callback_logo.png" align="right" alt="" width="120"/>


[![R CMD check](https://github.com/lcrawlab/PCKnockoffs/actions/workflows/check-standard.yml/badge.svg)](https://github.com/lcrawlab/PCKnockoffs/actions/workflows/check-standard.yml)

[![Docker Image CI](https://github.com/lcrawlab/PCKnockoffs/actions/workflows/docker-image.yml/badge.svg)](https://github.com/lcrawlab/PCKnockoffs/actions/workflows/docker-image.yml)


## Introduction

Analysis of single-cell RNA sequencing data frequently involves the detection of groups of similar cells via a clustering algorithm followed by differential gene expression analysis to detect marker genes that differ between those groups.
Due to the fact that the same data are used for both of these steps (referred to as data double-dipping), P-values computed via the usual tests can be poorly calibrated and even conservative procedures like the Bonferroni correction are not guarenteed control the number of false positives. 
KCC, is a P-value-free framework that controls the false discovery rate and the generalized family-wise error rate when carrying out clustering followed by differential gene expression analysis. The KCC method only requires a small computational overhead for common single-cell analysis workflows and can identify clusters without ''over-clustering'' the data. The KCC library is compatible with the [Seurat](https://satijalab.org/seurat/) library.

## The Method

## Installation

You can install the lastest development version by using the [devtools](https://CRAN.R-project.org/package=devtools) library. To install this package with devtools, use this command:

```r
devtools::install_github("lcrawlab/KCC")
```


## Tutorial

```r
library(Seurat)
library(kcc)

pbmc.data <- Read10X(data.dir = "")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)
pbmc <- FindNeighbors(pbmc)
pbmc <- RunUMAP(pbmc, dims = 1:10)

pbmc_default <- FindClusters(pbmc)
pbmc_kcc <- FindClustersKCC(pbmc)

DimPlot(pbmc_default) + DimPlot(pbmc_kcc)
```

## Questions and Feedback
For questions or concerns with KCC, please contact
[Alan DenAdel](mailto:alan_denadel@brown.edu).

## References

If you use `KCC` in your work, please cite the `KCC` publication as follows:

> **paper title**
>
> Alan DenAdel and Lorin Crawford
>
>_journal_ data doi: []().

