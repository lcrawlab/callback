# PCKnockoffs

## Introduction

Analysis of single-cell RNA sequencing data frequently involves the detection of groups of similar cells via a clustering algorithm followed by differential gene expression analysis to detect marker genes that differ between those groups.
Due to the fact that the same data are used for both of these steps (referred to as data double-dipping), $p$-values computed via the usual tests can be poorly calibrated and even conservative procedures like the Bonferroni correction are not guarenteed control the number of false positives. 
PCKnockoffs, is a $ p$-value-free framework that controls the false discovery rate and the generalized family-wise error rate when carrying out clustering followed by differential gene expression analysis. The PCKnockoffs method only requires a small computational overhead for common single-cell analysis workflows and can identify when ''over-clustering'' has been performed on single cell data. The PCKnockoffs library is compatible with the [Seurat](https://satijalab.org/seurat/) library.

## The Method

## Installation

You can install the lastest development version by using the [devtools](https://CRAN.R-project.org/package=devtools) library. To install this package with devtools, use this command:

```r
devtools::install_github("lcrawlab/PCKnockoffs")
```


## Tutorials

## Questions and Feedback
For questions or concerns with the MAPIT functions, please contact
[Alan DenAdel](mailto:alan_denadel@brown.edu).

## References


