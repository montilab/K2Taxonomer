## K2Taxonomer

### Introduction

*K2Taxonomer* is an R package built around a "top-down" recursive partitioning framework to perform unsupervised learning of nested “taxonomy-like” subgroups from high-throughput -omics data. This framework was devised to be flexibly applicable to different data structures, supporting the analysis of both bulk and single-cell data sets. In addition to implementing the algorithm, the package includes functionality to annotate estimated subgroups using gene- and pathway-level analyses.

The recursive partitioning approach utilized by `K2Taxonomer` presents advantages over conventional unsupervised approaches, including:

- Identification of robust partitions of a set of observations by aggregating ensembles of partition estimates from repeated perturbations of the data.
- Partition-specific feature selection, preventing the need to perform feature selection on the whole data set prior to running the algorithm.
- Tailoring of analyses to specific data structures through the use of different clustering algorithms for partition estimation.

The package documentation describes applications of `K2Taxonomer` to both single-cell and bulk gene expression data. For analyses of single-cell gene expression data `K2Taxonomer` is designed to characterize nested subgroups of previously identified cell types, such as those previously estimated by scRNAseq clustering analysis.

### Cite

Reed, Eric R, and Stefano Monti. “Multi-Resolution Characterization of Molecular Taxonomies in Bulk and Single-Cell Transcriptomics Data.” _Nucleic Acids Research_ 49, no. 17 (July 6, 2021): e98. https://doi.org/10.1093/nar/gkab552.

### Documentation

Articles describing `K2Taxonomer` workflows can be found on the package's [GitHub Page](https://montilab.github.io/K2Taxonomer/).

### Requirements

- R (>= 4.0)

### Installation

You may install `K2Taxonomer` from GitHub directly using the `devtools` R package or clone the repository and download from source.

```r
install.packages("devtools")
devtools::install_github("montilab/K2Taxonomer")
```



