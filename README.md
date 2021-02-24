## K2Taxonomer

### Introduction

*K2Taxonomer* is an R package built around a "top-down" recursive partitioning framework to perform unsupervised learning of nested “taxonomy-like” subgroups from high-throughput -omics data. This framework was devised to be flexible to different data structures, enabling its applicability to analyze both bulk and single-cell data sets. In addition to implementing the algorithm, this package includes functionality to annotate estimated subgroups using gene- and pathway-level analyses.

The recursive partitioning approach utilized by `K2Taxonomer` presents advantages over conventional unsupervised approaches, including:

- Identification of robust partitions of a set of observations by aggregating ensembles of partition estimates from repeated perturbations of the data.
- Performing partition-specific feature selection, preventing the need to perform feature selection on the whole data set prior to running the algorithm.
- Tailoring of analyses to specific data structures through the use of different clustering algorithms for performing partition estimation.

The documentation of this package describes how `K2Taxonomer` can be applied to either bulk or single-cell gene expression data. For analyses of single-cell gene expression data `K2Taxonomer` is designed to characterize nested subgroups of previously identified cell types, such as those previously estimated by scRNAseq clustering analysis.

A preprint of the manuscript describing `K2Taxonomer` is publicly available [here](https://www.biorxiv.org/content/10.1101/2020.11.05.370197v1/).

### Documentation

#### See the GitHub pages site
https://montilab.github.io/K2Taxonomer/

### Requirements

- R (>= 4.0)

### Installation

#### Install from GitHub

You may install `K2Taxonomer` from GitHub directly using the `devtools` R package or clone the repository and download from source. Typical download time is around 1 minute.

```r
install.packages("devtools")
devtools::install_github("montilab/K2Taxonomer")
```

### Usage

Here we demonstrate the basic functionality of `K2Taxonomer`, which is described in more detail in the vignette, [Running K2Taxonomer](https://montilab.github.io/K2Taxonomer/articles/RunningK2Taxonomer.html).

An alternative workflow for running `K2Taxonomer` for subgrouping cell type labels using single-cell expression data is described in the vignette, [Running K2Taxonomer on single-cell RNA sequencing data](https://montilab.github.io/K2Taxonomer/articles/K2Taxonomer_singlecell.html).

#### Load packages and read in gene expression data

```r
## K2Taxonomer package
library(K2Taxonomer)

## For creating and manipulating ExpressionSet objects
library(Biobase)

## Read in ExpressionSet object
data(sample.ExpressionSet)
```
#### Initialize `K2` object

```r
K2res <- K2preproc(sample.ExpressionSet)
```

#### Run K2Taxonomer algorithm

```r
K2res <- K2tax(K2res,
               stabThresh=0.5)
```

#### Run differential analysis on all subgroups

```r
K2res <- runDGEmods(K2res)
```

#### Run enrichment analysis on toy gene sets

```r
genesetsMadeUp <- list(
  GS1=genes[1:50],
  GS2=genes[51:100],
  GS3=genes[101:150]
)

K2res <- runGSEmods(K2res,
                     genesets=genesetsMadeUp,
                     qthresh=0.1)
```

#### Run single-sample enrichment on toy gene sets with *GSVA*

```r
K2res <- runGSVAmods(K2res,
                      ssGSEAalg="gsva",
                      ssGSEAcores=1,
                      verbose=FALSE)
```

#### Run differential analysis on single-sample enrichment

```r
K2res <- runDSSEmods(K2res)
```

#### Create dashboard of results

```r
K2dashboard(K2res,
            analysis_name="Example",
            output_dir=".")
```
