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

#### Load packages and read in gene expression data

```r
library(K2Taxonomer)
```

#### Required data sets

K2Taxonomer requires two data inputs

  - An object comprising expression and observation data. This must be one of three object classes: `ExpressionSet`, `Seurat`, or `SingleCellExperiment`.
  - An object comprising a named list of gene signatures
  
##### Expression and observational data

This example was written for a seurat object which includes the following

  - An "integrated" data slot which contains batch corrected scaled data used in Seurat clustering.
  - An "RNA" slot containing the un-integreated expression data
  - A column called "seurat_clusters", which contains the cluster labels.
  
##### Gene sets

These objects are simply a named list of vectors containing gene identifiers.
For example,

```r
GENESETS <- list(
  GS1 = c("LYZ", "AIF1", "S100A11", "FCER1G", "SAT1", "LST1", "DUSP1", "S100A4", "CTSS", "SERPINA1"),
  GS2 = c("STMN1", "MYBL2", "HIST1H4C", "RPLP0", "RPSA", "TYMS", "NUSAP1", "HMGB1", "LDHB", "C12orf75")
)
```

#### Initialize K2Taxonomer

```r
RNGkind("L'Ecuyer-CMRG")
set.seed(1)
K2res <- K2preproc(seu,
                   cohorts="seurat_clusters",
                   seuAssay = "integrated",
                   seuAssayDS = "RNA",
                   featMetric="F",
                   logCounts=TRUE,
                   clustFunc="cKmeansDownsampleSqrt",
                   useCors=8,
                   DGEmethod = "mast",
                   genesets = GENESETS,
                   ScoreGeneSetMethod = "AUCELL")
```

#### Run K2T algorithm

```r
K2res <- K2tax(K2res)
```

#### Run differential expression analysis to identify markers of each partition

```r
K2res <- runDGEmods(K2res)
```

#### Run gene set enrichment based on significantly differently expression genes

```r
K2res <- runFISHERmods(K2res)
```

#### Run gene set scoring with specified algorithm

ScoreGeneSetMethod from `K2preproc()`. This can be either "AUCELL" or "GSVA".

```r
K2res <- runScoreGeneSets(K2res)
```

#### Run Difference gene set scoring

```r
K2res <- runDSSEmods(K2res)
```

#### Create dashboard

```r
K2dashboard(K2res)
```

### Functions for results visualization

#### Plot dendrogram of K2tax() output

```r
plot(K2dendro(K2res))
```

##### Create interactive dendrogram

```r
K2visNetwork(K2res)
```

#### Create table of differential gene expression results

```r
DGEtable <- getDGETable(K2res)
```

##### Create interactive table of differential gene expression results

```r
getDGEInter(K2res, minDiff = 1, node = c("A", "B"))
```

#### Create table of enrichment results

```r
ENRtable <- getEnrichmentTable(K2res)
```

##### Create interactive table of enrichment results

```r
getEnrichmentInter(K2res, nodes = c("A", "D"))
```

#### Create plots of gene and enrichment scores per cohort at specific nodes

##### Genes

```r
plotGenePathway(K2res, feature = "MALAT1", node = "A")
```

##### Enrichment Scores

```r
plotGenePathway(K2res, feature = "GS1", node = "A", type = "gMat")
```

