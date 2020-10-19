## K2Taxonomer

### Documentation

#### See the GitHub pages site
https://montilab.github.io/K2Taxonomer/

### Requirements

- R (>= 3.5)

### Installation

#### Install from GitHub

```r
devtools::install_github("montilab/K2Taxonomer")
```

#### Clone GitHub and install from source

```sh
## In console
git clone https://github.com/montilab/K2Taxonomer
```

```r
## In R
install.packages("path/to/K2Taxonomer", repos=NULL, type="source")
```

### Usage

Full example available in vignette: ["Running K2Taxonomer"](https://montilab.github.io/K2Taxonomer/articles/RunningK2Taxonomer.html).

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
               stabThresh = 0.5)
```

#### Run differential analysis on all subgroups

```r
K2res <- runDGEmods(K2res)
```

#### Run enrichment analysis on toy gene sets

```r
genesetsMadeUp <- list(
  GS1 = genes[1:50],
  GS2 = genes[51:100],
  GS3 = genes[101:150]
)

K2res <- runGSEmods(K2res, 
                     genesets = genesetsMadeUp,
                     qthresh = 0.1)
```

#### Run single-sample enrichment on toy gene sets with *GSVA*

```r
K2res <- runGSVAmods(K2res, 
                      ssGSEAalg = "gsva",
                      ssGSEAcores = 1,
                      verbose = FALSE)
```

#### Run differential analysis on single-sample enrichment

```r
K2res <- runDSSEmods(K2res)
```

#### Create dashboard of results

```r
K2dashboard(K2res,
            analysis_name = "Example",
            output_dir = ".")
```
