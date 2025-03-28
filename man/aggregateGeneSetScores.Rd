% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregateGeneSetScores.R
\name{aggregateGeneSetScores}
\alias{aggregateGeneSetScores}
\title{Take difference of two paired bi-directional gene set scores}
\usage{
aggregateGeneSetScores(K2res, aggList)
}
\arguments{
\item{K2res}{A K2 class object.}

\item{aggList}{A named list where each item is a character vector of length, 
2, comprising the name of the 'up' gene set, and the name of the 'down'
gene set.}
}
\value{
An object of class K2.
}
\description{
Replaces gene set results from paired up- and down- gene sets with the difference
of the up-regulated genes and down-regulated genes
}
\examples{
## Read in ExpressionSet object
library(Biobase)
data(sample.ExpressionSet)

## Pre-process and create K2 object
K2res <- K2preproc(sample.ExpressionSet)

## Run K2 Taxonomer algorithm
K2res <- K2tax(K2res,
            stabThresh=0.5)

## Run differential analysis on each partition
K2res <- runDGEmods(K2res)

## Create dummy set of gene sets
DGEtable <- getDGETable(K2res)
genes <- unique(DGEtable$gene)
genesetsMadeUp <- list(
    GS1=genes[1:50],
    GS2=genes[51:100],
    GS3=genes[101:150])

## Run gene set hyperenrichment
K2res <- runGSEmods(K2res,
                genesets=genesetsMadeUp,
                qthresh=0.1)

## Run GSVA on genesets
K2res <- runGSVAmods(K2res,
                ssGSEAalg='gsva',
                ssGSEAcores=1,
                verbose=FALSE)

## Aggregate paired gene sets
aggList <- list(c('GS12', 'GS1', 'GS2'))
K2res <- aggregateGeneSetscores(K2resaggList, K2res)

}
\references{
\insertRef{reed_2020}{K2Taxonomer}
   \insertRef{gsva}{K2Taxonomer}
}
\keyword{clustering}
