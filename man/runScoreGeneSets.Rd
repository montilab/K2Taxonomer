% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runScoreGeneSets.R
\name{runScoreGeneSets}
\alias{runScoreGeneSets}
\title{Enrichment scoring of gene sets on expression data}
\usage{
runScoreGeneSets(K2res, ScoreGeneSetMethod = NULL, useCors = NULL)
}
\arguments{
\item{ScoreGeneSetMethod}{Character. Method for gene set scoring. Use one of either 'GSVA' (default) or 'AUCELL'.}

\item{useCors}{Numeric. Number of cores to use for parallelizable processes.}
}
\value{
An object of class K2.
}
\description{
Performs single-sample enrichment scoring using GSVA or AUCell. When AUCell
is specified, output values reflect the AUC levels scaled by 100 and log2
transformed.
}
\references{
\insertRef{reed_2020}{K2Taxonomer}
 \insertRef{gsva}{K2Taxonomer}
}
\keyword{clustering}
