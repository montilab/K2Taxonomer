% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getGenePathways.R
\name{getGenePathways}
\alias{getGenePathways}
\title{Map features (genes) to feature list (genesets)}
\usage{
getGenePathways(genesets)
}
\arguments{
\item{genesets}{A named list of features in row names of dataMatrix.}
}
\value{
A named vector
}
\description{
This function will create a named vector of collapsed names of named
genesets list, separated by '; '
}
\examples{
## Read in ExpressionSet object
library(Biobase)
data(sample.ExpressionSet)

## Create dummy set of gene sets
genes <- rownames(sample.ExpressionSet)
genesetsMadeUp <- list(
    GS1=genes[1:50],
    GS2=genes[51:100],
    GS3=genes[101:150])

head(getGenePathways(genesetsMadeUp))

}
\keyword{clustering}
