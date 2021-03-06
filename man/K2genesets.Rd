% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\name{K2genesets}
\alias{K2genesets}
\alias{K2genesets<-}
\alias{K2genesets,K2-method}
\alias{K2genesets<-,K2-method}
\title{Named list of genesets used in K2Taxonomer run}
\usage{
K2genesets(K2res)

K2genesets(K2res) <- value

\S4method{K2genesets}{K2}(K2res)

\S4method{K2genesets}{K2}(K2res) <- value
}
\arguments{
\item{K2res}{K2 class object.}

\item{value}{Named list. A list with an element for each ge}
}
\value{
Named list. A list with an element for each gene set.
}
\description{
Retrieve or assign named list of genesets for K2Taxonomer run.
}
\examples{
data(K2res)
head(K2genesets(K2res))

}
