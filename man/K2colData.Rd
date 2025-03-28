% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\name{K2colData}
\alias{K2colData}
\alias{K2colData<-}
\alias{K2colData,K2-method}
\alias{K2colData<-,K2-method}
\title{Data frame of phenotypic information from K2 object}
\usage{
K2colData(K2res)

K2colData(K2res) <- value

\S4method{K2colData}{K2}(K2res)

\S4method{K2colData}{K2}(K2res) <- value
}
\arguments{
\item{K2res}{K2 class object.}

\item{value}{Data frame. Observation information from K2Taxonomer run.}
}
\value{
Data frame. Observation information from K2Taxonomer run.
}
\description{
Retrieve or assign observation information in K2Taxonomer run.
}
\examples{
data(K2res)
head(K2colData(K2res))

}
