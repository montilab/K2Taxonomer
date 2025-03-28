% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R
\name{K2eMatDS}
\alias{K2eMatDS}
\alias{K2eMatDS<-}
\alias{K2eMatDS,K2-method}
\alias{K2eMatDS<-,K2-method}
\title{Expression matrix object used for differential expression}
\usage{
K2eMatDS(K2res)

K2eMatDS(K2res) <- value

\S4method{K2eMatDS}{K2}(K2res)

\S4method{K2eMatDS}{K2}(K2res) <- value
}
\arguments{
\item{K2res}{K2 class object.}

\item{value}{Matrix.}
}
\value{
Matrix.
}
\description{
Retrieve or assign Matrix object.
}
\examples{
data(K2res)
K2eMatDS(K2res)

}
