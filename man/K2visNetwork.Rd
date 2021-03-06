% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/K2visNetwork.R
\name{K2visNetwork}
\alias{K2visNetwork}
\title{Interactive K2 dendrogram}
\usage{
K2visNetwork(K2res)
}
\arguments{
\item{K2res}{A list object. The output of runK2tax().}
}
\value{
An interactive dendrogram created by `visNetwork::visNetwork()`.
}
\description{
Create an interactive dendrogram of the K2 Taxonomer results
}
\examples{

## Read in K2 Taxonomer results
data(K2res)

## Generate interactive dendrogram
K2visNetwork(K2res)

}
\references{
\insertRef{reed_2020}{K2Taxonomer}
}
