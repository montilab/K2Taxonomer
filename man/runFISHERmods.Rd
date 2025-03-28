% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runFISHERmods.R
\name{runFISHERmods}
\alias{runFISHERmods}
\title{Perform Fisher test overrepresentation analysis of each subgroup}
\usage{
runFISHERmods(
  K2res,
  genesets = NULL,
  qthresh = NULL,
  cthresh = NULL,
  ntotal = NULL
)
}
\arguments{
\item{genesets}{Named list. Feature sets to be includes in enrichment-based analyses.}

\item{qthresh}{Numeric. A value between 0 and 1 indicating the FDR cuttoff to define feature sets.}

\item{cthresh}{Numeric. A positive value for the coefficient cuttoff to define feature sets.}

\item{ntotal}{Numeric. A positive value to use as the background feature count. 20000 by default.}
}
\value{
An object of class K2.
}
\description{
Adds overrepresentation testing results using the output of runDGEmods().
}
\references{
\insertRef{reed_2020}{K2Taxonomer}
   \insertRef{bh}{K2Taxonomer}
}
\keyword{clustering}
