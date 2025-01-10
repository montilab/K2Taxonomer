#' Wrapper for constrained K-means on data subsampled to the square root of the 
#' number of observations in each cohort.
#'
#' This function is a wrapper for the constrained Kmeans algorithm using
#' lcvqe() from the conclust package. This function will subset each
#' cohort down to that with the smallest number of observations.This
#' function is not meant to be run individually, but as a 'clustFunc'
#' argument for running K2tax().
#' @param labels Vector of cohort values
#' @param features List of features (genes) to include
#' @return A character string of concatenated 1's and 2's pertaining to the
#' cohort assignments.
#' @references
#'  \insertRef{reed_2020}{K2Taxonomer}
#'  \insertRef{cKm}{K2Taxonomer}
#' @inheritParams K2tax
#' @export
#' @import conclust

## Wrapper to subsample
cKmeansDownsampleSqrt <- function(labels, features, K2res) {
  
  if("maxIter" %in% names(K2meta(K2res)$clustList)) {
    MI <- K2meta(K2res)$clustList$maxIter
  } else {
    MI <- 25
  }
  
  labs <- as.character(K2colData(K2res)[, K2meta(K2res)$cohorts])
  obsKeep <- labs %in% labels
  
  labsSub <- labs[obsKeep]
  eMatSub <- K2eMat(K2res)[features, obsKeep]
  
  ## Subsample the data
  sVec <- unlist(lapply(unique(labsSub), function(x) {
    wsamps <- which(labsSub == x)
    sample(wsamps, sqrt(length(wsamps)))
  }))
  eMatSub <- eMatSub[, sVec]
  labsSub <- labsSub[sVec]
  
  ## Set constraints
  mustLink <- outer(labsSub, labsSub, "==")
  mustLink[upper.tri(mustLink, diag=TRUE)] <- FALSE
  mustLink <- which(mustLink, arr.ind=TRUE)
  clink <- sample(nrow(eMatSub), 1)
  
  ## Cluster data
  dClust <- factor(lcvqe(t(eMatSub), k=2, mustLink=mustLink,
                      cantLink=matrix(c(clink, clink), nrow=1), maxIter=MI),
                levels=c(1, 2))
  
  ## Get label-level clusters
  dMat <- as.matrix(table(dClust, labsSub))[, labels]
  
  modVec <- apply(dMat, 2, which.max)
  mods <- paste(modVec, collapse="")
  
  return(mods)
}
