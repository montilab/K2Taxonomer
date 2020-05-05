#' Wrapper for hierarchical clustering
#'
#' This fuction is a wrapper for the constrained Kmeans algorithm using 
#' `lcvqe` from the `conclust` package. This function will subset each
#' cohort down to that with the smallest number of observations.
#' 
#' @param dataMatrix An P x N numeric matrix of data
#' @param clustList List of objects to use for clustering procedure.
#' @return A character string of concatenated 1's and 2's pertaining to the cluster assignment of each column in dataMatrix.
#' @keywords clustering
#' @export
#' @import conclust
#' @examples
#' cKmeans_wrapper_subsample(dataMatrix, clusList = list())
#' 

# Create wrapper to subsample
cKmeans_wrapper_subsample <- function(dataMatrix, clustList) {
  
  eMatSub <- clustList$eMat[rownames(dataMatrix), clustList$labs %in% colnames(dataMatrix)]
  labsSub <- clustList$labs[clustList$labs %in% colnames(dataMatrix)]
  
  # Subsample the data
  minSize <- min(table(labsSub))
  sVec <- unlist(lapply(unique(labsSub), function(x) sample(which(labsSub == x), minSize)))
  eMatSub <- eMatSub[,sVec]
  labsSub <- labsSub[sVec]
  
  # Get constraints
  mustLink <- outer(labsSub, labsSub, "==")
  mustLink[upper.tri(mustLink, diag = TRUE)] <- FALSE
  mustLink <- which(mustLink, arr.ind = TRUE)
  
  # Cluster data
  dClust = factor(lcvqe(t(eMatSub), 
                        k = 2, 
                        mustLink = mustLink, 
                        cantLink = matrix(c(1, 1), 
                                          nrow = 1), maxIter = clustList$maxIter), levels = c(1, 2))
  
  # Get label-level clusters
  dMat <- as.matrix(table(dClust, labsSub))[,colnames(dataMatrix)]
  modVec <- apply(dMat, 2, which.max)
  mods <- paste(modVec, collapse = "")
  
  return(mods)
}
