#' Wrapper for constrained K-means
#'
#' This fuction is a wrapper for the constrained Kmeans algorithm using
#' `lcvqe` from the `conclust` package.
#'
#' @param dataMatrix An P x N numeric matrix of data
#' @param clustList List of objects to use for clustering procedure.
#' @return A character string of concatenated 1's and 2's pertaining to the
#' cluster assignment of each column in dataMatrix.
#' @references
#'  \insertRef{reed_2020}{K2Taxonomer}
#'  \insertRef{cKm}{K2Taxonomer}
#' @keywords clustering
#' @export
#' @import conclust
#'

cKmeansWrapper <- function(dataMatrix, clustList) {
    
    clustList$labs <- as.character(clustList$labs)
    
    eMatSub <- clustList$eMat[rownames(dataMatrix), clustList$labs %in% colnames(dataMatrix)]
    labsSub <- clustList$labs[clustList$labs %in% colnames(dataMatrix)]
    
    ## Get constraints
    mustLink <- outer(labsSub, labsSub, "==")
    mustLink[upper.tri(mustLink, diag = TRUE)] <- FALSE
    mustLink <- which(mustLink, arr.ind = TRUE)
    
    ## Cluster data
    dClust = factor(lcvqe(t(eMatSub), k = 2, mustLink = mustLink, cantLink = matrix(c(1, 
        1), nrow = 1), maxIter = clustList$maxIter), levels = c(1, 2))
    
    ## Get label-level clusters
    dMat <- as.matrix(table(dClust, labsSub))[, colnames(dataMatrix)]
    modVec <- apply(dMat, 2, which.max)
    mods <- paste(modVec, collapse = "")
    
    return(mods)
}
