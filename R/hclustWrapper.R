#' Wrapper for hierarchical clustering
#'
#' This fuction is a wrapper for the hclust, dist, and cutree functions. It
#' outputs a string of a concatenated vector of 1's and 2's indexed by column.
#' @param dataMatrix An P x N numeric matrix of data
#' @param clustList List of objects to use for clustering procedure.
#' @return A character string of concatenated 1's and 2's pertaining to the
#' cluster assignment of each column in dataMatrix.
#' @keywords clustering
#' @export
#'

hclustWrapper <- function(dataMatrix, clustList) {
    
    if (length(clustList) == 0) {
        clustList <- list(aggMethod = "ward.D2", distMetric = "euclidean")
    }
    
    dDist <- dist(t(dataMatrix), method = clustList$distMetric)
    dClust <- hclust(dDist, method = clustList$aggMethod)
    modVec <- as.character(cutree(dClust, k = 2))
    mods <- paste(modVec, collapse = "")
    return(mods)
}
