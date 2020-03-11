#' Wrapper for hierarchical clustering
#'
#' This fuction is a wrapper for the hclust, dist, and cutree functions. It outputs a string of a concatenated vector of 1's and 2's indexed by column.
#' @param dataMatrix An P x N numeric matrix of data
#' @param clustList List of objects to use for clustering procedure (Note used in this function).
#' @return A character string of concatenated 1's and 2's pertaining to the cluster assignment of each column in dataMatrix.
#' @keywords clustering
#' @export
#' @examples
#' hclust_ward(dataMatrix)
#' 

hclust_ward <- function(dataMatrix, clustList = NULL){
  dDist <- dist(t(dataMatrix))
  dClust <- hclust(dDist, method = "ward.D2")
  modVec <- as.character(cutree(dClust, k = 2))
  mods <- paste(modVec, collapse = "")
  return(mods)
}