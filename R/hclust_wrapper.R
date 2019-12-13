#' Wrapper for hierarchical clustering
#'
#' This fuction is a wrapper for the hclust, dist, and cutree functions. It outputs a string of a concatenated vector of 1's and 2's indexed by column.
#' @param dataMatrix An P x N numeric matrix of data
#' @param distMetric A character string. One of the "method" arguments in the dist() function
#' @param aggMethod A character string.  One of the "method" arguments in the hclust() function
#' @return A character string of concatenated 1's and 2's pertaining to the cluster assignment of each column in dataMatrix.
#' @keywords clustering
#' @export
#' @examples
#' hclust_wrapper(dataMatrix, distMetric = "euclidean", aggMethod = "ward.D2")
#' 

hclust_wrapper <- function(dataMatrix, distMetric = "euclidean", aggMethod = "ward.D2"){
  dDist <- dist(t(dataMatrix), method = distMetric)
  dClust <- hclust(dDist, method = aggMethod)
  modVec <- as.character(cutree(dClust, k = 2))
  mods <- paste(modVec, collapse = "")
  return(mods)
}