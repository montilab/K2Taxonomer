#' Wrapper for hierarchical clustering
#'
#' This fuction is a wrapper for the hclust, dist, and cutree functions. It
#' outputs a string of a concatenated vector of 1's and 2's indexed by column.
#' This function is not meant to be run individually, but as a 'clustFunc'
#' argument for running K2tax().
#'
#' @param labels Vector of cohort values
#' @param features List of features (genes) to include
#' @return A character string of concatenated 1's and 2's pertaining to the
#' cohort assignments.
#' @references
#'  \insertRef{reed_2020}{K2Taxonomer}
#' @keywords clustering
#' @inheritParams K2tax
#' @export

hclustWrapper <- function(labels, features, K2res) {

    if ("aggMethod" %in% names(K2meta(K2res)$clustList)) {
      aggMethod <- K2meta(K2res)$clustList$aggMethod
    } else {
      aggMethod <- "ward.D2"
    }
  
    if ("distMetric" %in% names(K2meta(K2res)$clustList)) {
      distMetric <- K2meta(K2res)$clustList$distMetric
    } else {
      distMetric <- "euclidean"
    }
  
    eMatSub <- K2data(K2res)[features, labels]

    dDist <- dist(t(eMatSub), method=distMetric)
    dClust <- hclust(dDist, method=aggMethod)
    modVec <- as.character(cutree(dClust, k=2))
    mods <- paste(modVec, collapse="")
    return(mods)
}
