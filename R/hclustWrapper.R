#' Wrapper for hierarchical clustering
#'
#' This fuction is a wrapper for the hclust, dist, and cutree functions. It
#' outputs a string of a concatenated vector of 1's and 2's indexed by column.
#' This function is not meant to be run individually, but as a 'clustFunc'
#' argument for running `K2preproc()`, `runK2Taxonomer()`, and `K2tax()`.
#'
#' @param dataMatrix An P x N numeric matrix of data.
#' @param clustList Named list of objects to use for clustering procedure. Or
#' `NULL` to use Euclidean distance with Ward's method.
#' \itemize{
#'  \item{'distMetric'}{Distance method to use. See ?dist.}
#'  \item{'aggMethod'}{Agglomerative method to use. See ?hclust.}
#' }
#' @return A character string of concatenated 1's and 2's pertaining to the
#' cluster assignment of each column in dataMatrix.
#' @references
#'  \insertRef{reed_2020}{K2Taxonomer}
#' @keywords clustering
#' @export
#' @examples
#'
#' ## Read in sample ExpressionSet
#' require(Biobase)
#' data('sample.ExpressionSet')
#'
#' ## Extrat subset of data matrix
#' dm <- exprs(sample.ExpressionSet)[seq_len(50),]
#'
#' ## Create list of objects for clustering procedure.
#' cL <- list(distMetric='euclidean', aggMethod='ward.D2')
#'
#' ## Generate K=2 split with hierarchical clustering
#' hclustWrapper(dm, cL)
#'

hclustWrapper <- function(dataMatrix, clustList) {

    if (length(clustList) == 0) {
        clustList <- list(aggMethod="ward.D2", distMetric="euclidean")
    }

    dDist <- dist(t(dataMatrix), method=clustList$distMetric)
    dClust <- hclust(dDist, method=clustList$aggMethod)
    modVec <- as.character(cutree(dClust, k=2))
    mods <- paste(modVec, collapse="")
    return(mods)
}
