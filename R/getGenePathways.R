#' Map features (genes) to feature list (genesets)
#'
#' This function will create a named vector of collapsed names of named
#' genesets list, separated by '; '
#' @param genesets A named list of features in row names of dataMatrix.
#' @keywords clustering
#' @return A named vector
#' @export
#' @examples
#' ## Read in ExpressionSet object
#' library(Biobase)
#' data(sample.ExpressionSet)
#' 
#' ## Create dummy set of gene sets
#' genes <- rownames(sample.ExpressionSet)
#' genesetsMadeUp <- list(
#'     GS1 = genes[1:50],
#'     GS2 = genes[51:100],
#'     GS3 = genes[101:150]
#' )
#' 
#' head(getGenePathways(genesetsMadeUp))
#'

getGenePathways <- function(genesets) {
    
    ## Get unique gene and geneset ids
    geneNames <- unique(unlist(genesets))
    genesetNames <- names(genesets)
    
    genesetMat <- do.call(cbind, lapply(genesets, function(x) geneNames %in% x))
    
    ## Get set of pathways in which each gene resides
    gene2Pathway <- apply(genesetMat, 1, function(x) paste(genesetNames[x], collapse = "; "))
    names(gene2Pathway) <- geneNames
    
    return(gene2Pathway)
}
