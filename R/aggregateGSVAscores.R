#' Take difference of two paired GSVA scores
#'
#' Replaces GSVA results from paired up- and down- gene sets with the difference
#' of the up-regulated genes and down-regulated genes
#' @param aggList A list where each item is a vector of 3 items:
#' the new name, the name of the 'up' gene set, and the name of the 'down'
#' gene set.
#' @param K2res An object of class K2. The output of runDGEmods().
#' @return An object of class K2.
#' @references
#'  \insertRef{reed_2020}{K2Taxonomer}
#'  \insertRef{gsva}{K2Taxonomer}
#' @keywords clustering
#' @export
#' @import GSVA
#' @import Biobase
#' @examples
#' ## Read in ExpressionSet object
#' library(Biobase)
#' data(sample.ExpressionSet)
#' 
#' ## Pre-process and create K2 object
#' K2res <- K2preproc(sample.ExpressionSet)
#' 
#' ## Run K2 Taxonomer algorithm
#' K2res <- K2tax(K2res,
#'                stabThresh = 0.5)
#' 
#' ## Run differential analysis on each partition
#' K2res <- runDGEmods(K2res)
#' 
#' ## Create dummy set of gene sets
#' DGEtable <- getDGETable(K2res)
#' genes <- unique(DGEtable$gene)
#' genesetsMadeUp <- list(
#'     GS1 = genes[1:50],
#'     GS2 = genes[51:100],
#'     GS3 = genes[101:150]
#' )
#'
#' ## Run gene set hyperenrichment
#' K2res <- runGSEmods(K2res, 
#'                     genesets = genesetsMadeUp,
#'                     qthresh = 0.1)
#'                     
#' ## Run GSVA on genesets           
#' K2res <- runGSVAmods(K2res, 
#'                      ssGSEAalg = "gsva",
#'                      ssGSEAcores = 1,
#'                      verbose = FALSE)
#'
#' ## Aggregate paired gene sets
#' aggList <- list(c("GS12", "GS1", "GS2"))
#' K2res <- aggregateGSVAscores(aggList, K2res)
#'

aggregateGSVAscores <- function(aggList, K2res) {
    
    ## Run checks
    .isK2(K2res)
    
    ## Only perform this method for gsva runs
    if (K2meta(K2res)$ssGSEAalg != "gsva") {
        stop("Enrichment score aggregation only
            supported when ssGSEAalg == 'gsva'")
    }
    
    ## Get gSet
    gSet <- K2gSet(K2res)
    
    ## For each item in aggList subtract negative signature score from positive
    gAgg <- do.call(rbind, lapply(aggList, function(agg) {
        gSubUp <- exprs(gSet)[agg[2], ]
        gSubDown <- exprs(gSet)[agg[3], ]
        gSubAgg <- gSubUp - gSubDown
    }))
    rownames(gAgg) <- unlist(lapply(aggList, function(x) x[1]))
    
    ## Remove original scores
    origNames <- unlist(lapply(aggList, function(x) x[-1]))
    gSet <- gSet[!rownames(gSet) %in% origNames, ]
    
    ## Add new scores
    gExprs <- rbind(exprs(gSet), gAgg)
    gNew <- ExpressionSet(assayData = gExprs)
    pData(gNew) <- pData(gSet)
    
    ## Replace gSet
    K2gSet(K2res) <- gNew
    
    return(K2res)
    
}
