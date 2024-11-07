#' Take difference of two paired bi-directional gene set scores
#'
#' Replaces gene set results from paired up- and down- gene sets with the difference
#' of the up-regulated genes and down-regulated genes
#' @param K2res An object of class K2. The output of runDGEmods().
#' @param aggList A named list where each item is a character vector of length, 
#' 2, comprising the name of the 'up' gene set, and the name of the 'down'
#' gene set.
#' @return An object of class K2.
#' @references
#'    \insertRef{reed_2020}{K2Taxonomer}
#'    \insertRef{gsva}{K2Taxonomer}
#' @keywords clustering
#' @export
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
#'             stabThresh=0.5)
#'
#' ## Run differential analysis on each partition
#' K2res <- runDGEmods(K2res)
#'
#' ## Create dummy set of gene sets
#' DGEtable <- getDGETable(K2res)
#' genes <- unique(DGEtable$gene)
#' genesetsMadeUp <- list(
#'     GS1=genes[1:50],
#'     GS2=genes[51:100],
#'     GS3=genes[101:150])
#'
#' ## Run gene set hyperenrichment
#' K2res <- runGSEmods(K2res,
#'                 genesets=genesetsMadeUp,
#'                 qthresh=0.1)
#'
#' ## Run GSVA on genesets
#' K2res <- runGSVAmods(K2res,
#'                 ssGSEAalg='gsva',
#'                 ssGSEAcores=1,
#'                 verbose=FALSE)
#'
#' ## Aggregate paired gene sets
#' aggList <- list(c('GS12', 'GS1', 'GS2'))
#' K2res <- aggregateGeneSetscores(K2resaggList, K2res)
#'

aggregateGeneSetScores <- function(K2res, aggList) {

    ## Run checks
    .isK2(K2res)

    ## Get gMat
    gMat <- K2gMat(K2res)
    
    ## Only perform this method for gsva runs
    if (nrow(gMat) == 0) {
      stop("No gene set score matrix found. First, run runScoreGeneSets().\n")
    }
    
    ## Only perform this method for gsva runs
    if (!is.list(aggList) | is.null(names(aggList)) | 
        sum(names(aggList) == "") > 0) {
      stop("Argument, aggList, must be a list with names for every vector.\n")
    }

    ## For each item in aggList subtract negative signature score
    ## from positive
    gAgg <- do.call(rbind, lapply(aggList, function(agg) {
        gSubUp <- gMat[agg[1], ]
        gSubDown <- gMat[agg[2], ]
        gSubAgg <- gSubUp - gSubDown
        return(gSubAgg)
    }))
    rownames(gAgg) <- names(aggList)

    ## Remove original scores
    origNames <- unlist(lapply(aggList, function(x) x))
    gMat <- gMat[!rownames(gMat) %in% origNames, , drop = FALSE]

    ## Add new scores
    gNew <- rbind(gMat, gAgg)

    ## Replace gMat
    K2gMat(K2res) <- gNew
    
    ## Add geneset to K2genesets
    gsOld <- K2genesets(K2res)
    gsAdd <- lapply(aggList, function(agg) {
      return(unique(c(gsOld[[agg[1]]], gsOld[[agg[2]]])))
    })
    names(gsAdd) <- names(aggList)
    gsNew <- unlist(list(gsOld, gsAdd), recursive = FALSE)
    K2genesets(K2res) <- gsNew
    K2gene2Pathway(K2res) <- getGenePathways(gsNew)

    return(K2res)

}
